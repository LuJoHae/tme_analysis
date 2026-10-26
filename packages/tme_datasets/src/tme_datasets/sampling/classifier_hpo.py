"""Pure functional inner-loop optimization for compositional transforms, feature selection, and classifiers."""

from __future__ import annotations

import time
import warnings
from enum import Enum
from typing import Mapping, Sequence
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from joblib import Parallel, delayed  # type: ignore[import-untyped]
from sklearn.ensemble import (  # type: ignore[import-untyped]
    ExtraTreesClassifier,
    HistGradientBoostingClassifier,
    RandomForestClassifier,
)
from sklearn.exceptions import ConvergenceWarning  # type: ignore[import-untyped]
from sklearn.feature_selection import VarianceThreshold  # type: ignore[import-untyped]
from sklearn.linear_model import LogisticRegression  # type: ignore[import-untyped]
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score  # type: ignore[import-untyped]
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold  # type: ignore[import-untyped]
from sklearn.svm import SVC  # type: ignore[import-untyped]

from ..logging import get_logger

logger = get_logger("sampling.classifier_hpo")


class TransformType(str, Enum):
    """Compositional data transformation types for probability simplex deconvolution fractions."""

    RAW = "raw"
    CLR = "clr"
    LOGIT = "logit"
    ALR = "alr"
    CLR_RATIOS = "clr_ratios"


class FeatureSelectorType(str, Enum):
    """Feature selection strategy for cell state fractions."""

    NONE = "none"
    STABILITY_SELECTION = "stability_selection"
    ELASTICNET = "elasticnet"
    VARIANCE = "variance"


class ClassifierType(str, Enum):
    """Predictive classifier algorithm family."""

    LOGISTIC_REGRESSION = "logistic_regression"
    SVM_LINEAR = "svm_linear"
    SVM_RBF = "svm_rbf"
    RANDOM_FOREST = "random_forest"
    HIST_GRADIENT_BOOSTING = "hist_gradient_boosting"
    EXTRA_TREES = "extra_trees"


class InnerClassifierConfig(BaseModel):
    """Immutable parameter configuration for downstream feature selection and classifier training."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    transform: TransformType = TransformType.CLR
    selector: FeatureSelectorType = FeatureSelectorType.STABILITY_SELECTION
    classifier: ClassifierType = ClassifierType.LOGISTIC_REGRESSION
    stability_cutoff: float = 0.70
    stability_b: int = 40
    c_regularization: float = 1.0
    l1_ratio: float = 0.5
    rf_n_estimators: int = 100
    rf_max_depth: int = 3
    rf_min_samples_leaf: int = 2
    eps: float = 1e-5
    n_jobs: int = 1


class InnerOptimizationResult(BaseModel):
    """Immutable result of downstream feature selection and classifier optimization."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    best_config: InnerClassifierConfig
    best_loco_auc: float
    best_loco_pr_auc: float
    selected_features: tuple[str, ...]
    cohort_aucs: Mapping[str, float]
    all_results_df: pl.DataFrame
    elapsed_seconds: float


def apply_compositional_transform(
    X: np.ndarray,
    transform: TransformType,
    eps: float = 1e-5,
) -> np.ndarray:
    """Apply mathematically rigorous transformation to simplex-constrained fraction matrices."""
    X_safe: np.ndarray = np.asarray(np.clip(X, 0.0, 1.0), dtype=np.float64)

    match transform:
        case TransformType.RAW:
            return X_safe

        case TransformType.CLR:
            # Centered Log-Ratio: clr(x_k) = log(x_k + eps) - mean_k(log(x_k + eps))
            log_x = np.log(X_safe + eps)
            geom_mean = np.mean(log_x, axis=1, keepdims=True)
            res_clr: np.ndarray = np.asarray(log_x - geom_mean, dtype=np.float64)
            return res_clr

        case TransformType.LOGIT:
            # Proportion Logit: logit(x) = log((x + eps) / (1 - x + eps))
            numer = X_safe + eps
            denom = np.clip(1.0 - X_safe + eps, eps, None)
            res_logit: np.ndarray = np.asarray(np.log(numer / denom), dtype=np.float64)
            return res_logit

        case TransformType.ALR:
            # Additive Log-Ratio: alr(x_k) = log((x_k + eps) / (x_ref + eps))
            # Uses last column as reference
            ref_col = X_safe[:, -1:] + eps
            log_numer = np.log(X_safe[:, :-1] + eps)
            res_alr: np.ndarray = np.asarray(log_numer - np.log(ref_col), dtype=np.float64)
            return res_alr

        case TransformType.CLR_RATIOS:
            # Centered Log-Ratio augmented with all pairwise log-ratios log(s_i / s_j) = clr(s_i) - clr(s_j)
            log_x = np.log(X_safe + eps)
            geom_mean = np.mean(log_x, axis=1, keepdims=True)
            res_clr = np.asarray(log_x - geom_mean, dtype=np.float64)
            n_features = res_clr.shape[1]
            if n_features <= 1:
                return res_clr
            diff_cols = [
                res_clr[:, i : i + 1] - res_clr[:, j : j + 1]
                for i in range(n_features)
                for j in range(i + 1, n_features)
            ]
            res_augmented: np.ndarray = np.asarray(np.hstack([res_clr] + diff_cols), dtype=np.float64)
            return res_augmented


def get_transformed_feature_names(
    feature_names: Sequence[str],
    transform: TransformType,
) -> tuple[str, ...]:
    """Derive output feature names following compositional transformations."""
    match transform:
        case TransformType.RAW | TransformType.CLR | TransformType.LOGIT:
            return tuple(feature_names)
        case TransformType.ALR:
            return tuple(feature_names[:-1])
        case TransformType.CLR_RATIOS:
            n_feats = len(feature_names)
            ratio_names = [
                f"{feature_names[i]}_div_{feature_names[j]}"
                for i in range(n_feats)
                for j in range(i + 1, n_feats)
            ]
            return tuple(feature_names) + tuple(ratio_names)


def select_features_train_test(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    feature_names: Sequence[str],
    selector_type: FeatureSelectorType,
    config: InnerClassifierConfig,
) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    """Execute feature selection fitted strictly on training fold to prevent leakage."""
    n_features = X_train.shape[1]
    if n_features <= 2 or selector_type == FeatureSelectorType.NONE:
        return X_train, X_test, tuple(feature_names)

    match selector_type:
        case FeatureSelectorType.VARIANCE:
            vt = VarianceThreshold(threshold=1e-4)
            try:
                X_tr_sub = vt.fit_transform(X_train)
                X_te_sub = vt.transform(X_test)
                support = vt.get_support()
                selected = tuple(f for f, s in zip(feature_names, support) if s)
                if len(selected) >= 1:
                    return X_tr_sub, X_te_sub, selected
            except Exception:
                pass
            return X_train, X_test, tuple(feature_names)

        case FeatureSelectorType.ELASTICNET:
            # Sparse L1 feature selection via logistic regression
            try:
                lr = LogisticRegression(
                    solver="saga",
                    l1_ratio=config.l1_ratio,
                    C=config.c_regularization,
                    max_iter=300,
                    random_state=42,
                )
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=ConvergenceWarning)
                    lr.fit(X_train, y_train)
                coef_mag = np.abs(lr.coef_[0])
                support = coef_mag > 1e-4
                if np.sum(support) >= 2:
                    selected = tuple(f for f, s in zip(feature_names, support) if s)
                    return X_train[:, support], X_test[:, support], selected
            except Exception:
                pass
            return X_train, X_test, tuple(feature_names)

        case FeatureSelectorType.STABILITY_SELECTION:
            try:
                from selective_inference.stability_selection.sklearn_adapter import StabilitySelector  # type: ignore[import-untyped]

                selector = StabilitySelector(
                    cutoff=config.stability_cutoff,
                    B=config.stability_b,
                    random_state=42,
                    n_jobs=config.n_jobs,
                )
                selector.fit(X_train, y_train)
                selected_indices = selector.selected_variables_
                if len(selected_indices) >= 2:
                    sub_names = tuple(feature_names[i] for i in selected_indices)
                    return X_train[:, selected_indices], X_test[:, selected_indices], sub_names
            except Exception as exc:
                logger.debug("StabilitySelector fallback: %s", exc)

            return X_train, X_test, tuple(feature_names)


def fit_and_predict(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    classifier_type: ClassifierType,
    config: InnerClassifierConfig,
) -> np.ndarray:
    """Fit classifier on training fold and return continuous prediction probabilities on test fold."""
    match classifier_type:
        case ClassifierType.LOGISTIC_REGRESSION:
            clf = LogisticRegression(
                C=config.c_regularization,
                max_iter=200,
                solver="lbfgs",
                class_weight="balanced",
                random_state=42,
            )
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=ConvergenceWarning)
                clf.fit(X_train, y_train)
            return np.asarray(clf.predict_proba(X_test)[:, 1], dtype=np.float64)

        case ClassifierType.SVM_LINEAR:
            svm = SVC(
                kernel="linear",
                C=config.c_regularization,
                class_weight="balanced",
                random_state=42,
            )
            svm.fit(X_train, y_train)
            decision = svm.decision_function(X_test)
            return np.asarray(1.0 / (1.0 + np.exp(-decision)), dtype=np.float64)

        case ClassifierType.SVM_RBF:
            svm = SVC(
                kernel="rbf",
                C=config.c_regularization,
                class_weight="balanced",
                random_state=42,
            )
            svm.fit(X_train, y_train)
            decision = svm.decision_function(X_test)
            return np.asarray(1.0 / (1.0 + np.exp(-decision)), dtype=np.float64)

        case ClassifierType.RANDOM_FOREST:
            rf = RandomForestClassifier(
                n_estimators=config.rf_n_estimators,
                max_depth=config.rf_max_depth,
                min_samples_leaf=config.rf_min_samples_leaf,
                class_weight="balanced",
                random_state=42,
                n_jobs=config.n_jobs,
            )
            rf.fit(X_train, y_train)
            return np.asarray(rf.predict_proba(X_test)[:, 1], dtype=np.float64)

        case ClassifierType.HIST_GRADIENT_BOOSTING:
            hgb = HistGradientBoostingClassifier(
                max_iter=config.rf_n_estimators,
                max_depth=config.rf_max_depth,
                min_samples_leaf=max(2, config.rf_min_samples_leaf),
                l2_regularization=config.c_regularization,
                class_weight="balanced",
                random_state=42,
            )
            hgb.fit(X_train, y_train)
            return np.asarray(hgb.predict_proba(X_test)[:, 1], dtype=np.float64)

        case ClassifierType.EXTRA_TREES:
            et = ExtraTreesClassifier(
                n_estimators=config.rf_n_estimators,
                max_depth=config.rf_max_depth,
                min_samples_leaf=config.rf_min_samples_leaf,
                class_weight="balanced",
                random_state=42,
                n_jobs=config.n_jobs,
            )
            et.fit(X_train, y_train)
            return np.asarray(et.predict_proba(X_test)[:, 1], dtype=np.float64)


def evaluate_inner_pipeline(
    cohort_fractions: Mapping[str, np.ndarray],
    cohort_responses: Mapping[str, np.ndarray],
    feature_names: Sequence[str],
    config: InnerClassifierConfig,
) -> Result[tuple[float, float, dict[str, float], tuple[str, ...]], str]:
    """Execute Leave-One-Cohort-Out evaluation for a single feature selection + classifier configuration.

    Returns:
        Success(mean_auc, mean_pr_auc, cohort_aucs, final_selected_features) or Failure(error message).
    """
    cohort_ids = list(cohort_fractions.keys())
    if len(cohort_ids) < 2:
        return Failure("Need at least 2 bulk cohorts for Leave-One-Cohort-Out CV.")

    # 1. Apply compositional transform per cohort
    transformed_cohorts: dict[str, np.ndarray] = {
        cid: apply_compositional_transform(cohort_fractions[cid], config.transform, config.eps)
        for cid in cohort_ids
    }

    cohort_aucs: dict[str, float] = {}
    cohort_pr_aucs: dict[str, float] = {}
    last_selected_features: tuple[str, ...] = tuple(feature_names)

    for test_cid in cohort_ids:
        # Assemble training data
        train_X_list = [transformed_cohorts[cid] for cid in cohort_ids if cid != test_cid]
        train_y_list = [cohort_responses[cid] for cid in cohort_ids if cid != test_cid]

        train_X = np.vstack(train_X_list)
        train_y = np.concatenate(train_y_list)
        test_X = transformed_cohorts[test_cid]
        test_y = cohort_responses[test_cid]

        if len(np.unique(train_y)) < 2 or len(np.unique(test_y)) < 2:
            cohort_aucs[test_cid] = 0.5
            cohort_pr_aucs[test_cid] = float(np.mean(test_y))
            continue

        # Feature selection strictly on training fold
        tr_sub, te_sub, selected_cols = select_features_train_test(
            X_train=train_X,
            y_train=train_y,
            X_test=test_X,
            feature_names=feature_names,
            selector_type=config.selector,
            config=config,
        )
        last_selected_features = selected_cols

        try:
            pred_probs = fit_and_predict(
                X_train=tr_sub,
                y_train=train_y,
                X_test=te_sub,
                classifier_type=config.classifier,
                config=config,
            )
            score_auc = float(roc_auc_score(test_y, pred_probs))
            prec, rec, _ = precision_recall_curve(test_y, pred_probs)
            score_pr = float(auc(rec, prec))
            cohort_aucs[test_cid] = score_auc
            cohort_pr_aucs[test_cid] = score_pr
        except Exception:
            cohort_aucs[test_cid] = 0.5
            cohort_pr_aucs[test_cid] = float(np.mean(test_y))

    mean_auc = float(np.mean(list(cohort_aucs.values())))
    mean_pr = float(np.mean(list(cohort_pr_aucs.values())))

    return Success((mean_auc, mean_pr, cohort_aucs, last_selected_features))


def generate_candidate_inner_configs(n_jobs: int = 1) -> tuple[InnerClassifierConfig, ...]:
    """Generate curated candidate configurations spanning transforms, selectors, and classifiers."""
    configs: list[InnerClassifierConfig] = [
        # Baseline: CLR + Logistic Regression
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
            c_regularization=1.0,
            n_jobs=n_jobs,
        ),
        # CLR + Relaxed L2 Logistic Regression (C=3.0)
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
            c_regularization=3.0,
            n_jobs=n_jobs,
        ),
        # CLR + Stability Selection + Logistic Regression
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.STABILITY_SELECTION,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
            stability_cutoff=0.60,
            c_regularization=2.0,
            n_jobs=n_jobs,
        ),
        # Logit + ElasticNet Selection + Logistic Regression
        InnerClassifierConfig(
            transform=TransformType.LOGIT,
            selector=FeatureSelectorType.ELASTICNET,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
            c_regularization=1.0,
            l1_ratio=0.3,
            n_jobs=n_jobs,
        ),
        # CLR + Random Forest (Standard depth 3)
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.RANDOM_FOREST,
            rf_n_estimators=100,
            rf_max_depth=3,
            rf_min_samples_leaf=2,
            n_jobs=n_jobs,
        ),
        # CLR + Random Forest (Deeper depth 4, min_leaf=1)
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.RANDOM_FOREST,
            rf_n_estimators=100,
            rf_max_depth=4,
            rf_min_samples_leaf=1,
            n_jobs=n_jobs,
        ),
        # CLR + HistGradientBoosting (LightGBM-style GBDT, depth 3)
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.HIST_GRADIENT_BOOSTING,
            rf_n_estimators=100,
            rf_max_depth=3,
            rf_min_samples_leaf=5,
            c_regularization=0.5,
            n_jobs=n_jobs,
        ),
        # CLR + HistGradientBoosting (Deeper depth 4, min_leaf=2, weak reg)
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.HIST_GRADIENT_BOOSTING,
            rf_n_estimators=100,
            rf_max_depth=4,
            rf_min_samples_leaf=2,
            c_regularization=0.1,
            n_jobs=n_jobs,
        ),
        # CLR + ExtraTrees (Depth 3, min_leaf=2)
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.EXTRA_TREES,
            rf_n_estimators=100,
            rf_max_depth=3,
            rf_min_samples_leaf=2,
            n_jobs=n_jobs,
        ),
        # CLR + ExtraTrees (Deeper depth 4, min_leaf=1)
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.EXTRA_TREES,
            rf_n_estimators=100,
            rf_max_depth=4,
            rf_min_samples_leaf=1,
            n_jobs=n_jobs,
        ),
        # CLR + ExtraTrees (Deeper depth 5, min_leaf=1)
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.EXTRA_TREES,
            rf_n_estimators=100,
            rf_max_depth=5,
            rf_min_samples_leaf=1,
            n_jobs=n_jobs,
        ),
        # CLR Pairwise Ratios + ElasticNet + Logistic Regression (Direct Immune Balance Contrasts)
        InnerClassifierConfig(
            transform=TransformType.CLR_RATIOS,
            selector=FeatureSelectorType.ELASTICNET,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
            c_regularization=2.0,
            l1_ratio=0.3,
            n_jobs=n_jobs,
        ),
        # CLR Pairwise Ratios + ExtraTrees (Depth 4, min_leaf=1)
        InnerClassifierConfig(
            transform=TransformType.CLR_RATIOS,
            selector=FeatureSelectorType.ELASTICNET,
            classifier=ClassifierType.EXTRA_TREES,
            rf_n_estimators=100,
            rf_max_depth=4,
            rf_min_samples_leaf=1,
            n_jobs=n_jobs,
        ),
        # CLR Pairwise Ratios + HistGradientBoosting
        InnerClassifierConfig(
            transform=TransformType.CLR_RATIOS,
            selector=FeatureSelectorType.ELASTICNET,
            classifier=ClassifierType.HIST_GRADIENT_BOOSTING,
            rf_n_estimators=100,
            rf_max_depth=4,
            rf_min_samples_leaf=2,
            c_regularization=0.1,
            n_jobs=n_jobs,
        ),
    ]
    return tuple(configs)


def run_inner_classifier_hpo(
    cohort_fractions: Mapping[str, np.ndarray],
    cohort_responses: Mapping[str, np.ndarray],
    feature_names: Sequence[str],
    candidate_configs: Sequence[InnerClassifierConfig] | None = None,
    n_jobs: int = 1,
) -> Result[InnerOptimizationResult, str]:
    """Execute high-speed inner loop search across transforms, feature selectors, and classifiers."""
    t0 = time.time()
    configs = candidate_configs if candidate_configs is not None else generate_candidate_inner_configs(n_jobs=n_jobs)

    results: list[dict[str, object]] = []
    best_auc = -1.0
    best_cfg = configs[0]
    best_pr = 0.0
    best_cohort_aucs: dict[str, float] = {}
    best_features: tuple[str, ...] = tuple(feature_names)

    for cfg in configs:
        res = evaluate_inner_pipeline(
            cohort_fractions=cohort_fractions,
            cohort_responses=cohort_responses,
            feature_names=feature_names,
            config=cfg,
        )
        match res:
            case Failure(err):
                logger.warning("Inner evaluation failed for config %s: %s", cfg, err)
                continue
            case Success((mean_auc, mean_pr, c_aucs, sel_feats)):
                results.append({
                    "transform": cfg.transform.value,
                    "selector": cfg.selector.value,
                    "classifier": cfg.classifier.value,
                    "c_regularization": cfg.c_regularization,
                    "mean_loco_auc": mean_auc,
                    "mean_loco_pr_auc": mean_pr,
                    "n_selected_features": len(sel_feats),
                })
                if mean_auc > best_auc:
                    best_auc = mean_auc
                    best_pr = mean_pr
                    best_cfg = cfg
                    best_cohort_aucs = c_aucs
                    best_features = sel_feats

    if not results:
        return Failure("All candidate classifier configurations failed to evaluate.")

    res_df = pl.DataFrame(results)
    elapsed = round(time.time() - t0, 3)

    return Success(
        InnerOptimizationResult(
            best_config=best_cfg,
            best_loco_auc=best_auc,
            best_loco_pr_auc=best_pr,
            selected_features=best_features,
            cohort_aucs=best_cohort_aucs,
            all_results_df=res_df,
            elapsed_seconds=elapsed,
        )
    )


def _evaluate_single_fold(
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    X_trans: np.ndarray,
    y: np.ndarray,
    feature_names: Sequence[str],
    config: InnerClassifierConfig,
) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    """Evaluate feature selection and classifier on a single train-test split.

    Returns:
        (test_idx, pred_probs, selected_features)
    """
    X_tr, y_tr = X_trans[train_idx], y[train_idx]
    X_te, y_te = X_trans[test_idx], y[test_idx]

    prior_prob = float(np.mean(y_tr)) if len(y_tr) > 0 else 0.5

    if len(np.unique(y_tr)) < 2 or len(np.unique(y_te)) < 2:
        return test_idx, np.full(len(test_idx), prior_prob, dtype=np.float64), tuple(feature_names)

    # 1. Feature selection strictly on training fold
    tr_sub, te_sub, selected_cols = select_features_train_test(
        X_train=X_tr,
        y_train=y_tr,
        X_test=X_te,
        feature_names=feature_names,
        selector_type=config.selector,
        config=config,
    )

    # 2. Fit and predict
    try:
        pred_probs = fit_and_predict(
            X_train=tr_sub,
            y_train=y_tr,
            X_test=te_sub,
            classifier_type=config.classifier,
            config=config,
        )
        return test_idx, pred_probs, selected_cols
    except Exception as exc:
        logger.debug("Fold prediction failed (%s); falling back to prior.", exc)
        return test_idx, np.full(len(test_idx), prior_prob, dtype=np.float64), tuple(feature_names)


def evaluate_single_cohort_stratified_cv(
    frac_mat: np.ndarray,
    y: np.ndarray,
    feature_names: Sequence[str],
    config: InnerClassifierConfig,
    n_splits: int = 5,
    n_repeats: int = 3,
    seed: int = 42,
) -> Result[tuple[float, float, tuple[str, ...]], str]:
    """Execute leakage-free Repeated Stratified K-Fold CV on a single cohort.

    Aggregates out-of-fold predicted probabilities across the complete cohort to eliminate
    discretization artifacts in small folds, computing pooled ROC-AUC and PR-AUC.
    Parallelizes across all fold evaluations concurrently utilizing available CPU workers.
    """
    if len(np.unique(y)) < 2:
        return Failure("Target vector y must contain at least 2 distinct classes.")

    min_class_count = int(np.min(np.bincount(y)))
    actual_n_splits = max(2, min(n_splits, min_class_count))
    if len(y) < actual_n_splits * 2:
        return Failure(f"Sample size {len(y)} too small for {actual_n_splits}-fold CV.")

    X_trans = apply_compositional_transform(frac_mat, config.transform, config.eps)
    active_feature_names = get_transformed_feature_names(feature_names, config.transform)

    if n_repeats > 1:
        cv = RepeatedStratifiedKFold(n_splits=actual_n_splits, n_repeats=n_repeats, random_state=seed)
    else:
        cv = StratifiedKFold(n_splits=actual_n_splits, shuffle=True, random_state=seed)

    splits = list(cv.split(X_trans, y))
    n_workers = min(max(1, config.n_jobs), len(splits))

    fold_results = Parallel(n_jobs=n_workers, prefer="threads")(
        delayed(_evaluate_single_fold)(
            train_idx=train_idx,
            test_idx=test_idx,
            X_trans=X_trans,
            y=y,
            feature_names=active_feature_names,
            config=config,
        )
        for train_idx, test_idx in splits
    )

    all_selected = [feat for res in fold_results for feat in res[2]]

    # Feature consensus: features selected in >= 30% of folds
    cutoff = len(splits) * 0.30
    from collections import Counter
    counts = Counter(all_selected)
    consensus_features = tuple(f for f in active_feature_names if counts[f] >= cutoff)
    if not consensus_features:
        consensus_features = tuple(active_feature_names)

    # Compute Pooled Out-of-Fold (OOF) AUC per repeat
    repeat_aucs: list[float] = []
    repeat_prs: list[float] = []

    for r in range(n_repeats):
        repeat_splits = fold_results[r * actual_n_splits : (r + 1) * actual_n_splits]
        oof_probs = np.full(len(y), 0.5, dtype=np.float64)
        for test_idx, pred_probs, _ in repeat_splits:
            oof_probs[test_idx] = pred_probs

        try:
            r_auc = float(roc_auc_score(y, oof_probs))
            prec, rec, _ = precision_recall_curve(y, oof_probs)
            r_pr = float(auc(rec, prec))
        except Exception:
            r_auc = 0.5
            r_pr = float(np.mean(y))

        repeat_aucs.append(r_auc)
        repeat_prs.append(r_pr)

    mean_auc = float(np.mean(repeat_aucs))
    mean_pr = float(np.mean(repeat_prs))

    return Success((mean_auc, mean_pr, consensus_features))


def compute_single_cohort_oof_predictions(
    frac_mat: np.ndarray,
    y: np.ndarray,
    feature_names: Sequence[str],
    config: InnerClassifierConfig,
    n_splits: int = 5,
    n_repeats: int = 3,
    seed: int = 42,
) -> Result[np.ndarray, str]:
    """Compute pooled out-of-fold predicted probabilities across Repeated Stratified CV splits."""
    if len(np.unique(y)) < 2:
        return Failure("Target vector y must contain at least 2 distinct classes.")

    min_class_count = int(np.min(np.bincount(y)))
    actual_n_splits = max(2, min(n_splits, min_class_count))
    if len(y) < actual_n_splits * 2:
        return Failure(f"Sample size {len(y)} too small for {actual_n_splits}-fold CV.")

    X_trans = apply_compositional_transform(frac_mat, config.transform, config.eps)
    active_feature_names = get_transformed_feature_names(feature_names, config.transform)

    if n_repeats > 1:
        cv = RepeatedStratifiedKFold(n_splits=actual_n_splits, n_repeats=n_repeats, random_state=seed)
    else:
        cv = StratifiedKFold(n_splits=actual_n_splits, shuffle=True, random_state=seed)

    splits = list(cv.split(X_trans, y))
    n_workers = min(max(1, config.n_jobs), len(splits))

    fold_results = Parallel(n_jobs=n_workers, prefer="threads")(
        delayed(_evaluate_single_fold)(
            train_idx=train_idx,
            test_idx=test_idx,
            X_trans=X_trans,
            y=y,
            feature_names=active_feature_names,
            config=config,
        )
        for train_idx, test_idx in splits
    )

    all_repeat_probs = np.zeros(len(y), dtype=np.float64)
    for r in range(n_repeats):
        repeat_splits = fold_results[r * actual_n_splits : (r + 1) * actual_n_splits]
        oof_probs = np.full(len(y), 0.5, dtype=np.float64)
        for test_idx, pred_probs, _ in repeat_splits:
            oof_probs[test_idx] = pred_probs
        all_repeat_probs += oof_probs

    mean_probs = all_repeat_probs / float(n_repeats)
    return Success(mean_probs)


def run_single_cohort_classifier_hpo(
    frac_mat: np.ndarray,
    y: np.ndarray,
    feature_names: Sequence[str],
    candidate_configs: Sequence[InnerClassifierConfig] | None = None,
    n_splits: int = 5,
    n_repeats: int = 3,
    seed: int = 42,
    n_jobs: int = 16,
) -> Result[InnerOptimizationResult, str]:
    """Execute high-speed inner loop search across transforms, feature selectors, and classifiers for a single cohort."""
    t0 = time.time()
    configs = candidate_configs if candidate_configs is not None else generate_candidate_inner_configs(n_jobs=n_jobs)

    results: list[dict[str, object]] = []
    best_auc = -1.0
    best_cfg = configs[0]
    best_pr = 0.0
    best_features: tuple[str, ...] = tuple(feature_names)

    for cfg in configs:
        res = evaluate_single_cohort_stratified_cv(
            frac_mat=frac_mat,
            y=y,
            feature_names=feature_names,
            config=cfg,
            n_splits=n_splits,
            n_repeats=n_repeats,
            seed=seed,
        )
        match res:
            case Failure(err):
                logger.warning("Single-cohort evaluation failed for config %s: %s", cfg, err)
                continue
            case Success((mean_auc, mean_pr, sel_feats)):
                results.append({
                    "transform": cfg.transform.value,
                    "selector": cfg.selector.value,
                    "classifier": cfg.classifier.value,
                    "c_regularization": cfg.c_regularization,
                    "mean_loco_auc": mean_auc,
                    "mean_loco_pr_auc": mean_pr,
                    "n_selected_features": len(sel_feats),
                })
                if mean_auc > best_auc:
                    best_auc = mean_auc
                    best_pr = mean_pr
                    best_cfg = cfg
                    best_features = sel_feats

    if not results:
        return Failure("All candidate classifier configurations failed to evaluate.")

    res_df = pl.DataFrame(results)
    elapsed = round(time.time() - t0, 3)

    return Success(
        InnerOptimizationResult(
            best_config=best_cfg,
            best_loco_auc=best_auc,
            best_loco_pr_auc=best_pr,
            selected_features=best_features,
            cohort_aucs={"target": best_auc},
            all_results_df=res_df,
            elapsed_seconds=elapsed,
        )
    )

