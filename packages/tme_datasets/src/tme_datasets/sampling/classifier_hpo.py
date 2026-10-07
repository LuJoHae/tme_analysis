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
from sklearn.ensemble import RandomForestClassifier  # type: ignore[import-untyped]
from sklearn.exceptions import ConvergenceWarning  # type: ignore[import-untyped]
from sklearn.feature_selection import VarianceThreshold  # type: ignore[import-untyped]
from sklearn.linear_model import LogisticRegression  # type: ignore[import-untyped]
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score  # type: ignore[import-untyped]
from sklearn.svm import SVC  # type: ignore[import-untyped]

from ..logging import get_logger

logger = get_logger("sampling.classifier_hpo")


class TransformType(str, Enum):
    """Compositional data transformation types for probability simplex deconvolution fractions."""

    RAW = "raw"
    CLR = "clr"
    LOGIT = "logit"
    ALR = "alr"


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
                min_samples_leaf=2,
                class_weight="balanced",
                random_state=42,
                n_jobs=config.n_jobs,
            )
            rf.fit(X_train, y_train)
            return np.asarray(rf.predict_proba(X_test)[:, 1], dtype=np.float64)


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
        # CLR + Stability Selection + Logistic Regression
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.STABILITY_SELECTION,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
            stability_cutoff=0.65,
            c_regularization=1.0,
            n_jobs=n_jobs,
        ),
        # Logit + ElasticNet Selection + Logistic Regression
        InnerClassifierConfig(
            transform=TransformType.LOGIT,
            selector=FeatureSelectorType.ELASTICNET,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
            c_regularization=0.5,
            l1_ratio=0.5,
            n_jobs=n_jobs,
        ),
        # CLR + Linear SVM
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.SVM_LINEAR,
            c_regularization=1.0,
            n_jobs=n_jobs,
        ),
        # Raw + Random Forest
        InnerClassifierConfig(
            transform=TransformType.RAW,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.RANDOM_FOREST,
            rf_n_estimators=100,
            rf_max_depth=3,
            n_jobs=n_jobs,
        ),
        # CLR + Random Forest
        InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.RANDOM_FOREST,
            rf_n_estimators=100,
            rf_max_depth=3,
            n_jobs=n_jobs,
        ),
        # Logit + Stability Selection + Logistic Regression
        InnerClassifierConfig(
            transform=TransformType.LOGIT,
            selector=FeatureSelectorType.STABILITY_SELECTION,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
            stability_cutoff=0.70,
            c_regularization=1.0,
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
