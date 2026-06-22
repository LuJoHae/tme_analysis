from sklearn.base import BaseEstimator, TransformerMixin
from dataclasses import dataclass
from typing import Any
import seaborn as sns
from numpy import dtype, generic, ndarray
from pandas import DataFrame, Series
import numpy as np
import pandas as pd
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.model_selection import cross_val_predict
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import precision_recall_curve, average_precision_score, roc_auc_score, roc_curve
from sklearn.model_selection import LeaveOneOut
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


@dataclass
class BootstrapCIResult:
    """Container for bootstrap confidence interval results from ROC and PR curves."""
    base_fpr: np.ndarray
    base_recall: np.ndarray
    pr_ci: np.ndarray
    precisions_lower: np.ndarray
    precisions_lower_small: np.ndarray
    precisions_upper: np.ndarray
    precisions_upper_small: np.ndarray
    roc_ci: np.ndarray
    tprs_lower: np.ndarray
    tprs_lower_small: np.ndarray
    tprs_upper: np.ndarray
    tprs_upper_small: np.ndarray


class AppendCategoricalAfterPCA(BaseEstimator, TransformerMixin):
    """Custom transformer that performs PCA on the first `n_pca_features` columns,
    then concatenates the remaining columns (one-hot encoded categorical metadata)
    unchanged."""

    def __init__(self, n_components, n_pca_features):
        self.n_components = n_components
        self.n_pca_features = n_pca_features

    def fit(self, X, y=None):
        self.pca_ = PCA(n_components=self.n_components)
        self.pca_.fit(X[:, :self.n_pca_features])
        return self

    def transform(self, X):
        x_pca = self.pca_.transform(X[:, :self.n_pca_features])
        x_extra = X[:, self.n_pca_features:]
        return np.hstack([x_pca, x_extra])


def rf_pipeline(
        x, y,
        axes,
        n_estimators=50,
        n_components=2,
        categorical_metadata=None,
        max_depth: int | None = None
):
    clf, x_input = construct_features_and_pipeline(
        x, categorical_metadata, n_components, n_estimators, max_depth=max_depth
    )

    classes_ = np.unique(y)
    loo = LeaveOneOut()
    y_proba = cross_val_predict(clf, x_input, y, cv=loo, method="predict_proba")[
        :, list(classes_).index("R")]

    # Cast to numpy array to prevent Pandas index alignment issues during bootstrapping
    y_true = np.array((y == 'R').astype(int))

    ci = bootstrap_ci_intervals(y_proba, y_true)

    plot_decision_regions_pca(clf, x_input, y, axes[2], n_components,
                              categorical_metadata,
                              x.shape[1] if categorical_metadata is not None else
                              x.shape[1])

    plot_rf_pipeline_results(axes, ci, y_proba, y_true)

    return y_proba, axes, clf


def bootstrap_ci_intervals(y_proba,
                           y_true: ndarray[tuple[Any, ...], dtype[generic[Any]]] |
                                   ndarray[tuple[Any, ...], dtype[
                                       Any]]) -> BootstrapCIResult:
    # --- Bootstrap Procedure ---
    n_bootstraps = 1000
    rng = np.random.RandomState(42)

    bootstrapped_roc_aucs, bootstrapped_pr_aucs = [], []
    tprs_boot, precisions_boot = [], []

    # Define fixed domains for interpolation
    base_fpr = np.linspace(0, 1, 101)
    base_recall = np.linspace(0, 1, 101)

    for _ in range(n_bootstraps):
        # Draw sample with replacement
        indices = rng.randint(0, len(y_proba), len(y_proba))

        # Skip samples that do not contain both classes (undefined AUC)
        if len(np.unique(y_true[indices])) < 2:
            continue

        y_true_b, y_proba_b = y_true[indices], y_proba[indices]

        # Scalar Metrics
        bootstrapped_roc_aucs.append(roc_auc_score(y_true_b, y_proba_b))
        bootstrapped_pr_aucs.append(average_precision_score(y_true_b, y_proba_b))

        # Curve Interpolation
        fpr_b, tpr_b, _ = roc_curve(y_true_b, y_proba_b)
        tpr_interp = np.interp(base_fpr, fpr_b, tpr_b)
        tpr_interp[0] = 0.0  # Boundary condition
        tprs_boot.append(tpr_interp)

        precision_b, recall_b, _ = precision_recall_curve(y_true_b, y_proba_b)
        # PR curve interpolates backwards (recall is monotonically decreasing in sklearn output)
        precision_interp = np.interp(base_recall, recall_b[::-1], precision_b[::-1])
        precisions_boot.append(precision_interp)

    # Compute 95% Confidence Intervals for scalars
    interval = [2.5, 97.5]
    roc_ci = np.percentile(bootstrapped_roc_aucs, interval)
    pr_ci = np.percentile(bootstrapped_pr_aucs, interval)
    interval_small = [15.87, 84.13]
    roc_ci_small = np.percentile(bootstrapped_roc_aucs, interval_small)
    pr_ci_small = np.percentile(bootstrapped_pr_aucs, interval_small)

    # Compute 95% Confidence Intervals for pointwise curves
    tprs_lower, tprs_upper = np.percentile(tprs_boot, interval, axis=0)
    precisions_lower, precisions_upper = np.percentile(precisions_boot, interval,
                                                       axis=0)
    tprs_lower_small, tprs_upper_small = np.percentile(tprs_boot, interval_small,
                                                       axis=0)
    precisions_lower_small, precisions_upper_small = np.percentile(precisions_boot,
                                                                   interval_small,
                                                                   axis=0)
    return BootstrapCIResult(
        base_fpr=base_fpr,
        base_recall=base_recall,
        pr_ci=pr_ci,
        precisions_lower=precisions_lower,
        precisions_lower_small=precisions_lower_small,
        precisions_upper=precisions_upper,
        precisions_upper_small=precisions_upper_small,
        roc_ci=roc_ci,
        tprs_lower=tprs_lower,
        tprs_lower_small=tprs_lower_small,
        tprs_upper=tprs_upper,
        tprs_upper_small=tprs_upper_small,
    )


def plot_rf_pipeline_results(axes, ci: BootstrapCIResult, y_proba,
                             y_true: ndarray[tuple[Any, ...], dtype[generic[Any]]] |
                                     ndarray[tuple[Any, ...], dtype[Any]]):
    # --- Visualization ---

    # Main ROC point estimate
    fpr, tpr, _ = roc_curve(y_true, y_proba)
    auc = roc_auc_score(y_true, y_proba)
    axes[0].plot(fpr, tpr,
                 label=f"ROC (AUC = {auc:.3f})\n95% CI: [{ci.roc_ci[0]:.3f}, {ci.roc_ci[1]:.3f}]")
    axes[0].fill_between(ci.base_fpr, ci.tprs_lower, ci.tprs_upper, color='blue',
                         alpha=0.1, label='95% CI Region')
    axes[0].fill_between(ci.base_fpr, ci.tprs_lower_small, ci.tprs_upper_small,
                         color='blue', alpha=0.1, label='68% CI Region')
    axes[0].plot([0, 1], [0, 1], linestyle="--", color="grey")
    axes[0].set_xlabel("False Positive Rate")
    axes[0].set_ylabel("True Positive Rate")
    axes[0].set_title("ROC Curve (R vs NR)")
    axes[0].legend(loc="lower right")

    # Main PR point estimate
    precision, recall, _ = precision_recall_curve(y_true, y_proba)
    pr_auc = average_precision_score(y_true, y_proba)
    baseline = y_true.mean()
    axes[1].plot(recall, precision,
                 label=f"PR (AUC = {pr_auc:.3f})\n95% CI: [{ci.pr_ci[0]:.3f}, {ci.pr_ci[1]:.3f}]")
    axes[1].fill_between(ci.base_recall, ci.precisions_lower, ci.precisions_upper,
                         color='blue', alpha=0.1, label='95% CI Region')
    axes[1].fill_between(ci.base_recall, ci.precisions_lower_small,
                         ci.precisions_upper_small, color='blue', alpha=0.1,
                         label='68% CI Region')
    axes[1].axhline(baseline, linestyle="--", color="grey",
                    label=f"Baseline = {baseline:.3f}")
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    axes[1].set_title("Precision-Recall Curve (R vs NR)")
    axes[1].legend(loc="lower left")


def construct_features_and_pipeline(
        x,
        categorical_metadata,
        n_components: int,
        n_estimators: int,
        max_depth: int | None = 0
) -> tuple[
    Pipeline, Series[Any] | DataFrame | Any
]:
    # One-hot encode categorical metadata and append AFTER PCA via the custom transformer
    if categorical_metadata is not None:
        if isinstance(categorical_metadata, pd.Series):
            categorical_metadata = categorical_metadata.to_frame()
        cat_dummies = pd.get_dummies(categorical_metadata).astype(float)
        cat_dummies = cat_dummies.loc[x.index]
        n_pca_features = x.shape[1]
        x_full = pd.concat([x, cat_dummies], axis=1)

        # StandardScaler is applied only to the numeric (pre-PCA) features;
        # one-hot columns are passed through unchanged.
        from sklearn.compose import ColumnTransformer
        preprocessor = ColumnTransformer(
            transformers=[
                ("scale", StandardScaler(), list(range(n_pca_features))),
            ],
            remainder="passthrough",
        )

        clf = make_pipeline(
            preprocessor,
            AppendCategoricalAfterPCA(n_components=n_components,
                                      n_pca_features=n_pca_features),
            RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth,
                                   class_weight="balanced", random_state=42)
        )
        x_input = x_full
    else:
        clf = make_pipeline(
            StandardScaler(),
            PCA(n_components=n_components),
            RandomForestClassifier(n_estimators=n_estimators, class_weight="balanced",
                                   random_state=42)
        )
        x_input = x
    return clf, x_input


def plot_decision_regions_pca(clf, x_input, y, ax, n_components, categorical_metadata,
                              n_pca_features):
    """Plot decision regions in PCA space (first 2 components only)."""
    # Fit the classifier on all data
    clf.fit(x_input, y)

    # Transform data to PCA space for visualization
    if categorical_metadata is not None:
        # Extract the PCA transformation from the pipeline
        preprocessor = clf.named_steps['columntransformer']
        x_scaled = preprocessor.transform(x_input)
        pca_transformer = clf.named_steps['appendcategoricalafterpca']
        x_pca_full = pca_transformer.transform(x_scaled)
        x_pca = x_pca_full[:, :n_components]
    else:
        scaler = clf.named_steps['standardscaler']
        pca = clf.named_steps['pca']
        x_scaled = scaler.transform(x_input)
        x_pca = pca.transform(x_scaled)

    # Create mesh grid
    h = 0.02
    x_min, x_max = x_pca[:, 0].min() - 1, x_pca[:, 0].max() + 1
    y_min, y_max = x_pca[:, 1].min() - 1, x_pca[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

    # For prediction, we need to pad with zeros for components > 2 and categorical features
    mesh_points = np.c_[xx.ravel(), yy.ravel()]
    if categorical_metadata is not None:
        # Pad with zeros for additional PCA components and categorical features
        n_extra_features = x_pca_full.shape[1] - n_components
        mesh_points = np.hstack(
            [mesh_points, np.zeros((mesh_points.shape[0], n_extra_features))])
    elif n_components > 2:
        # Pad with zeros for additional PCA components
        mesh_points = np.hstack(
            [mesh_points, np.zeros((mesh_points.shape[0], n_components - 2))])

    # Predict on mesh
    Z = clf.named_steps['randomforestclassifier'].predict(mesh_points)
    Z_encoded = (Z == 'R').astype(int)
    Z_encoded = Z_encoded.reshape(xx.shape)

    # Plot decision boundary
    cmap_light = ListedColormap(['#FFAAAA', '#AAAAFF'])
    cmap_bold = ListedColormap(['#FF0000', '#0000FF'])
    ax.contourf(xx, yy, Z_encoded, alpha=0.3, cmap=cmap_light)

    # Plot data points
    y_encoded = (y == 'R').astype(int)
    scatter = ax.scatter(x_pca[:, 0], x_pca[:, 1], c=y_encoded, cmap=cmap_bold,
                         edgecolor='k', s=50, alpha=0.7)

    ax.set_xlabel('First Principal Component')
    ax.set_ylabel('Second Principal Component')
    ax.set_title('Decision Regions in PCA Space')
    ax.legend(handles=scatter.legend_elements()[0], labels=['NR', 'R'], loc='best')


def plot_subtype_responser_fractions(data, axes, title):
    sns.countplot(
        x="MFP",
        hue="Response",
        data=data,
        ax=axes[0]
    )
    axes[0].set_title(title + ": Responses")

    sns.histplot(
        data=data,
        x='MFP',
        hue='Response',
        multiple='fill',
        shrink=0.8,
        ax=axes[1]
    )
    axes[1].set_title(title + ": Response Fractions")