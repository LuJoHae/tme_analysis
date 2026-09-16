#!/usr/bin/env python3
"""
Calculate unweighted, sample-weighted, and uncertainty-weighted AUC metrics
for Bagaev-constrained deconvolution vs direct expression models.

Adheres strictly to the Functional Python Coding Style Rules:
- Pure functions and immutable frozen data models
- Returns Result monad for error handling
- Polars for declarative dataframe transformations
- Exhaustive pattern matching and type annotations
"""

from enum import Enum
from pathlib import Path
from typing import Final, Sequence, assert_never
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success


class WeightingMethod(str, Enum):
    UNWEIGHTED = "Unweighted"
    SAMPLE_SIZE = "Sample Size (N)"
    INVERSE_VARIANCE = "Inverse Variance (1/sigma^2)"
    INVERSE_STD = "Inverse Std (1/sigma)"


class CohortScope(str, Enum):
    ALL_COHORTS = "All 9 Cohorts (Pan-Cancer)"
    MELANOMA = "Melanoma (4 Cohorts)"


class PredictionMode(str, Enum):
    DECONVOLUTED_FRACTIONS = "Deconvoluted Fractions"
    DIRECT_EXPRESSION = "Direct Expression"


class CohortMeta(BaseModel):
    model_config = ConfigDict(frozen=True)
    cohort: str
    cancer_type: str
    sample_size: int
    is_melanoma: bool


class MetricSummaryRow(BaseModel):
    model_config = ConfigDict(frozen=True)
    cohort: str
    mode: str
    metric: str
    mean_val: float
    std_val: float
    sample_size: int
    is_melanoma: bool


class WeightedResult(BaseModel):
    model_config = ConfigDict(frozen=True)
    scope: CohortScope
    mode: PredictionMode
    metric: str
    weighting: WeightingMethod
    weighted_mean: float
    total_samples: int
    n_cohorts: int


class BenchmarkComparison(BaseModel):
    model_config = ConfigDict(frozen=True)
    scope: str
    metric: str
    weighting: str
    deconv_auc: float
    direct_auc: float
    delta_deconv_minus_direct: float


# Registry of known iAtlas cohorts with sample counts matching experimental cohorts
COHORT_REGISTRY: Final[tuple[CohortMeta, ...]] = (
    CohortMeta(cohort="Hugo-iAtlas", cancer_type="Melanoma", sample_size=27, is_melanoma=True),
    CohortMeta(cohort="Riaz-iAtlas", cancer_type="Melanoma", sample_size=98, is_melanoma=True),
    CohortMeta(cohort="Liu-iAtlas", cancer_type="Melanoma", sample_size=122, is_melanoma=True),
    CohortMeta(cohort="Gide-iAtlas", cancer_type="Melanoma", sample_size=91, is_melanoma=True),
    CohortMeta(cohort="Rosenberg-iAtlas", cancer_type="Bladder (BLCA)", sample_size=298, is_melanoma=False),
    CohortMeta(cohort="Padron-iAtlas", cancer_type="Pancreatic (PAAD)", sample_size=85, is_melanoma=False),
    CohortMeta(cohort="Anders-iAtlas", cancer_type="Breast (BRCA)", sample_size=31, is_melanoma=False),
    CohortMeta(cohort="McDermott-iAtlas", cancer_type="Renal (RCC)", sample_size=247, is_melanoma=False),
    CohortMeta(cohort="Choueiri-iAtlas", cancer_type="Renal (ccRCC)", sample_size=16, is_melanoma=False),
)

EPSILON: Final[float] = 1e-8


# Pure Core Functions


def find_cohort_meta(cohort_name: str) -> Maybe[CohortMeta]:
    """Finds cohort metadata from registry declaratively."""
    matches = tuple(c for c in COHORT_REGISTRY if c.cohort == cohort_name)
    return Some(matches[0]) if len(matches) > 0 else Nothing


def load_prediction_dataset(filepath: Path) -> Result[pl.DataFrame, str]:
    """Pure IO wrapper loading and validating the prediction results dataset."""
    return (
        Success(filepath)
        .bind(
            lambda p: Success(pl.read_csv(p))
            if p.is_file()
            else Failure(f"File not found: {p}")
        )
        .bind(
            lambda df: Success(df)
            if {"Cohort", "Mode", "ROC_AUC", "PR_AUC"}.issubset(set(df.columns))
            else Failure("Required columns missing from prediction results CSV.")
        )
    )


def compute_cohort_metrics(
    df: pl.DataFrame, metric_col: str
) -> Result[tuple[MetricSummaryRow, ...], str]:
    """Aggregates seed-level predictions into cohort-level mean and std dev."""
    # Filter only individual registered cohorts
    registered_names = tuple(c.cohort for c in COHORT_REGISTRY)
    df_filtered = df.filter(pl.col("Cohort").is_in(list(registered_names)))

    if df_filtered.is_empty():
        return Failure(f"No records matching registered cohorts for metric {metric_col}")

    aggregated = (
        df_filtered.group_by(["Cohort", "Mode"])
        .agg([
            pl.col(metric_col).mean().alias("mean_val"),
            pl.col(metric_col).std().alias("std_val"),
        ])
    )

    rows: tuple[MetricSummaryRow, ...] = tuple(
        MetricSummaryRow(
            cohort=row["Cohort"],
            mode=row["Mode"],
            metric=metric_col,
            mean_val=row["mean_val"],
            std_val=row["std_val"] if row["std_val"] is not None else 0.0,
            sample_size=find_cohort_meta(row["Cohort"]).value_or(
                CohortMeta(cohort=row["Cohort"], cancer_type="Unknown", sample_size=0, is_melanoma=False)
            ).sample_size,
            is_melanoma=find_cohort_meta(row["Cohort"]).value_or(
                CohortMeta(cohort=row["Cohort"], cancer_type="Unknown", sample_size=0, is_melanoma=False)
            ).is_melanoma,
        )
        for row in aggregated.iter_rows(named=True)
    )
    return Success(rows)


def calculate_weight(row: MetricSummaryRow, method: WeightingMethod) -> float:
    """Calculates weight for a single cohort observation according to method."""
    match method:
        case WeightingMethod.UNWEIGHTED:
            return 1.0
        case WeightingMethod.SAMPLE_SIZE:
            return float(row.sample_size)
        case WeightingMethod.INVERSE_VARIANCE:
            var = (row.std_val ** 2) if row.std_val > EPSILON else (EPSILON ** 2)
            return 1.0 / var
        case WeightingMethod.INVERSE_STD:
            std = row.std_val if row.std_val > EPSILON else EPSILON
            return 1.0 / std
        case _ as unreachable:
            assert_never(unreachable)


def compute_weighted_average(
    rows: Sequence[MetricSummaryRow],
    scope: CohortScope,
    mode: PredictionMode,
    metric_name: str,
    method: WeightingMethod,
) -> Result[WeightedResult, str]:
    """Pure calculation of weighted mean across a filtered subset of cohort rows."""
    subset = tuple(
        r
        for r in rows
        if r.mode == mode.value
        and (
            scope == CohortScope.ALL_COHORTS
            or (scope == CohortScope.MELANOMA and r.is_melanoma)
        )
    )

    if len(subset) == 0:
        return Failure(f"Empty cohort subset for {scope.value} / {mode.value}")

    weights = tuple(calculate_weight(r, method) for r in subset)
    total_weight = sum(weights)
    if total_weight <= 0.0:
        return Failure(f"Total weight is non-positive ({total_weight})")

    weighted_sum = sum(w * r.mean_val for w, r in zip(weights, subset))
    weighted_mean = weighted_sum / total_weight
    total_samples = sum(r.sample_size for r in subset)

    return Success(
        WeightedResult(
            scope=scope,
            mode=mode,
            metric=metric_name,
            weighting=method,
            weighted_mean=weighted_mean,
            total_samples=total_samples,
            n_cohorts=len(subset),
        )
    )


def extract_pooled_benchmark_results(
    df: pl.DataFrame, metric_col: str
) -> tuple[tuple[str, str, float], ...]:
    """Extracts reference mean values for Combined pooled cohorts from dataset."""
    combined_names = ("Combined-All", "Combined-Melanoma", "Combined-RCC")
    df_comb = df.filter(pl.col("Cohort").is_in(list(combined_names)))
    if df_comb.is_empty():
        return ()

    agg = (
        df_comb.group_by(["Cohort", "Mode"])
        .agg([pl.col(metric_col).mean().alias("mean_val")])
    )
    return tuple(
        (str(row["Cohort"]), str(row["Mode"]), float(row["mean_val"]))
        for row in agg.iter_rows(named=True)
    )


def run_all_aggregations(
    df: pl.DataFrame,
    metric_col: str,
) -> Result[tuple[WeightedResult, ...], str]:
    """Coordinates the aggregation of all scopes, modes, and weighting schemes."""
    return compute_cohort_metrics(df, metric_col).map(
        lambda cohort_rows: tuple(
            result.unwrap()
            for scope in (CohortScope.ALL_COHORTS, CohortScope.MELANOMA)
            for mode in (
                PredictionMode.DECONVOLUTED_FRACTIONS,
                PredictionMode.DIRECT_EXPRESSION,
            )
            for method in (
                WeightingMethod.UNWEIGHTED,
                WeightingMethod.SAMPLE_SIZE,
                WeightingMethod.INVERSE_VARIANCE,
                WeightingMethod.INVERSE_STD,
            )
            for result in (
                compute_weighted_average(cohort_rows, scope, mode, metric_col, method),
            )
            if isinstance(result, Success)
        )
    )


def build_comparison_table(
    weighted_results: Sequence[WeightedResult],
) -> tuple[BenchmarkComparison, ...]:
    """Pairs Deconvoluted Fractions and Direct Expression to compute comparative deltas."""
    # Find matching pairs
    comparisons: tuple[BenchmarkComparison, ...] = tuple(
        BenchmarkComparison(
            scope=deconv.scope.value,
            metric=deconv.metric,
            weighting=deconv.weighting.value,
            deconv_auc=deconv.weighted_mean,
            direct_auc=direct.weighted_mean,
            delta_deconv_minus_direct=deconv.weighted_mean - direct.weighted_mean,
        )
        for deconv in weighted_results
        if deconv.mode == PredictionMode.DECONVOLUTED_FRACTIONS
        for direct in weighted_results
        if direct.mode == PredictionMode.DIRECT_EXPRESSION
        and direct.scope == deconv.scope
        and direct.metric == deconv.metric
        and direct.weighting == deconv.weighting
    )
    return comparisons


def convert_results_to_polars(
    results: Sequence[WeightedResult],
) -> pl.DataFrame:
    """Converts sequence of result models to a clean Polars DataFrame."""
    return pl.DataFrame([
        {
            "Scope": r.scope.value,
            "Mode": r.mode.value,
            "Metric": r.metric,
            "Weighting": r.weighting.value,
            "Weighted_Mean": r.weighted_mean,
            "Total_Samples": r.total_samples,
            "N_Cohorts": r.n_cohorts,
        }
        for r in results
    ])


# Imperative Shell


def format_markdown_table(comparisons: Sequence[BenchmarkComparison]) -> str:
    """Renders comparisons as a GitHub-flavored Markdown table."""
    header = (
        "| Cancer Scope | Metric | Weighting Scheme | Deconv Fractions | Direct Expression | Delta (Deconv - Direct) |\n"
        "| :--- | :--- | :--- | :---: | :---: | :---: |"
    )
    rows = "\n".join(
        f"| {c.scope} | {c.metric} | {c.weighting} | {c.deconv_auc:.4f} | {c.direct_auc:.4f} | {c.delta_deconv_minus_direct:+.4f} |"
        for c in comparisons
    )
    return f"{header}\n{rows}"


def format_cohort_breakdown_table(rows: Sequence[MetricSummaryRow]) -> str:
    """Renders individual cohort metrics and error bars."""
    header = (
        "| Cohort | Cancer Type | N Samples | Mode | Metric | Mean | Std (Error Bar) |\n"
        "| :--- | :--- | :---: | :--- | :--- | :---: | :---: |"
    )
    lines = "\n".join(
        f"| {r.cohort} | {r.cancer_type if hasattr(r, 'cancer_type') else find_cohort_meta(r.cohort).value_or(CohortMeta(cohort=r.cohort, cancer_type='-', sample_size=0, is_melanoma=False)).cancer_type} | {r.sample_size} | {r.mode} | {r.metric} | {r.mean_val:.4f} | ±{r.std_val:.4f} |"
        for r in sorted(rows, key=lambda x: (x.cohort, x.mode))
    )
    return f"{header}\n{lines}"


def process_pipeline(
    df: pl.DataFrame, output_csv_path: Path
) -> Result[pl.DataFrame, str]:
    """Processes DataFrame aggregations, writes output CSV, and displays summary."""
    roc_results = run_all_aggregations(df, "ROC_AUC")
    pr_results = run_all_aggregations(df, "PR_AUC")

    if isinstance(roc_results, Failure):
        return Failure(f"ROC_AUC aggregation failed: {roc_results.failure()}")
    if isinstance(pr_results, Failure):
        return Failure(f"PR_AUC aggregation failed: {pr_results.failure()}")

    roc_list = roc_results.unwrap()
    pr_list = pr_results.unwrap()
    all_results = roc_list + pr_list
    df_out = convert_results_to_polars(all_results)

    # Ensure directory exists and write output
    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    df_out.write_csv(output_csv_path)
    print(f"Saved results summary to: {output_csv_path}\n")

    # Generate and print console summaries
    roc_comparisons = build_comparison_table(roc_list)
    pr_comparisons = build_comparison_table(pr_list)

    print("================================================================================")
    print("ROC_AUC COMPARISONS ACROSS WEIGHTING SCHEMES")
    print("================================================================================")
    print(format_markdown_table(roc_comparisons))
    print()

    print("================================================================================")
    print("PR_AUC (PRECISION-RECALL AUC) COMPARISONS ACROSS WEIGHTING SCHEMES")
    print("================================================================================")
    print(format_markdown_table(pr_comparisons))
    print()

    # Also display reference pooled cohorts
    pooled_roc = extract_pooled_benchmark_results(df, "ROC_AUC")
    print("================================================================================")
    print("REFERENCE POOLED CROSS-VALIDATION MODELS (ROC_AUC)")
    print("================================================================================")
    print("| Pooled Cohort | Mode | Pooled CV ROC_AUC |")
    print("| :--- | :--- | :---: |")
    for c_name, m_name, val in pooled_roc:
        print(f"| {c_name} | {m_name} | {val:.4f} |")
    print()

    return Success(df_out)


def execute_pipeline(
    input_path: Path, output_csv_path: Path
) -> Result[pl.DataFrame, str]:
    """Orchestrates loading, calculation, saving, and reporting."""
    print(f"Loading data from: {input_path}")
    return load_prediction_dataset(input_path).bind(
        lambda df: process_pipeline(df, output_csv_path)
    )


def main() -> None:
    """Entry point for command-line execution."""
    base_dir = Path(__file__).resolve().parent.parent
    input_file = base_dir / "output/bagaev-constrained-deconv/bagaev_prediction_results.csv"
    output_file = base_dir / "output/bagaev-constrained-deconv/bagaev_auc_weighting_summary.csv"

    match execute_pipeline(input_file, output_file):
        case Success(_):
            print("Pipeline completed successfully.")
        case Failure(err):
            print(f"Pipeline error: {err}")


if __name__ == "__main__":
    main()
