"""CLI Script to verify milopy compatibility across all preprocessed single-cell .h5ad datasets."""

import argparse
from pathlib import Path
import altair as alt
import polars as pl
from returns.result import Success, Failure
from single_cell_immuno_datasets import DataDirectories, verify_milopy_compatibility, MilopyCheckResult


def create_readiness_chart(df: pl.DataFrame, out_svg_path: Path) -> Path:
    """Generates Altair SVG bar chart showing cell counts per dataset colored by milopy readiness status."""
    pandas_df = df.to_pandas()
    
    chart = alt.Chart(pandas_df).mark_bar().encode(
        x=alt.X("accession:N", title="Dataset Accession", sort="-y"),
        y=alt.Y("n_obs:Q", title="Total Cells (n_obs)"),
        color=alt.Color("status_summary:N", title="Milopy Readiness", scale=alt.Scale(
            domain=["PASS_MILOPY_READY", "FAIL_SINGLE_SAMPLE", "FAIL_NO_DESIGN_CONTRAST", "FAIL_INSUFFICIENT_CELLS", "FAIL_DRY_RUN_ERROR"],
            range=["#2ca02c", "#ff7f0e", "#1f77b4", "#d62728", "#9467bd"]
        )),
        tooltip=["accession", "n_obs", "n_vars", "n_samples", "n_conditions", "status_summary", "error_detail"]
    ).properties(
        title="ICB Single-Cell Datasets: milopy DA Compatibility & Cell Counts",
        width=600,
        height=350
    )
    
    out_svg_path.parent.mkdir(parents=True, exist_ok=True)
    chart.save(str(out_svg_path))
    return out_svg_path


def generate_markdown_report(results: list[MilopyCheckResult], out_md_path: Path) -> Path:
    """Generates markdown report summarizing milopy compatibility for all preprocessed datasets."""
    lines = [
        "# milopy Neighborhood Differential Abundance Readiness Report",
        "",
        "## Preprocessed ICB Single-Cell Datasets milopy Verification Summary",
        "",
        "| # | Dataset Accession | Cells (n_obs) | Genes (n_vars) | Sparse X? | PCA / KNN? | Sample Col (n_samples) | Design Col (n_cond) | milopy Dry-Run | Overall Status | Error Details |",
        "|:---:|:---|:---:|:---:|:---:|:---:|:---|:---|:---:|:---|:---|",
    ]
    
    for idx, r in enumerate(results, 1):
        pca_knn_str = f"{'Yes' if r.has_pca else 'No'} / {'Yes' if r.has_knn else 'No'}"
        dry_str = "PASS" if r.milopy_dry_run_pass else "FAIL"
        sparse_str = "Yes" if r.has_sparse_x else "No"
        err_str = r.error_detail if r.error_detail else "None"
        
        status_badge = f"**{r.status_summary}**" if "PASS" in r.status_summary else f"`{r.status_summary}`"
        
        lines.append(
            f"| {idx} | {r.accession} | {r.n_obs:,} | {r.n_vars:,} | {sparse_str} | {pca_knn_str} | {r.sample_col} ({r.n_samples}) | {r.design_col} ({r.n_conditions}) | {dry_str} | {status_badge} | {err_str} |"
        )
        
    lines.extend([
        "",
        "---",
        "",
        "### milopy Compatibility Requirements & Guidance",
        "- **PASS_MILOPY_READY**: Matrix is sparse, contains $\\ge 2$ samples and $\\ge 2$ design conditions, and successfully ran `milopy` neighborhood building & GLM testing.",
        "- **FAIL_SINGLE_SAMPLE**: Only 1 unique sample detected. milopy requires $\\ge 2$ samples for GLM differential abundance testing.",
        "- **FAIL_NO_DESIGN_CONTRAST**: Outcome column contains < 2 condition levels.",
        "",
    ])
    
    out_md_path.parent.mkdir(parents=True, exist_ok=True)
    out_md_path.write_text("\n".join(lines), encoding="utf-8")
    return out_md_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify milopy compatibility across all preprocessed single-cell datasets.")
    parser.add_argument("--base-dir", type=str, default="/storage/halu/data", help="Base directory (default: /storage/halu/data)")
    args = parser.parse_args()

    dirs = DataDirectories.with_base(Path(args.base_dir))
    h5ad_files = sorted(list(dirs.preprocessed_dir.glob("*.h5ad")))

    if not h5ad_files:
        print(f"Error: No .h5ad files found in {dirs.preprocessed_dir}")
        return

    print(f"Checking milopy compatibility for {len(h5ad_files)} preprocessed datasets in {dirs.preprocessed_dir}...")
    results: list[MilopyCheckResult] = []

    for f in h5ad_files:
        match verify_milopy_compatibility(f):
            case Success(r):
                print(f" - [{r.accession}]: {r.status_summary} ({r.n_obs} cells, {r.n_samples} samples, {r.n_conditions} conditions)")
                results.append(r)
            case Failure(err):
                print(f" - [{f.name}]: Verification error: {err}")

    # Create Polars summary dataframe
    df = pl.DataFrame([r.model_dump() for r in results])
    
    report_md = dirs.reports_dir / "milopy_compatibility_report.md"
    generate_markdown_report(results, report_md)
    print(f"Saved milopy compatibility report -> {report_md}")

    svg_chart = dirs.reports_dir / "milopy_readiness_summary.svg"
    create_readiness_chart(df, svg_chart)
    print(f"Saved milopy readiness SVG chart -> {svg_chart}")


if __name__ == "__main__":
    main()
