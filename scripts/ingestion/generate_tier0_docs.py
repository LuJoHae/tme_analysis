#!/usr/bin/env python3
"""Generate comprehensive Markdown documentation for Tier 0 Benchmark Core Cohorts."""

from pathlib import Path
import polars as pl


def main() -> None:
    df = pl.read_parquet("data/registry/discovered_solid_tumor_sc_datasets.parquet")
    t0 = df.filter(pl.col("tier") == "Tier 0 (Benchmark Core)").sort(
        ["indication", "tier0_score"], descending=[False, True]
    )

    md: list[str] = []
    md.append("# Tier 0: Premier Benchmark Core scRNA-seq Cancer Cohorts\n")
    md.append(
        "This document details the **22 ultra-select Tier 0 cohorts** selected from the 352-cohort discovery catalog. "
        "These cohorts represent high-powered clinical immunotherapy trials with verified response outcomes, "
        "high patient numbers ($N \\ge 9$, median $N = 51$), and broad cellular microenvironment representation across "
        "9 human solid tumor indications.\n"
    )

    unselected_cnt = len(t0.filter(pl.col("cell_selection_strategy").str.contains("Unselected")))
    pct_unselected = (unselected_cnt / len(t0)) * 100

    md.append("## Summary Statistics")
    md.append(f"- **Total Selected Cohorts**: {len(t0)}")
    md.append(f"- **Total Clinically Annotated Patients/Samples**: {t0['sample_or_patient_count'].sum():,}")
    md.append(f"- **Total Single-Cell/Nucleus Transcriptomes**: ~{t0['cell_count_estimate'].sum():,}")
    md.append(f"- **Unselected Whole-TME Suspensions**: {unselected_cnt} cohorts ({pct_unselected:.0f}%)\n")

    md.append("## Cohort Master Table\n")
    md.append(
        "| Accession | Indication | Score | Patients | Cells | Modality & Tech | Cell Selection | Clinical Response Details |"
    )
    md.append(
        "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |"
    )

    for row in t0.iter_rows(named=True):
        acc = row["accession"]
        ind = row["indication"]
        score = row["tier0_score"]
        pts = row["sample_or_patient_count"]
        cells = f"{row['cell_count_estimate']:,}"
        mod = f"{row['modality']} ({row['technology']})"
        strat = row["cell_selection_strategy"]
        resp = row["response_details"].replace("|", "/")
        md.append(f"| **{acc}** | {ind} | **{score}/100** | {pts} | {cells} | {mod} | {strat} | {resp} |")

    md.append("\n---\n")
    md.append("## Detailed Cohort Profiles\n")

    for idx, row in enumerate(t0.iter_rows(named=True), 1):
        acc = row["accession"]
        ind = row["indication"]
        score = row["tier0_score"]
        pts = row["sample_or_patient_count"]
        cells = f"{row['cell_count_estimate']:,}"
        title = row["title"]
        strat = row["cell_selection_strategy"]
        quote = row["cell_selection_quote"]
        url = row["accession_url"]
        summary = row["summary"]
        mat = row["matrix_files"]
        rat = row["tier0_rationale"]

        md.append(f"### {idx}. {acc} — {ind} (Score: {score}/100)")
        md.append(f"**Title**: {title}\n")
        md.append(f"- **Accession URL**: [{acc}]({url})")
        md.append(f"- **Indication**: {ind}")
        md.append(f"- **Sample/Patient Count**: {pts} patients / biological specimens")
        md.append(f"- **Cell Count Estimate**: ~{cells} single cells/nuclei")
        md.append(f"- **Cell Selection Strategy**: `{strat}`")
        md.append(f"- **Protocol Excerpt**: > *\"{quote}\"*")
        md.append(f"- **Scoring Rationale**: `{rat}`")
        md.append(f"- **Available Raw Matrices**: `{mat}`\n")
        md.append(f"**Study Abstract / Clinical Setting**:\n{summary}\n")

    out_p = Path("docs/tier0_benchmark_core_cohorts.md")
    out_p.write_text("\n".join(md))
    print(f"Written complete profile to {out_p}")


if __name__ == "__main__":
    main()
