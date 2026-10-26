#!/usr/bin/env python3
"""Programmatic Scoring & Prioritization Ranker for scRNA-seq Cancer Cohorts.

Implements the 100-point multi-dimensional rubric and 4 hard gates approved in
plan_premier_cohort_selection.md to systematically score and select
Tier 0 (Benchmark Core) datasets across solid tumor indications.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import polars as pl


@dataclass(frozen=True)
class CohortEvaluation:
    """Evaluation result for a single scRNA-seq cohort."""

    accession: str
    indication: str
    title: str
    technology: str
    cell_selection: str
    sample_or_patient_count: int
    cell_count_estimate: int
    passed_hard_gates: bool
    gate_failure_reason: str
    clinical_score: int       # max 30
    scale_score: int          # max 20
    kinetics_score: int       # max 15
    breadth_score: int        # max 15
    matrix_score: int         # max 10
    impact_score: int         # max 10
    total_score: int          # max 100
    rationale: str


def evaluate_cohort(row: dict[str, object]) -> CohortEvaluation:
    """Evaluate a single cohort against hard gates and scoring dimensions."""
    acc = str(row["accession"])
    ind = str(row["indication"])
    title = str(row.get("title", ""))
    summary = str(row.get("summary", ""))
    combined_text = f"{title} {summary}".lower()
    resp_details = str(row.get("response_details", "")).lower()
    cell_strat = str(row.get("cell_selection_strategy", ""))
    mat_files = str(row.get("matrix_files", "")).lower()
    doi = str(row.get("doi", ""))
    is_resp_annot = bool(row.get("clinical_response_annotated", False))

    raw_sample_cnt = row.get("sample_or_patient_count", 0)
    sample_cnt = int(raw_sample_cnt) if raw_sample_cnt is not None else 0

    raw_cell_cnt = row.get("cell_count_estimate", 0)
    cell_cnt = int(raw_cell_cnt) if raw_cell_cnt is not None else 0

    # 1. Evaluate Hard Gates
    gate_reasons: list[str] = []

    # Gate 1: In vivo human tumor (exclude cell line-only, organoid-only, mouse xenograft-only)
    in_vitro_markers = ["cell line", "cell lines", "organoid", "organoids", "xenograft", "zebrafish"]
    if any(m in combined_text for m in in_vitro_markers) and sample_cnt < 8:
        gate_reasons.append("In vitro model / cell line without sufficient patient biopsy cohort")

    # Gate 2: Statistical power minimum (minimum >= 8 patients / samples; exclude single-patient case reports)
    case_report_markers = ["in a 16-month-old patient", "case report", "in a single patient", "of a single patient"]
    if any(m in combined_text for m in case_report_markers):
        gate_reasons.append("Single-patient case study with multiple tissue slices (N=1 patient)")
    elif sample_cnt < 8:
        gate_reasons.append(f"Sample size below power threshold ({sample_cnt} < 8)")

    # Gate 3: Clinical outcome annotation
    if not is_resp_annot and "response" not in resp_details and "treatment" not in resp_details:
        gate_reasons.append("Missing verifiable clinical outcome linkage")

    passed_gates = len(gate_reasons) == 0
    failure_reason = "; ".join(gate_reasons) if not passed_gates else "Passed all hard gates"

    # 2. Score Dimensions
    # Dimension 1: Clinical Rigor & RECIST Ground Truth (Max 30)
    recist_keywords = [
        "recist", "cr", "pr", "sd", "pd", "progression-free", "pfs", "overall survival",
        "complete response", "partial response", "irrecist", "objective response", "clinical benefit"
    ]
    binary_keywords = ["responder", "non-responder", "response", "resistance", "durable", "efficacy"]

    if any(k in combined_text or k in resp_details for k in recist_keywords):
        clinical_score = 30
    elif any(k in combined_text or k in resp_details for k in binary_keywords):
        clinical_score = 20
    else:
        clinical_score = 10

    # Dimension 2: Patient Scale & Statistical Power (Max 20)
    if sample_cnt >= 30:
        scale_score = 20
    elif sample_cnt >= 15:
        scale_score = 15
    elif sample_cnt >= 8:
        scale_score = 10
    else:
        scale_score = 4

    if cell_cnt >= 50000 and scale_score < 20:
        scale_score = min(20, scale_score + 2)

    # Dimension 3: Biopsy Kinetics (Max 15)
    longitudinal_keywords = [
        "longitudinal", "paired", "sequential", "pre- and post", "pre- and on",
        "before and after", "pre/post", "on-treatment", "time points", "serial"
    ]
    baseline_keywords = ["baseline", "pre-treatment", "treatment-naïve", "treatment-naive", "prior to"]

    if any(k in combined_text for k in longitudinal_keywords):
        kinetics_score = 15
    elif any(k in combined_text for k in baseline_keywords):
        kinetics_score = 10
    elif "resistant" in combined_text or "progression" in combined_text:
        kinetics_score = 8
    else:
        kinetics_score = 5

    # Dimension 4: Cellular TME Breadth (Max 15)
    if "Unselected" in cell_strat:
        breadth_score = 15
    elif "Nuclei" in cell_strat:
        breadth_score = 13
    elif "CD45+" in cell_strat or "Leukocyte" in cell_strat:
        breadth_score = 10
    else:
        breadth_score = 6

    # Dimension 5: Assay & Matrix Availability (Max 10)
    if ".h5ad" in mat_files or ".mtx.gz" in mat_files or ".h5" in mat_files:
        matrix_score = 10
    elif ".tsv.gz" in mat_files or ".csv.gz" in mat_files or ".txt.gz" in mat_files:
        matrix_score = 9
    elif "raw.tar" in mat_files:
        matrix_score = 8
    else:
        matrix_score = 6

    # Dimension 6: Impact & Indication Balance (Max 10)
    impact_score = 5  # Base representation
    if doi or "pmid" in doi.lower() or "pmid" in combined_text:
        impact_score += 5

    total_score = clinical_score + scale_score + kinetics_score + breadth_score + matrix_score + impact_score

    rationale_parts = [
        f"Clin:{clinical_score}/30",
        f"Scale:{scale_score}/20 (N={sample_cnt})",
        f"Kinetics:{kinetics_score}/15",
        f"Breadth:{breadth_score}/15 ({cell_strat})",
        f"Matrix:{matrix_score}/10",
        f"Impact:{impact_score}/10",
    ]
    rationale = " | ".join(rationale_parts)

    return CohortEvaluation(
        accession=acc,
        indication=ind,
        title=title,
        technology=str(row.get("technology", "")),
        cell_selection=cell_strat,
        sample_or_patient_count=sample_cnt,
        cell_count_estimate=cell_cnt,
        passed_hard_gates=passed_gates,
        gate_failure_reason=failure_reason,
        clinical_score=clinical_score,
        scale_score=scale_score,
        kinetics_score=kinetics_score,
        breadth_score=breadth_score,
        matrix_score=matrix_score,
        impact_score=impact_score,
        total_score=total_score,
        rationale=rationale,
    )


def rank_cohorts(
    registry_path: Path,
    target_count: int = 18,
) -> tuple[list[CohortEvaluation], list[CohortEvaluation]]:
    """Rank Tier 1 cohorts and select the balanced Tier 0 cohort group."""
    df = pl.read_parquet(registry_path)
    tier1_df = df.filter(pl.col("tier").str.starts_with("Tier 1"))

    evals = [evaluate_cohort(row) for row in tier1_df.iter_rows(named=True)]

    # Separate into passed gates vs failed gates
    passed = [e for e in evals if e.passed_hard_gates]
    failed = [e for e in evals if not e.passed_hard_gates]

    # Sort passed candidates by total_score descending, then sample count descending FIRST
    passed.sort(key=lambda x: (x.total_score, x.sample_or_patient_count, x.cell_count_estimate), reverse=True)

    # Deduplicate studies with identical core title within same indication (e.g. SubSeries and SuperSeries)
    seen_titles: set[str] = set()
    deduped_passed: list[CohortEvaluation] = []
    for e in passed:
        norm_title = f"{e.indication}_{''.join(c for c in e.title.lower() if c.isalnum())[:55]}"
        if norm_title in seen_titles:
            continue
        seen_titles.add(norm_title)
        deduped_passed.append(e)

    passed = deduped_passed

    # Balanced selection: Ensure top 1-2 per indication, then fill remaining by top score respecting max_per_indication
    max_per_indication = 3
    indications = sorted(list(set(e.indication for e in passed)))
    selected: list[CohortEvaluation] = []
    selected_accs: set[str] = set()
    ind_counts: dict[str, int] = {ind: 0 for ind in indications}

    # Pass 1: Select top cohort for each represented indication
    for ind in indications:
        ind_cohorts = [e for e in passed if e.indication == ind]
        if ind_cohorts:
            top_ind = ind_cohorts[0]
            selected.append(top_ind)
            selected_accs.add(top_ind.accession)
            ind_counts[ind] += 1

    # Pass 2: Select second top cohort for key high-prevalence indications (Melanoma, NSCLC, Breast, CRC, ccRCC, HCC, HNSCC)
    key_indications = {"Melanoma", "NSCLC", "Breast", "CRC", "ccRCC", "HCC", "HNSCC"}
    for ind in key_indications:
        if ind_counts.get(ind, 0) < max_per_indication:
            ind_cohorts = [e for e in passed if e.indication == ind and e.accession not in selected_accs]
            if ind_cohorts:
                top_ind = ind_cohorts[0]
                selected.append(top_ind)
                selected_accs.add(top_ind.accession)
                ind_counts[ind] += 1

    # Pass 3: Fill remaining slots up to target_count using highest total score with max_per_indication cap
    remaining = [e for e in passed if e.accession not in selected_accs]
    remaining.sort(key=lambda x: (x.total_score, x.sample_or_patient_count), reverse=True)
    for e in remaining:
        if len(selected) >= target_count:
            break
        if ind_counts[e.indication] < max_per_indication:
            selected.append(e)
            selected_accs.add(e.accession)
            ind_counts[e.indication] += 1

    # Sort final selected list by Indication, then score descending
    selected.sort(key=lambda x: (x.indication, -x.total_score))

    return selected, evals


def update_registry_with_tier0(
    registry_path: Path,
    tsv_path: Path,
    tier0_accessions: Sequence[str],
    eval_dict: dict[str, CohortEvaluation],
) -> None:
    """Update parquet and TSV registry with Tier 0 classification and scores."""
    df = pl.read_parquet(registry_path)

    # Compute new tier and score columns
    new_tiers: list[str] = []
    scores: list[int] = []
    rationales: list[str] = []

    for row in df.iter_rows(named=True):
        acc = str(row["accession"])
        curr_tier = str(row["tier"])
        if acc in tier0_accessions:
            new_tiers.append("Tier 0 (Benchmark Core)")
        else:
            new_tiers.append(curr_tier)

        ev = eval_dict.get(acc)
        if ev:
            scores.append(ev.total_score)
            rationales.append(ev.rationale)
        else:
            scores.append(0)
            rationales.append("Not evaluated (Tier 2 Baseline Atlas)")

    updated_df = df.with_columns(
        pl.Series("tier", new_tiers),
        pl.Series("tier0_score", scores),
        pl.Series("tier0_rationale", rationales),
    )

    updated_df.write_parquet(registry_path)
    updated_df.write_csv(tsv_path, separator="\t")
    print(f"Updated registry: {registry_path} and {tsv_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Rank scRNA-seq cohorts to select Tier 0 Premier Core.")
    parser.add_argument(
        "--registry",
        type=Path,
        default=Path("data/registry/discovered_solid_tumor_sc_datasets.parquet"),
        help="Path to Parquet registry.",
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        default=Path("data/registry/discovered_solid_tumor_sc_datasets.tsv"),
        help="Path to TSV registry.",
    )
    parser.add_argument(
        "--target-count",
        type=int,
        default=18,
        help="Target number of Tier 0 benchmark cohorts (default: 18).",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Write updated Tier 0 classifications directly into the registry.",
    )

    args = parser.parse_args()

    selected, all_evals = rank_cohorts(registry_path=args.registry, target_count=args.target_count)
    eval_dict = {e.accession: e for e in all_evals}

    print(f"\n=======================================================")
    print(f"  Tier 0 (Premier Benchmark Core) Selection Results")
    print(f"  Total Tier 1 Evaluated: {len(all_evals)}")
    print(f"  Passed Hard Gates:      {sum(1 for e in all_evals if e.passed_hard_gates)}")
    print(f"  Selected for Tier 0:    {len(selected)} cohorts")
    print(f"=======================================================\n")

    print(f"{'Accession':<12} | {'Indication':<10} | {'Score':<5} | {'Patients':<8} | {'Cells':<8} | {'Cell Selection':<35} | Title")
    print("-" * 120)
    for s in selected:
        short_title = (s.title[:45] + "...") if len(s.title) > 45 else s.title
        print(f"{s.accession:<12} | {s.indication:<10} | {s.total_score:<5} | {s.sample_or_patient_count:<8} | {s.cell_count_estimate:<8} | {s.cell_selection[:35]:<35} | {short_title}")

    if args.apply:
        tier0_accs = [s.accession for s in selected]
        update_registry_with_tier0(
            registry_path=args.registry,
            tsv_path=args.tsv,
            tier0_accessions=tier0_accs,
            eval_dict=eval_dict,
        )


if __name__ == "__main__":
    main()
