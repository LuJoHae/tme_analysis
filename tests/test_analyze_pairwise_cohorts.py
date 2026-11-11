"""Unit and Hypothesis property-based tests for scripts/milopy_analysis/analyze_pairwise_cohorts.py.

Verifies:
- Exhaustive PairNoveltyClass tagging and Polars univariate DA join invariants via Hypothesis
- Lesion-level evaluation unit resolution (handling discordant pre/post response biopsies per patient)
- Index-invariant canonical lineage-pair recurrence aggregation across cohorts
- Pure end-to-end execution via dependency-injected CohortDataLoader and FigureExporter protocols
- Valid Nature Methods SVG compilation and XML well-formedness via vl-convert
"""

from __future__ import annotations

from pathlib import Path
import xml.etree.ElementTree as ET

from hypothesis import given, settings, strategies as st
import numpy as np
import polars as pl
from returns.maybe import Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp
import vl_convert as vlc

from scripts.milopy_analysis.analyze_pairwise_cohorts import (
    AltairChart,
    CachedCohortInputs,
    CohortPairwiseRunConfig,
    CohortSortOrder,
    PairNoveltyClass,
    aggregate_recurrent_celltype_pairs,
    annotate_pairs_with_univariate_da,
    build_cohort_quadrant_volcano_chart,
    build_cross_cohort_pairwise_comparison_chart,
    build_evaluation_design,
    build_recurrent_celltype_interactions_chart,
    classify_pair_novelty,
    parse_cli_args,
    render_markdown_summary_report,
    run_all_cohorts_pipeline,
)


def _synthetic_cohort_inputs(
    accession: str,
    swap_uv: bool = False,
    seed: int = 42,
) -> CachedCohortInputs:
    """Constructs a synthetic cached cohort with a clear cross-lineage +- / -+ synergy."""
    rng = np.random.default_rng(seed)
    n_evals = 60  # 30 R, 30 NR
    cells_per_eval = 10
    n_cells = n_evals * cells_per_eval
    n_nhoods = 6

    pt_ids: list[str] = []
    resps: list[str] = []
    cell_types: list[str] = []
    rows: list[int] = []
    cols: list[int] = []

    # Responders have CD8_T present & Treg missing; Non-responders have both or Treg only
    idx_cd8 = 1 if swap_uv else 0
    idx_treg = 0 if swap_uv else 1

    for j in range(n_evals):
        is_r = j < 30
        # Test discordant pre/post patient IDs on j=0 and j=30 sharing patient_id "P_DISCORDANT"
        pt = "P_DISCORDANT" if j in (0, 30) else f"PT_{j:03d}"
        resp = "responder" if is_r else "non-responder"
        if is_r:
            # 26 R in +-, 2 R in ++, 1 R in -+, 1 R in --
            has_cd8 = j < 28
            has_treg = j in (26, 27, 28)
        else:
            # 2 NR in +-, 10 NR in ++ (30..39), 8 NR in -+ (40..47), 10 NR in -- (48..57)
            has_cd8 = (30 <= j < 40) or (58 <= j < 60)
            has_treg = 30 <= j < 48

        for c in range(cells_per_eval):
            cid = j * cells_per_eval + c
            pt_ids.append(pt)
            resps.append(resp)
            if c < 3:
                cell_types.append("CD8_T")
                if has_cd8:
                    rows.append(cid)
                    cols.append(idx_cd8)
            elif c < 6:
                cell_types.append("Treg")
                if has_treg:
                    rows.append(cid)
                    cols.append(idx_treg)
            else:
                cell_types.append("Macrophage")
                nh_k = 2 + (c - 6)
                if nh_k < n_nhoods and bool(rng.random() < 0.45):
                    rows.append(cid)
                    cols.append(nh_k)

    nhoods_mat = sp.csc_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, cols)),
        shape=(n_cells, n_nhoods),
    )
    cells_df = pl.DataFrame(
        {
            "cell_id": [f"cell_{i}" for i in range(n_cells)],
            "patient_id": pt_ids,
            "clinical_response": resps,
            "cell_type": cell_types,
        }
    )
    # Mark both nhoods 0 and 1 as univariate NS (FDR=0.40) to verify de_novo tagging
    da_df = pl.DataFrame(
        {
            "Nhood": list(range(n_nhoods)),
            "logFC": [0.3, -0.4, 0.05, -0.02, 0.1, -0.1],
            "FDR": [0.40, 0.35, 0.90, 0.92, 0.85, 0.88],
        }
    )
    return CachedCohortInputs(
        accession=accession,
        nhoods_mat=nhoods_mat,
        cells_df=cells_df,
        univariate_da_df=da_df,
    )


@settings(max_examples=40, deadline=None)
@given(
    fdr_u=st.floats(min_value=0.0, max_value=1.0, allow_nan=False),
    fdr_v=st.floats(min_value=0.0, max_value=1.0, allow_nan=False),
    thresh=st.floats(min_value=0.01, max_value=0.50, allow_nan=False),
)
def test_prop_novelty_classification_matches_univariate_join(
    fdr_u: float,
    fdr_v: float,
    thresh: float,
) -> None:
    sig_u = fdr_u < thresh
    sig_v = fdr_v < thresh
    expected = classify_pair_novelty(sig_u, sig_v)

    pairs_df = pl.DataFrame({"nhood_u": [0], "nhood_v": [1], "is_significant": [True]})
    da_df = pl.DataFrame({"Nhood": [0, 1], "logFC": [1.2, -0.8], "FDR": [fdr_u, fdr_v]})
    annotated = annotate_pairs_with_univariate_da(pairs_df, da_df, fdr_threshold=thresh)

    assert annotated.height == 1
    assert annotated["novelty_class"][0] == expected.value
    assert bool(annotated["is_univariate_sig_u"][0]) == sig_u
    assert bool(annotated["is_univariate_sig_v"][0]) == sig_v


def test_lesion_evaluation_design_splits_discordant_patient_responses() -> None:
    inputs = _synthetic_cohort_inputs("GSE120575", swap_uv=False, seed=11)
    design_res = build_evaluation_design(inputs)
    assert isinstance(design_res, Success)
    design = design_res.unwrap()

    # Patient "P_DISCORDANT" has both responder and non-responder biopsies -> 60 distinct evaluations
    assert design.n_responders == 30
    assert design.n_non_responders == 30
    assert "P_DISCORDANT__responder" in design.eval_ids
    assert "P_DISCORDANT__non-responder" in design.eval_ids
    assert design.count_matrix.shape == (6, 60)
    assert np.all(design.library_sizes == 10.0)


def test_end_to_end_pure_pipeline_recurrence_and_svg_xml(tmp_path: Path) -> None:
    cohort_a = _synthetic_cohort_inputs("GSE120575", swap_uv=False, seed=21)
    cohort_b = _synthetic_cohort_inputs("GSE243013", swap_uv=True, seed=22)
    store = {"GSE120575": cohort_a, "GSE243013": cohort_b}

    def memory_loader(accession: str, input_dir: Path) -> Result[CachedCohortInputs, str]:
        _ = input_dir
        return (
            Success(store[accession])
            if accession in store
            else Failure(f"Unknown cohort {accession}")
        )

    exported_svgs: list[str] = []

    def memory_svg_exporter(chart: AltairChart, base_path: Path) -> Result[tuple[Path, Path], str]:
        svg_str = str(vlc.vegalite_to_svg(chart.to_dict()))
        # Parse XML strictly to prove well-formedness
        root = ET.fromstring(svg_str)
        assert root.tag.endswith("svg")
        exported_svgs.append(svg_str)
        return Success((base_path.with_suffix(".svg"), base_path.with_suffix(".png")))

    cfg = parse_cli_args(
        [
            "--cohorts",
            "GSE120575,GSE243013",
            "--output-dir",
            str(tmp_path / "pairwise_out"),
            "--sort-order",
            CohortSortOrder.DISCOVERY_COUNT.value,
        ]
    )
    assert cfg.cohort_filter == Some(("GSE120575", "GSE243013"))

    pipe_res = run_all_cohorts_pipeline(
        config=cfg,
        loader=memory_loader,
        exporter=memory_svg_exporter,
        persist_outputs=True,
    )
    assert isinstance(pipe_res, Success)
    bundles, summary_df, recurrent_df = pipe_res.unwrap()

    assert len(bundles) == 2
    assert summary_df.height == 2
    # Both cohorts had CD8_T+ / Treg- swapped in index order (0,1 vs 1,0); canonicalization must merge them!
    cd8_treg = recurrent_df.filter(
        (pl.col("lineage_primary") == "CD8_T")
        & (pl.col("lineage_secondary") == "Treg")
        & (pl.col("regime") == "+- (First Present / Second Missing)")
    )
    assert cd8_treg.height >= 1
    assert int(cd8_treg["n_cohorts"][0]) == 2
    assert int(cd8_treg["total_de_novo_pairs"][0]) >= 2

    # 2 per-cohort volcano SVGs + 2 cross-cohort report SVGs = 4 validated XML SVGs
    assert len(exported_svgs) == 4

    # Verify persisted Parquet and Markdown files
    assert (tmp_path / "pairwise_out" / "GSE120575" / "quadrant_pairs.parquet").exists()
    assert (tmp_path / "pairwise_out" / "reports" / "pairwise_cohorts_summary.parquet").exists()
    md_path = tmp_path / "pairwise_out" / "reports" / "milopy_pairwise_cohorts_summary.md"
    assert md_path.exists()
    md_content = render_markdown_summary_report(summary_df, recurrent_df, bundles, cfg)
    assert "CD8_T" in md_content and "Treg" in md_content

    # Also verify empty chart branches compile cleanly to XML SVG
    vol_empty = build_cohort_quadrant_volcano_chart(
        bundles[0].quadrant_df.clear(),
        bundles[0].summary,
        0.10,
    ).unwrap()
    assert ET.fromstring(str(vlc.vegalite_to_svg(vol_empty.to_dict()))).tag.endswith("svg")
    comp_chart = build_cross_cohort_pairwise_comparison_chart(summary_df).unwrap()
    assert ET.fromstring(str(vlc.vegalite_to_svg(comp_chart.to_dict()))).tag.endswith("svg")
    rec_empty = build_recurrent_celltype_interactions_chart(bundles, recurrent_df.clear()).unwrap()
    assert ET.fromstring(str(vlc.vegalite_to_svg(rec_empty.to_dict()))).tag.endswith("svg")
