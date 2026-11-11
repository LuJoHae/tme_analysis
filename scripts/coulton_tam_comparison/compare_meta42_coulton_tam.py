#!/usr/bin/env python3
"""Main CLI & imperative shell comparing Consensus-Meta42 states to Coulton_TAM_Atlas substates.

Architecturally follows Functional Core, Imperative Shell (.agents/rules/code-style-guide.md):
- Loads `Coulton_TAM_Atlas` and the 11 clinical bulk ICB cohorts (`N = 1,058`) strictly via
  the `tme_datasets` public API (`load_coulton_tam_atlas`, `load_dataset`, `get_data_paths`).
- Streams sparse CSR batches (`25,000` cells/chunk) into pure out-of-core accumulators so peak
  memory stays `< 4.5 GB` on 16 GB machines.
- Delegates all similarity algebra, marker discovery, penalized Newton-Raphson bulk ICB
  meta-analysis, and Markdown synthesis to `scripts.coulton_tam_comparison.core`.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import gc
import os
from pathlib import Path
import sys
from typing import Final, Literal, assert_never

import numpy as np
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp

_REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tme_datasets import (  # noqa: E402
    get_data_paths,
    load_coulton_tam_atlas,
    load_dataset,
)

from scripts.coulton_tam_comparison.core import (  # noqa: E402
    accumulate_coulton_batch,
    create_empty_coulton_accumulator,
    deduplicate_gene_columns,
    deserialize_coulton_centroids_from_polars,
    finalize_coulton_centroids,
    run_comparison_core_pure,
    serialize_coulton_centroids_to_polars,
)
from scripts.coulton_tam_comparison.models import (  # noqa: E402
    AssaySource,
    BULK_ICB_COHORT_SPECS,
    BulkIcbCohortBundle,
    BulkIcbSuiteLoader,
    COULTON_ORDERED_CLUSTERS,
    ComparisonConfig,
    ComparisonOutputs,
    CoultonAtlasStreamer,
    CoultonCentroidBundle,
    FloatMatrix,
    Meta42ArtifactLoader,
    Meta42ReferenceBundle,
    VectorChartExporter,
)
from scripts.coulton_tam_comparison.render_figures import (  # noqa: E402
    default_vector_exporter,
    export_all_comparison_figures,
    export_main_paper_typst_caption,
)


def _resolve_repo_path(path: Path) -> Path:
    """Resolve relative paths cleanly against repository root."""
    return path if path.is_absolute() else (_REPO_ROOT / path)


def linearize_and_standardize_bulk_matrix(raw_mat: FloatMatrix) -> FloatMatrix:
    """Pure bulk RNA-seq linearization (`2^X - 1` if log2-scaled), TPM scaling, and `log2(TPM+1)` z-score."""
    clean_mat = np.nan_to_num(
        np.asarray(raw_mat, dtype=np.float64),
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )
    max_val = float(np.max(clean_mat)) if clean_mat.size > 0 else 0.0
    was_log2 = max_val < 50.0
    lin_mat = (
        np.maximum(np.power(2.0, np.clip(clean_mat, 0.0, 32.0)) - 1.0, 0.0)
        if was_log2
        else np.maximum(clean_mat, 0.0)
    )
    row_sums = np.maximum(lin_mat.sum(axis=1, keepdims=True), 1e-12)
    tpm_mat = (lin_mat / row_sums) * 1e6
    log2_mat = np.log2(tpm_mat + 1.0)
    g_mean = log2_mat.mean(axis=0, keepdims=True)
    g_std = np.maximum(log2_mat.std(axis=0, keepdims=True), 1e-6)
    return np.asarray((log2_mat - g_mean) / g_std, dtype=np.float64)


def _assay_to_provider_literal(assay: AssaySource) -> Literal["RNA", "SCT"]:
    """Map AssaySource enum exhaustively to tme_datasets TamAssayChoice literal."""
    match assay:
        case AssaySource.RNA_LOG1P:
            return "RNA"
        case AssaySource.SCT_DATA:
            return "SCT"
        case _ as unreachable:
            assert_never(unreachable)


def stream_coulton_centroids_via_api(
    assay: AssaySource,
    batch_size: int = 25_000,
) -> Result[CoultonCentroidBundle, str]:
    """Load a single Coulton_TAM_Atlas assay selectively via tme_datasets and stream CSR batches."""
    assay_literal = _assay_to_provider_literal(assay)
    match load_coulton_tam_atlas(
        assay=assay_literal,
        auto_download=False,
        include_all_layers=False,
        repo_root=Nothing,
    ):
        case Failure(err):
            return Failure(f"Failed to load Coulton_TAM_Atlas ({assay_literal}): {err}")
        case Success(loaded_adata):
            adata = loaded_adata

    obs_cols = set(adata.obs.columns)
    state_col = (
        "cell_state"
        if "cell_state" in obs_cols
        else ("short.label" if "short.label" in obs_cols else "")
    )
    if not state_col:
        return Failure("Coulton_TAM_Atlas.obs is missing 'cell_state' / 'short.label'")
    if "tissue" not in obs_cols:
        return Failure("Coulton_TAM_Atlas.obs is missing 'tissue'")

    cluster_to_code = {lbl: idx for idx, lbl in enumerate(COULTON_ORDERED_CLUSTERS)}
    raw_states = [str(x).strip() for x in adata.obs[state_col].tolist()]
    cluster_codes = np.asarray(
        [cluster_to_code.get(st, -1) for st in raw_states],
        dtype=np.int64,
    )
    is_tumor_mask = np.asarray(
        [str(t).strip() == "Tumor" for t in adata.obs["tissue"].tolist()],
        dtype=np.bool_,
    )

    var_cols = set(adata.var.columns)
    gene_names_raw = (
        [str(x).strip() for x in adata.var["gene_name"].tolist()]
        if "gene_name" in var_cols
        else [str(x).strip() for x in adata.var_names.tolist()]
    )
    fallback_ids = (
        [str(x).strip() for x in adata.var["original_id"].tolist()]
        if "original_id" in var_cols
        else [str(x).strip() for x in adata.var_names.tolist()]
    )
    resolved_symbols = [
        gn if (gn and gn.upper() not in {"NAN", "NONE", "<NA>"}) else fb
        for gn, fb in zip(gene_names_raw, fallback_ids, strict=False)
    ]

    if assay == AssaySource.SCT_DATA and "in_sct_assay" in var_cols:
        active_col_mask = np.asarray(adata.var["in_sct_assay"].to_numpy(), dtype=np.bool_)
    else:
        active_col_mask = np.ones(adata.n_vars, dtype=np.bool_)

    all_cols_active = bool(np.all(active_col_mask))
    active_col_indices = np.flatnonzero(active_col_mask)
    active_symbols = tuple(resolved_symbols[int(i)] for i in active_col_indices)

    acc = create_empty_coulton_accumulator(
        assay=assay,
        cluster_labels=COULTON_ORDERED_CLUSTERS,
        gene_symbols=active_symbols,
    )

    n_cells = int(adata.n_obs)
    for start in range(0, n_cells, batch_size):
        stop = min(start + batch_size, n_cells)
        csr_batch = sp.csr_matrix(adata.X[start:stop])
        csr_sub = csr_batch if all_cols_active else csr_batch[:, active_col_indices]
        match accumulate_coulton_batch(
            acc=acc,
            batch_csr=csr_sub,
            cluster_codes=cluster_codes[start:stop],
            is_tumor_mask=is_tumor_mask[start:stop],
        ):
            case Failure(b_err):
                return Failure(b_err)
            case Success(next_acc):
                acc = next_acc

    finalized = finalize_coulton_centroids(acc)
    deduped = deduplicate_gene_columns(raw_symbols=active_symbols, bundle=finalized)
    return Success(deduped)


def load_or_compute_coulton_centroids(
    config: ComparisonConfig,
    streamer: CoultonAtlasStreamer = stream_coulton_centroids_via_api,
) -> Result[tuple[CoultonCentroidBundle, ...], str]:
    """Reuse cached Coulton centroids if valid, or stream out-of-core from Coulton_TAM_Atlas."""
    out_dir = _resolve_repo_path(config.output_dir)
    cache_path = out_dir / "coulton_substate_centroids.parquet"
    required_assays: tuple[AssaySource, ...] = (
        (AssaySource.RNA_LOG1P, AssaySource.SCT_DATA)
        if config.include_sct_sensitivity
        else (AssaySource.RNA_LOG1P,)
    )

    if config.reuse_cached_centroids and cache_path.is_file():
        cached_df = pl.read_parquet(cache_path)
        match deserialize_coulton_centroids_from_polars(cached_df):
            case Success(cached_bundles):
                present = {b.assay for b in cached_bundles}
                if all(req in present for req in required_assays):
                    filtered = tuple(b for b in cached_bundles if b.assay in required_assays)
                    return Success(filtered)
            case Failure(_):
                pass

    computed: list[CoultonCentroidBundle] = []
    for req_assay in required_assays:
        print(
            f"[INFO] Streaming Coulton_TAM_Atlas ({req_assay.value}) centroids in "
            f"{config.batch_size:,}-cell CSR chunks..."
        )
        match streamer(assay=req_assay, batch_size=config.batch_size):
            case Failure(err):
                return Failure(err)
            case Success(bundle):
                computed.append(bundle)
        gc.collect()

    out_dir.mkdir(parents=True, exist_ok=True)
    serialize_coulton_centroids_to_polars(computed).write_parquet(cache_path)
    return Success(tuple(computed))


def load_meta42_artifacts_from_disk(
    config: ComparisonConfig,
) -> Result[Meta42ReferenceBundle, str]:
    """Load Consensus-Meta42 reference profiles and ICI effect parquets from manuscript outputs."""
    phi_p = _resolve_repo_path(config.phi_meta_path)
    conc_p = _resolve_repo_path(config.concordance_3b_path)
    hub_p = _resolve_repo_path(config.hub_summary_path)
    eff_p = _resolve_repo_path(config.state_effects_path)
    deg_p = _resolve_repo_path(config.top_degs_path)

    for label, path in (
        ("phi_meta", phi_p),
        ("concordance_3b", conc_p),
        ("hub_summary", hub_p),
        ("state_effects", eff_p),
        ("top_degs", deg_p),
    ):
        if not path.is_file():
            return Failure(f"Missing required Meta42 artifact ({label}): {path}")

    return Success(
        Meta42ReferenceBundle(
            phi_meta_df=pl.read_parquet(phi_p),
            concordance_3b_df=pl.read_parquet(conc_p),
            hub_summary_df=pl.read_parquet(hub_p),
            state_effects_df=pl.read_parquet(eff_p),
            top_degs_df=pl.read_parquet(deg_p),
        )
    )


def _load_single_bulk_icb_cohort(
    cohort_id: str,
    cancer_type: str,
    short_name: str,
    preprocessed_dir: Path,
) -> Result[BulkIcbCohortBundle, str]:
    """Load one bulk RNA-seq ICB cohort strictly via `tme_datasets.load_dataset`."""
    match load_dataset(
        cohort_id,
        base_dir=preprocessed_dir,
        auto_download=False,
        normalize_ensembl=False,
    ):
        case Failure(err):
            return Failure(f"Failed to load bulk ICB cohort '{cohort_id}': {err}")
        case Success(loaded_adata):
            adata = loaded_adata

    if "response_binary" not in adata.obs.columns:
        return Failure(f"Bulk cohort '{cohort_id}' lacks 'response_binary' in .obs")

    y_raw = np.asarray(adata.obs["response_binary"].to_numpy(), dtype=np.float64)
    valid_mask = np.isfinite(y_raw) & ((y_raw == 0.0) | (y_raw == 1.0))
    n_valid = int(np.sum(valid_mask))
    if n_valid < 5:
        return Failure(f"Bulk cohort '{cohort_id}' has only {n_valid} evaluable patients")

    x_dense = (
        np.asarray(adata.X.toarray(), dtype=np.float64)
        if hasattr(adata.X, "toarray")
        else np.asarray(adata.X, dtype=np.float64)
    )
    x_valid = np.asarray(x_dense[valid_mask, :], dtype=np.float64)
    y_valid = np.asarray(y_raw[valid_mask], dtype=np.float64)

    raw_genes = (
        [str(g) for g in adata.var["gene_name"].tolist()]
        if "gene_name" in adata.var.columns
        else [str(g) for g in adata.var_names.tolist()]
    )
    keep_idx: list[int] = []
    dedup_genes: list[str] = []
    seen_genes: set[str] = set()
    for g_i, g_raw in enumerate(raw_genes):
        g_clean = g_raw.strip().upper()
        if g_clean and g_clean not in {"NAN", "NONE", "<NA>"} and g_clean not in seen_genes:
            seen_genes.add(g_clean)
            dedup_genes.append(g_clean)
            keep_idx.append(g_i)

    log2_z_mat = linearize_and_standardize_bulk_matrix(x_valid[:, keep_idx])
    gene_to_idx = {g: idx for idx, g in enumerate(dedup_genes)}

    return Success(
        BulkIcbCohortBundle(
            cohort_id=cohort_id,
            cancer_type=cancer_type,
            short_name=short_name,
            n_samples=n_valid,
            n_responders=int(np.sum(y_valid == 1.0)),
            response_binary=y_valid,
            genes=tuple(dedup_genes),
            gene_to_idx=gene_to_idx,
            log2_z_mat=log2_z_mat,
        )
    )


def load_bulk_icb_suite_via_api() -> Result[tuple[BulkIcbCohortBundle, ...], str]:
    """Load all 11 clinical bulk RNA-seq ICB cohorts (`N = 1,058`) strictly via `tme_datasets`."""
    preproc_dir = get_data_paths().preprocessed_dir
    cohorts: list[BulkIcbCohortBundle] = []
    for cid, cancer_type, short_name in BULK_ICB_COHORT_SPECS:
        match _load_single_bulk_icb_cohort(
            cohort_id=cid,
            cancer_type=cancer_type,
            short_name=short_name,
            preprocessed_dir=preproc_dir,
        ):
            case Failure(err):
                return Failure(err)
            case Success(bundle):
                cohorts.append(bundle)
    return Success(tuple(cohorts))


def persist_comparison_outputs(
    outputs: ComparisonOutputs,
    config: ComparisonConfig,
    exporter: VectorChartExporter = default_vector_exporter,
) -> Result[tuple[Path, ...], str]:
    """Write all Polars Parquet deliverables, Markdown report, and SVG/PNG figures to disk."""
    out_dir = _resolve_repo_path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    parquet_specs: tuple[tuple[str, pl.DataFrame], ...] = (
        ("coulton_substate_centroids.parquet", outputs.coulton_centroids_df),
        ("coulton_top_markers.parquet", outputs.coulton_top_markers_df),
        ("meta42_vs_coulton_full_similarity.parquet", outputs.full_similarity_df),
        ("meta42_macrophage_nearest_coulton_substates.parquet", outputs.nearest_substates_df),
        ("meta42_coulton_sensitivity_summary.parquet", outputs.sensitivity_summary_df),
        ("coulton_substate_ici_outcomes.parquet", outputs.coulton_ici_outcomes_df),
        ("meta42_macrophage_indication_effects.parquet", outputs.meta42_indication_effects_df),
    )

    written: list[Path] = []
    for fname, df in parquet_specs:
        dest = out_dir / fname
        df.write_parquet(dest)
        written.append(dest)

    report_path = out_dir / "meta42_coulton_tam_comparison_report.md"
    report_path.write_text(outputs.markdown_report, encoding="utf-8")
    written.append(report_path)

    if config.render_plots:
        match export_all_comparison_figures(
            full_sim_df=outputs.full_similarity_df,
            nearest_df=outputs.nearest_substates_df,
            coulton_ici_df=outputs.coulton_ici_outcomes_df,
            figures_dir=out_dir / "figures",
            exporter=exporter,
        ):
            case Failure(fig_err):
                return Failure(fig_err)
            case Success(fig_pairs):
                for svg_p, png_p in fig_pairs:
                    written.extend((svg_p, png_p))
        match export_main_paper_typst_caption(
            full_sim_df=outputs.full_similarity_df,
            nearest_df=outputs.nearest_substates_df,
            coulton_ici_df=outputs.coulton_ici_outcomes_df,
            results_dir=out_dir,
        ):
            case Failure(cap_err):
                return Failure(cap_err)
            case Success(cap_p):
                written.append(cap_p)

    return Success(tuple(written))


def execute_comparison_pipeline(
    config: ComparisonConfig,
    streamer: CoultonAtlasStreamer = stream_coulton_centroids_via_api,
    meta42_loader: Meta42ArtifactLoader = load_meta42_artifacts_from_disk,
    bulk_loader: BulkIcbSuiteLoader = load_bulk_icb_suite_via_api,
    exporter: VectorChartExporter = default_vector_exporter,
) -> Result[tuple[ComparisonOutputs, tuple[Path, ...]], str]:
    """Dependency-injected pipeline coordinating IO boundaries and pure core analysis."""
    match meta42_loader(config):
        case Failure(m_err):
            return Failure(m_err)
        case Success(loaded_meta42):
            meta42_bundle = loaded_meta42
        case _:
            return Failure("Unexpected Result state in meta42_loader")

    match load_or_compute_coulton_centroids(config=config, streamer=streamer):
        case Failure(c_err):
            return Failure(c_err)
        case Success(loaded_coulton):
            coulton_bundles = loaded_coulton
        case _:
            return Failure("Unexpected Result state in load_or_compute_coulton_centroids")

    match bulk_loader():
        case Failure(b_err):
            return Failure(b_err)
        case Success(loaded_bulk):
            bulk_cohorts = loaded_bulk
        case _:
            return Failure("Unexpected Result state in bulk_loader")

    match run_comparison_core_pure(
        meta42_bundle=meta42_bundle,
        coulton_bundles=coulton_bundles,
        bulk_cohorts=bulk_cohorts,
        config=config,
    ):
        case Failure(core_err):
            return Failure(core_err)
        case Success(computed_outputs):
            outputs = computed_outputs
        case _:
            return Failure("Unexpected Result state in run_comparison_core_pure")

    match persist_comparison_outputs(outputs=outputs, config=config, exporter=exporter):
        case Failure(p_err):
            return Failure(p_err)
        case Success(written_paths):
            return Success((outputs, written_paths))
        case _:
            return Failure("Unexpected Result state in persist_comparison_outputs")


def parse_cli_config(argv: Sequence[str]) -> ComparisonConfig:
    """Parse CLI flags into an immutable `ComparisonConfig`."""
    parser = argparse.ArgumentParser(
        description=(
            "Compare Consensus-Meta42 macrophage/TME states against Coulton_TAM_Atlas "
            "substates and quantify ICI checkpoint therapy outcomes."
        )
    )
    parser.add_argument(
        "--manuscript-out-dir",
        type=Path,
        default=Path(os.environ.get("TME_PIPELINE_OUTPUT_ROOT", "output/manuscript_pipeline")),
        help="Root directory containing 3b and 3c Consensus-Meta42 upstream parquets.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/results/coulton_tam_comparison"),
        help="Destination directory for Parquet tables, Markdown report, and SVG/PNG figures.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=25_000,
        help="Cell chunk size when streaming Coulton_TAM_Atlas sparse CSR matrices.",
    )
    parser.add_argument(
        "--recompute-centroids",
        action="store_true",
        help="Force recomputation of Coulton_TAM_Atlas centroids even if cached parquet exists.",
    )
    parser.add_argument(
        "--rna-only",
        action="store_true",
        help="Skip SCT@data sensitivity streaming and evaluate RNA@data only.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip Altair SVG/PNG rendering.",
    )
    parsed = parser.parse_args(list(argv))
    ms_root: Path = parsed.manuscript_out_dir
    dir_3b = ms_root / "3b_pan_cohort_unified_deconv_benchmark"
    dir_3c = ms_root / "3c_iatlas_instaprism_milopy_benchmark"
    return ComparisonConfig(
        phi_meta_path=dir_3b / "unified_reference_phi_meta.parquet",
        concordance_3b_path=dir_3b / "unified_all_states_detailed_concordance.parquet",
        hub_summary_path=dir_3c / "consensus_metastate_all42_hub_summary.parquet",
        state_effects_path=dir_3c / "milopy_iatlas_state_effects.parquet",
        top_degs_path=dir_3c / "consensus_metastate_all42_top_degs.parquet",
        output_dir=parsed.output_dir,
        batch_size=parsed.batch_size,
        include_sct_sensitivity=not bool(parsed.rna_only),
        reuse_cached_centroids=not bool(parsed.recompute_centroids),
        render_plots=not bool(parsed.skip_plots),
    )


def main(argv: Maybe[Sequence[str]] = Nothing) -> int:
    """Imperative CLI entrypoint."""
    raw_argv: Sequence[str] = argv.unwrap() if isinstance(argv, Some) else sys.argv[1:]
    config = parse_cli_config(raw_argv)

    match execute_comparison_pipeline(config=config):
        case Failure(err):
            print(f"[ERROR] Comparison pipeline failed: {err}", file=sys.stderr)
            return 1
        case Success((outputs, written_paths)):
            n_near = outputs.nearest_substates_df.height
            print(
                f"[OK] Completed Meta42 <-> Coulton_TAM_Atlas comparison "
                f"({n_near} ranked Top-3 state-substate pairs; {len(written_paths)} files written)."
            )
            for p in written_paths:
                print(f"  -> {p}")
            return 0
        case _:
            return 1


if __name__ == "__main__":
    raise SystemExit(main())
