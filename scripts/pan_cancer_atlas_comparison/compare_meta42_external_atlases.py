#!/usr/bin/env python3
"""Main CLI & imperative shell comparing Consensus-Meta42 to Pan-Cancer T-Cell, CAF & TAM Atlases.

Architecturally follows Functional Core, Imperative Shell (.agents/rules/code-style-guide.md):
- Loads all 6 T-cell & CAF atlases (`Gao_CAF_Atlas`, `Cords_CAF_Atlas`, `Buechler_CAF_Atlas`,
  `Zheng_TCell_Atlas`, `Chu_TCell_Atlas`, `Shi_PanCancer_Atlas`), `Coulton_TAM_Atlas`, and the
  11 clinical bulk ICB cohorts (`N = 1,058`) strictly via the `tme_datasets` public API.
- Streams sparse CSR batches (`25,000` cells/chunk) into pure out-of-core accumulators so peak
  RAM remains `< 4.0 GB` on 16 GB machines.
- Delegates all similarity algebra, cross-atlas equivalence, DerSimonian-Laird bulk ICB
  meta-analysis, conditional stromal gating, and reporting to `scripts.pan_cancer_atlas_comparison.core`.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import gc
from pathlib import Path
import sys
from typing import Final, assert_never

import anndata as ad  # type: ignore[import-untyped]
import numpy as np
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp

_REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tme_datasets import (  # noqa: E402
    ShiCompartmentFilter,
    get_data_paths,
    load_buechler_caf_atlas,
    load_chu_tcell_atlas,
    load_cords_caf_atlas,
    load_coulton_tam_atlas,
    load_dataset,
    load_gao_caf_atlas,
    load_shi_pancancer_atlas,
    load_zheng_tcell_atlas,
)

from scripts.pan_cancer_atlas_comparison.core import (  # noqa: E402
    accumulate_atlas_batch,
    adapt_coulton_cached_parquet_to_bundle,
    build_chu_curated_centroid_bundle,
    create_empty_atlas_accumulator,
    deduplicate_atlas_gene_columns,
    deserialize_atlas_centroids_from_polars,
    filter_atlas_min_cells,
    finalize_atlas_centroids,
    run_pan_atlas_comparison_pure,
    serialize_atlas_centroids_to_polars,
)
from scripts.pan_cancer_atlas_comparison.models import (  # noqa: E402
    AtlasCentroidBundle,
    BULK_ICB_COHORT_SPECS,
    BulkIcbCohortBundle,
    ExternalAtlasId,
    FloatMatrix,
    Meta42ReferenceBundle,
    PanAtlasConfig,
    ProfileOrigin,
    VectorChartExporter,
)
from scripts.pan_cancer_atlas_comparison.render_figures import (  # noqa: E402
    default_vector_exporter,
    export_all_pan_atlas_figures,
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


def _extract_adata_gene_symbols(adata: ad.AnnData) -> tuple[str, ...]:
    """Resolve canonical HGNC symbols from `.var['gene_name']` with `.var_names` fallback."""
    var_cols = set(adata.var.columns)
    primary = (
        [str(x).strip() for x in adata.var["gene_name"].tolist()]
        if "gene_name" in var_cols
        else [str(x).strip() for x in adata.var_names.tolist()]
    )
    fallback = (
        [str(x).strip() for x in adata.var["original_id"].tolist()]
        if "original_id" in var_cols
        else [str(x).strip() for x in adata.var_names.tolist()]
    )
    return tuple(
        gn if (gn and gn.upper() not in {"NAN", "NONE", "<NA>", "NA"}) else fb
        for gn, fb in zip(primary, fallback, strict=False)
    )


def _infer_is_tumor_mask(adata: ad.AnnData) -> np.ndarray[tuple[int], np.dtype[np.bool_]]:
    """Determine tumor-tissue cell membership across harmonized `.obs` columns."""
    obs_cols = set(adata.obs.columns)
    if "sample_type" in obs_cols:
        vals = [str(v).strip().lower() for v in adata.obs["sample_type"].tolist()]
        mask = np.asarray(
            [v in {"tumor", "metastasis", "primary", "tumour", "t", "m"} for v in vals],
            dtype=np.bool_,
        )
        if bool(np.any(mask)):
            return mask
    if "tissue" in obs_cols:
        t_vals = [str(v).strip().lower() for v in adata.obs["tissue"].tolist()]
        t_mask = np.asarray(
            [v in {"tumor", "tumour", "metastasis"} for v in t_vals],
            dtype=np.bool_,
        )
        if bool(np.any(t_mask)):
            return t_mask
    return np.ones(int(adata.n_obs), dtype=np.bool_)


def _stream_adata_into_bundle(
    adata: ad.AnnData,
    atlas_id: ExternalAtlasId,
    state_col_candidates: Sequence[str],
    row_filter_mask: Maybe[np.ndarray[tuple[int], np.dtype[np.bool_]]],
    batch_size: int,
    min_cells: int,
) -> Result[AtlasCentroidBundle, str]:
    """Stream sparse CSR chunks from a loaded `AnnData` into an `AtlasCentroidBundle`."""
    obs_cols = set(adata.obs.columns)
    chosen_col = next((c for c in state_col_candidates if c in obs_cols), "")
    if not chosen_col:
        return Failure(
            f"{atlas_id.value} .obs lacks expected cluster column among {list(state_col_candidates)}"
        )

    raw_states = [str(x).strip() for x in adata.obs[chosen_col].tolist()]
    row_keep = (
        row_filter_mask.unwrap()
        if isinstance(row_filter_mask, Some)
        else np.ones(int(adata.n_obs), dtype=np.bool_)
    )

    valid_labels: list[str] = []
    seen_labels: set[str] = set()
    for st_val, keep_flag in zip(raw_states, row_keep.tolist(), strict=False):
        if (
            keep_flag
            and st_val
            and st_val.lower() not in {"na", "nan", "none", "unassigned", "unknown", "23_na"}
            and st_val not in seen_labels
        ):
            seen_labels.add(st_val)
            valid_labels.append(st_val)

    ordered_labels = tuple(sorted(valid_labels))
    if not ordered_labels:
        return Failure(f"No valid cluster labels found in {atlas_id.value} ({chosen_col})")

    label_to_code = {lbl: idx for idx, lbl in enumerate(ordered_labels)}
    cluster_codes = np.asarray(
        [
            label_to_code.get(st_val, -1) if keep_flag else -1
            for st_val, keep_flag in zip(raw_states, row_keep.tolist(), strict=False)
        ],
        dtype=np.int64,
    )
    tumor_mask = _infer_is_tumor_mask(adata)
    raw_symbols = _extract_adata_gene_symbols(adata)

    acc = create_empty_atlas_accumulator(
        atlas_id=atlas_id,
        origin=ProfileOrigin.EMPIRICAL_SINGLE_CELL_H5AD,
        cluster_labels=ordered_labels,
        gene_symbols=raw_symbols,
    )

    n_cells = int(adata.n_obs)
    for start in range(0, n_cells, batch_size):
        stop = min(start + batch_size, n_cells)
        codes_slice = cluster_codes[start:stop]
        if not bool(np.any(codes_slice >= 0)):
            continue
        csr_batch = sp.csr_matrix(adata.X[start:stop])
        match accumulate_atlas_batch(
            acc=acc,
            batch_csr=csr_batch,
            cluster_codes=codes_slice,
            is_tumor_mask=tumor_mask[start:stop],
        ):
            case Failure(err):
                return Failure(err)
            case Success(next_acc):
                acc = next_acc

    finalized = finalize_atlas_centroids(acc)
    deduped = deduplicate_atlas_gene_columns(raw_symbols=raw_symbols, bundle=finalized)
    return Success(filter_atlas_min_cells(deduped, min_cells=min_cells))


def stream_single_external_atlas(
    atlas_id: ExternalAtlasId,
    config: PanAtlasConfig,
    fallback_gene_universe: Sequence[str],
) -> Result[AtlasCentroidBundle, str]:
    """Load and stream one external atlas strictly via `tme_datasets` public API."""
    match atlas_id:
        case ExternalAtlasId.GAO_CAF:
            match load_gao_caf_atlas(auto_download=False, include_all_layers=False):
                case Failure(err):
                    return Failure(err)
                case Success(ad_gao):
                    return _stream_adata_into_bundle(
                        ad_gao,
                        atlas_id,
                        ("author_cluster_id", "Cluster", "cell_state"),
                        Nothing,
                        config.batch_size,
                        config.min_cells_per_substate,
                    )
                case _:
                    return Failure("Unexpected Result variant from load_gao_caf_atlas")

        case ExternalAtlasId.CORDS_CAF:
            match load_cords_caf_atlas(auto_download=False, include_all_layers=False):
                case Failure(err):
                    return Failure(err)
                case Success(ad_cords):
                    return _stream_adata_into_bundle(
                        ad_cords,
                        atlas_id,
                        ("author_cluster_id", "cell_state", "cell_subtype"),
                        Nothing,
                        config.batch_size,
                        config.min_cells_per_substate,
                    )
                case _:
                    return Failure("Unexpected Result variant from load_cords_caf_atlas")

        case ExternalAtlasId.BUECHLER_CAF:
            match load_buechler_caf_atlas(auto_download=False, include_all_layers=False):
                case Failure(err):
                    return Failure(err)
                case Success(ad_buechler):
                    return _stream_adata_into_bundle(
                        ad_buechler,
                        atlas_id,
                        ("author_cluster_id", "Cluster", "cell_state"),
                        Nothing,
                        config.batch_size,
                        config.min_cells_per_substate,
                    )
                case _:
                    return Failure("Unexpected Result variant from load_buechler_caf_atlas")

        case ExternalAtlasId.ZHENG_TCELL:
            match load_zheng_tcell_atlas(auto_download=False, include_all_layers=False):
                case Failure(err):
                    return Failure(err)
                case Success(ad_zheng):
                    return _stream_adata_into_bundle(
                        ad_zheng,
                        atlas_id,
                        ("author_cluster_id", "cell_state", "meta.cluster"),
                        Nothing,
                        config.batch_size,
                        config.min_cells_per_substate,
                    )
                case _:
                    return Failure("Unexpected Result variant from load_zheng_tcell_atlas")

        case ExternalAtlasId.CHU_TCELL:
            match load_chu_tcell_atlas(auto_download=False, include_all_layers=False):
                case Success(ad_chu):
                    return _stream_adata_into_bundle(
                        ad_chu,
                        atlas_id,
                        ("author_cluster_id", "cell_state", "cell_subtype"),
                        Nothing,
                        config.batch_size,
                        config.min_cells_per_substate,
                    )
                case Failure(_):
                    print(
                        "[INFO] Chu_TCell_Atlas raw .rds not staged locally; constructing "
                        "ProfileOrigin.CURATED_SIGNATURE_PRIOR bundle across 39 Chu T-cell states."
                    )
                    return Success(build_chu_curated_centroid_bundle(fallback_gene_universe))
                case _:
                    return Failure("Unexpected Result variant from load_chu_tcell_atlas")

        case ExternalAtlasId.SHI_TCELL | ExternalAtlasId.SHI_CAF:
            match load_shi_pancancer_atlas(
                compartment_filter=ShiCompartmentFilter.T_AND_CAF,
                auto_download=False,
                include_all_layers=False,
            ):
                case Failure(err):
                    return Failure(err)
                case Success(ad_shi):
                    comp_vals = (
                        [str(x).strip() for x in ad_shi.obs["lineage_compartment"].tolist()]
                        if "lineage_compartment" in ad_shi.obs.columns
                        else ["T_CELL"] * int(ad_shi.n_obs)
                    )
                    target_set = (
                        {"T_CELL", "NK_ILC"}
                        if atlas_id == ExternalAtlasId.SHI_TCELL
                        else {"STROMAL_CAF"}
                    )
                    mask_shi = np.asarray([c in target_set for c in comp_vals], dtype=np.bool_)
                    return _stream_adata_into_bundle(
                        ad_shi,
                        atlas_id,
                        ("cell_state", "author_cluster_id", "cell_subtype"),
                        Some(mask_shi),
                        config.batch_size,
                        config.min_cells_per_substate,
                    )
                case _:
                    return Failure("Unexpected Result variant from load_shi_pancancer_atlas")

        case ExternalAtlasId.COULTON_TAM:
            coulton_pq = _resolve_repo_path(config.coulton_centroids_path)
            if coulton_pq.is_file():
                return adapt_coulton_cached_parquet_to_bundle(pl.read_parquet(coulton_pq))
            match load_coulton_tam_atlas(assay="RNA", auto_download=False, include_all_layers=False):
                case Failure(err):
                    return Failure(err)
                case Success(ad_tam):
                    return _stream_adata_into_bundle(
                        ad_tam,
                        atlas_id,
                        ("cell_state", "short.label"),
                        Nothing,
                        config.batch_size,
                        config.min_cells_per_substate,
                    )
                case _:
                    return Failure("Unexpected Result variant from load_coulton_tam_atlas")

        case _ as unreachable:
            assert_never(unreachable)


def load_or_compute_all_atlas_centroids(
    config: PanAtlasConfig,
    meta42_bundle: Meta42ReferenceBundle,
) -> Result[tuple[AtlasCentroidBundle, ...], str]:
    """Load cached external atlas centroids or stream each available atlas via `tme_datasets`."""
    out_dir = _resolve_repo_path(config.output_dir)
    cache_pq = out_dir / "external_atlas_substate_centroids.parquet"
    target_order: tuple[ExternalAtlasId, ...] = (
        ExternalAtlasId.ZHENG_TCELL,
        ExternalAtlasId.CHU_TCELL,
        ExternalAtlasId.SHI_TCELL,
        ExternalAtlasId.GAO_CAF,
        ExternalAtlasId.CORDS_CAF,
        ExternalAtlasId.BUECHLER_CAF,
        ExternalAtlasId.SHI_CAF,
        ExternalAtlasId.COULTON_TAM,
    )

    cached_map: dict[ExternalAtlasId, AtlasCentroidBundle] = {}
    if config.reuse_cached_centroids and cache_pq.is_file():
        match deserialize_atlas_centroids_from_polars(pl.read_parquet(cache_pq)):
            case Success(existing_bundles):
                for eb in existing_bundles:
                    cached_map[eb.atlas_id] = eb
            case Failure(_):
                pass

    gene_universe = tuple(c for c in meta42_bundle.phi_meta_df.columns if c != "cluster")
    bundles: list[AtlasCentroidBundle] = []

    for aid in target_order:
        if aid in cached_map and (
            cached_map[aid].origin == ProfileOrigin.EMPIRICAL_SINGLE_CELL_H5AD
            or aid == ExternalAtlasId.CHU_TCELL
        ):
            print(f"[INFO] Reusing cached centroids for {aid.value} ({len(cached_map[aid].cluster_labels)} clusters)")
            bundles.append(cached_map[aid])
            continue

        print(f"[INFO] Streaming {aid.value} via tme_datasets ({config.batch_size:,} cells/batch)...")
        match stream_single_external_atlas(aid, config=config, fallback_gene_universe=gene_universe):
            case Success(bundle):
                print(
                    f"[OK] {aid.value}: {len(bundle.cluster_labels)} clusters, "
                    f"{int(bundle.all_cell_counts.sum()):,} cells, {len(bundle.gene_symbols):,} genes"
                )
                bundles.append(bundle)
            case Failure(err):
                print(f"[WARN] Skipping {aid.value}: {err}", file=sys.stderr)
        gc.collect()

    if not bundles:
        return Failure("No external reference atlases could be loaded or constructed")

    out_dir.mkdir(parents=True, exist_ok=True)
    serialize_atlas_centroids_to_polars(bundles).write_parquet(cache_pq)
    return Success(tuple(bundles))


def _resolve_artifact_with_fallback(primary: Path, fallback_dir: Path) -> Path:
    """Resolve primary artifact path with fallback to `output/results/iatlas_instaprism_benchmark`."""
    resolved = _resolve_repo_path(primary)
    if resolved.is_file():
        return resolved
    alt = _resolve_repo_path(fallback_dir / primary.name)
    return alt if alt.is_file() else resolved


def load_meta42_artifacts_from_disk(config: PanAtlasConfig) -> Result[Meta42ReferenceBundle, str]:
    """Load Consensus-Meta42 reference profiles and Stage 3b/3c ICI effect parquets."""
    legacy_3c_dir = Path("output/results/iatlas_instaprism_benchmark")
    phi_p = _resolve_repo_path(config.phi_meta_path)
    conc_p = _resolve_repo_path(config.concordance_3b_path)
    hub_p = _resolve_artifact_with_fallback(config.hub_summary_path, legacy_3c_dir)
    eff_p = _resolve_artifact_with_fallback(config.state_effects_path, legacy_3c_dir)
    deg_p = _resolve_artifact_with_fallback(config.top_degs_path, legacy_3c_dir)

    for label, path in (
        ("phi_meta", phi_p),
        ("concordance_3b", conc_p),
        ("hub_summary", hub_p),
        ("state_effects", eff_p),
        ("top_degs", deg_p),
    ):
        if not path.is_file():
            return Failure(f"Missing required Consensus-Meta42 artifact ({label}): {path}")

    return Success(
        Meta42ReferenceBundle(
            phi_meta_df=pl.read_parquet(phi_p),
            concordance_3b_df=pl.read_parquet(conc_p),
            hub_summary_df=pl.read_parquet(hub_p),
            state_effects_df=pl.read_parquet(eff_p),
            top_degs_df=pl.read_parquet(deg_p),
        )
    )


def load_bulk_icb_suite_via_api() -> Result[tuple[BulkIcbCohortBundle, ...], str]:
    """Load all 11 clinical bulk RNA-seq ICB cohorts (`N = 1,058`) strictly via `tme_datasets.load_dataset`."""
    prep_dir = get_data_paths().preprocessed_dir
    loaded: list[BulkIcbCohortBundle] = []

    for cohort_id, cancer_type, short_name in BULK_ICB_COHORT_SPECS:
        match load_dataset(
            cohort_id,
            base_dir=prep_dir,
            auto_download=False,
            normalize_ensembl=False,
        ):
            case Failure(err):
                return Failure(f"Failed to load bulk ICB cohort '{cohort_id}': {err}")
            case Success(adata):
                pass

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
        genes = _extract_adata_gene_symbols(adata)
        gene_to_idx: dict[str, int] = {}
        for idx, g_sym in enumerate(genes):
            if g_sym and g_sym not in gene_to_idx:
                gene_to_idx[g_sym] = idx

        loaded.append(
            BulkIcbCohortBundle(
                cohort_id=cohort_id,
                cancer_type=cancer_type,
                short_name=short_name,
                n_samples=n_valid,
                n_responders=int(np.sum(y_valid == 1.0)),
                response_binary=y_valid,
                genes=genes,
                gene_to_idx=gene_to_idx,
                log2_z_mat=linearize_and_standardize_bulk_matrix(x_valid),
            )
        )

    return Success(tuple(loaded))


def parse_args(argv: Sequence[str]) -> PanAtlasConfig:
    """Parse CLI arguments into a frozen `PanAtlasConfig`."""
    parser = argparse.ArgumentParser(
        description="Validate Consensus-Meta42 against Pan-Cancer T-Cell, CAF, and Myeloid Triad Atlases."
    )
    parser.add_argument(
        "--manuscript-out-dir",
        type=Path,
        default=Path("output/manuscript_pipeline"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/results/pan_cancer_atlas_comparison"),
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=Path("output/figures/pan_cancer_atlas_comparison"),
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=25_000,
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Re-stream all external atlas centroids from .h5ad even if cached parquet exists.",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Compute parquet deliverables and Markdown report without rendering Altair figures.",
    )
    ns = parser.parse_args(list(argv))
    ms_root: Path = ns.manuscript_out_dir
    return PanAtlasConfig(
        output_dir=ns.output_dir,
        figures_dir=ns.figures_dir,
        phi_meta_path=ms_root / "3b_pan_cohort_unified_deconv_benchmark/unified_reference_phi_meta.parquet",
        concordance_3b_path=ms_root / "3b_pan_cohort_unified_deconv_benchmark/unified_all_states_detailed_concordance.parquet",
        hub_summary_path=ms_root / "3c_iatlas_instaprism_milopy_benchmark/consensus_metastate_all42_hub_summary.parquet",
        state_effects_path=ms_root / "3c_iatlas_instaprism_milopy_benchmark/milopy_iatlas_state_effects.parquet",
        top_degs_path=ms_root / "3c_iatlas_instaprism_milopy_benchmark/consensus_metastate_all42_top_degs.parquet",
        batch_size=int(ns.batch_size),
        reuse_cached_centroids=not bool(ns.no_cache),
        render_plots=not bool(ns.skip_plots),
    )


def run_cli(
    config: PanAtlasConfig,
    exporter: VectorChartExporter = default_vector_exporter,
) -> Result[Path, str]:
    """Execute full Imperative Shell -> Pure Functional Core -> Parquet/SVG persistence."""
    match load_meta42_artifacts_from_disk(config):
        case Failure(err):
            return Failure(err)
        case Success(meta42_bundle):
            pass

    match load_or_compute_all_atlas_centroids(config, meta42_bundle=meta42_bundle):
        case Failure(err):
            return Failure(err)
        case Success(atlas_bundles):
            pass

    match load_bulk_icb_suite_via_api():
        case Failure(err):
            return Failure(err)
        case Success(bulk_cohorts):
            pass

    print(
        f"[INFO] Running pure pan-atlas Triad synthesis across {len(atlas_bundles)} atlases "
        f"and {sum(c.n_samples for c in bulk_cohorts):,} ICB patients..."
    )
    match run_pan_atlas_comparison_pure(
        meta42_bundle=meta42_bundle,
        atlas_bundles=atlas_bundles,
        bulk_cohorts=bulk_cohorts,
        config=config,
    ):
        case Failure(err):
            return Failure(err)
        case Success(outputs):
            pass

    out_dir = _resolve_repo_path(config.output_dir)
    fig_dir = _resolve_repo_path(config.figures_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs.atlas_centroids_df.write_parquet(out_dir / "external_atlas_substate_centroids.parquet")
    outputs.atlas_top_markers_df.write_parquet(out_dir / "external_atlas_top_markers.parquet")
    outputs.meta42_vs_atlases_similarity_df.write_parquet(
        out_dir / "meta42_vs_external_atlases_similarity.parquet"
    )
    outputs.nearest_tcell_matches_df.write_parquet(out_dir / "meta42_tcell_nearest_substates.parquet")
    outputs.nearest_caf_matches_df.write_parquet(out_dir / "meta42_caf_nearest_substates.parquet")
    outputs.cross_atlas_equivalence_df.write_parquet(
        out_dir / "cross_atlas_substate_equivalence.parquet"
    )
    outputs.atlas_ici_outcomes_df.write_parquet(
        out_dir / "external_atlas_substate_ici_outcomes.parquet"
    )
    outputs.stromal_gating_df.write_parquet(out_dir / "caf_tcell_stromal_gating_icb.parquet")
    outputs.triad_block_summary_df.write_parquet(out_dir / "triad_block_specificity_summary.parquet")

    report_path = out_dir / "external_atlas_comparison_report.md"
    report_path.write_text(outputs.markdown_report, encoding="utf-8")

    if config.render_plots:
        match export_all_pan_atlas_figures(
            full_sim_df=outputs.meta42_vs_atlases_similarity_df,
            nearest_tcell_df=outputs.nearest_tcell_matches_df,
            nearest_caf_df=outputs.nearest_caf_matches_df,
            atlas_ici_df=outputs.atlas_ici_outcomes_df,
            triad_summary_df=outputs.triad_block_summary_df,
            figures_dir=fig_dir,
            exporter=exporter,
        ):
            case Failure(f_err):
                return Failure(f_err)
            case Success(fig_paths):
                for svg_p, png_p in fig_paths:
                    print(f"[OK] Exported {svg_p} & {png_p}")

    return Success(report_path)


def main(argv: Maybe[Sequence[str]] = Nothing) -> int:
    """CLI entrypoint."""
    raw_argv: Sequence[str] = argv.unwrap() if isinstance(argv, Some) else sys.argv[1:]
    cfg = parse_args(raw_argv)
    match run_cli(cfg):
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            return 1
        case Success(report_path):
            print(f"[SUCCESS] Pan-Cancer Triad Atlas comparison complete -> {report_path}")
            return 0
        case _:
            return 1


if __name__ == "__main__":
    raise SystemExit(main())
