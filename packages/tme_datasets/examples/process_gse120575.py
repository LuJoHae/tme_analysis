"""End-to-end downstream processing pipeline for Sade-Feldman et al. (GSE120575).

Demonstrates:
1. Loading GSE120575 and caching to H5AD for instantaneous future loading.
2. Clinical metadata standardization (response, timepoint, patient ID, therapy).
3. Confounding gene filtering (MT-, RPS-, RPL-, HLA-).
4. Variance-stabilizing Analytic Pearson Residuals (SCTransform).
5. TME major lineage scoring with Polars.
6. In-silico immune checkpoint target knockout (PDCD1, CTLA4).
7. Dirichlet pseudobulk mixture simulation with known ground-truth fractions.

Run with:
    uv run python packages/tme_datasets/examples/process_gse120575.py
"""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import polars as pl
from returns.maybe import Some
from returns.result import Failure, Success

from tme_datasets import (
    PerturbationConfig,
    PseudobulkConfig,
    SCTransformConfig,
    SCTransformFlavor,
    filter_confounding_genes,
    get_tme_major_lineage_collection,
    in_silico_knockout,
    load_dataset,
    normalize_sctransform,
    score_geneset_zscore,
    simulate_pseudobulk,
)


def harmonize_gse120575_metadata(adata: ad.AnnData) -> ad.AnnData:
    """Standardize raw GEO metadata columns into clean clinical fields."""
    new_adata = adata.copy()

    # 1. Binarize response (Responder: 1, Non-responder: 0)
    if "characteristics: response" in new_adata.obs.columns:
        resp_map = {"Responder": 1, "Non-responder": 0}
        new_adata.obs["response_binary"] = new_adata.obs["characteristics: response"].map(resp_map)

    # 2. Extract patient ID and timepoint (Pre vs Post)
    pat_col = "characteristics: patinet ID (Pre=baseline; Post= on treatment)"
    if pat_col in new_adata.obs.columns:
        pat_series = new_adata.obs[pat_col].astype(str)
        new_adata.obs["timepoint"] = pat_series.apply(
            lambda s: "Pre" if "Pre" in s else ("Post" if "Post" in s else "Unknown")
        )
        new_adata.obs["patient_id"] = pat_series.str.extract(r"(P\d+)")

    # 3. Therapy
    if "characteristics: therapy" in new_adata.obs.columns:
        new_adata.obs["therapy"] = new_adata.obs["characteristics: therapy"]

    return new_adata


def main() -> None:
    print("=" * 70)
    print("GSE120575 Sade-Feldman et al. (Cell 2018) Single-Cell Pipeline")
    print("=" * 70)

    cache_file = Path("data/preprocessed/GSE120575.h5ad")

    # Step 1: Ingestion / Cache Check
    if cache_file.exists():
        print(f"\n[1/7] Loading cached AnnData from {cache_file}...")
        adata = ad.read_h5ad(cache_file)
    else:
        print("\n[1/7] Loading raw GSE120575 via tme_datasets...")
        match load_dataset("GSE120575"):
            case Success(loaded_adata):
                adata = loaded_adata
                cache_file.parent.mkdir(parents=True, exist_ok=True)
                adata.obs_names.name = "cell_id"
                adata.var_names.name = "gene_id"
                adata.write_h5ad(cache_file)
                print(f"Saved cached dataset to {cache_file}")
            case Failure(err):
                print(f"Failed to load GSE120575: {err}")
                return

    print(f"Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Step 2: Clinical Feature Harmonization
    print("\n[2/7] Harmonizing clinical observations...")
    adata = harmonize_gse120575_metadata(adata)
    if "response_binary" in adata.obs:
        resp_counts = adata.obs["response_binary"].value_counts().to_dict()
        print(f"Response distribution (1=Responder, 0=Non-responder): {resp_counts}")
    if "timepoint" in adata.obs:
        print(f"Timepoint distribution: {adata.obs['timepoint'].value_counts().to_dict()}")

    # Step 3: Confounding Gene Filtering
    print("\n[3/7] Filtering technical and confounding genes (MT-, RPS-, RPL-, HLA-)...")
    clean_adata = filter_confounding_genes(adata).unwrap()
    print(f"Filtered matrix shape: {clean_adata.shape} (dropped {adata.n_vars - clean_adata.n_vars} genes)")

    # Step 4: Variance Stabilization (SCTransform Analytic Pearson Residuals)
    print("\n[4/7] Normalizing via Analytic Pearson Residuals (SCTransform)...")
    # Work on a subsample if running on interactive dev workstation for speed
    work_adata = clean_adata if clean_adata.n_obs <= 2000 else clean_adata[:2000].copy()
    sct_cfg = SCTransformConfig(
        flavor=SCTransformFlavor.ANALYTIC,
        n_top_genes=Some(2000),
        clip_residuals=True,
        use_layer_as_x=False,
    )
    norm_adata = normalize_sctransform(work_adata, sct_cfg).unwrap()
    print(f"Calculated Pearson residuals layer: {norm_adata.layers['pearson_residuals'].shape}")
    print(f"Selected HVGs: {norm_adata.var['highly_variable'].sum()}")

    # Step 5: TME Major Lineage Scoring
    print("\n[5/7] Scoring major tumor microenvironment lineage signatures...")
    lineage_coll = get_tme_major_lineage_collection()
    print(f"Lineage gene sets: {list(lineage_coll.gene_sets.keys())}")
    scores_df = score_geneset_zscore(norm_adata, lineage_coll).unwrap()
    print(f"Lineage scores DataFrame (Polars shape): {scores_df.shape}")
    print(scores_df.head(3))

    # Assign coarse pseudo cell-type based on top lineage score
    numeric_cols = [c for c in scores_df.columns if c not in ("cell_id", "sample_id")]
    matrix_scores = scores_df.select(numeric_cols).to_numpy().astype(np.float64)
    max_idx = np.argmax(matrix_scores, axis=1)
    norm_adata.obs["cell_type"] = [numeric_cols[i] for i in max_idx]
    if norm_adata.obs["cell_type"].nunique() < 2 and "timepoint" in norm_adata.obs:
        norm_adata.obs["cell_type"] = norm_adata.obs["timepoint"]
    print(f"Cell-type assignments preview:\n{norm_adata.obs['cell_type'].value_counts()}")

    # Step 6: In-Silico Immune Checkpoint Knockout
    print("\n[6/7] Simulating in-silico knockout of immune checkpoints (PDCD1, CTLA4)...")
    ko_cfg = PerturbationConfig(target_genes=("PDCD1", "CTLA4"))
    ko_adata = in_silico_knockout(norm_adata, ko_cfg).unwrap()
    print(f"Knockout perturbation complete on {ko_adata.shape[0]} cells.")

    # Step 7: Pseudobulk Deconvolution Simulation
    print("\n[7/7] Generating synthetic pseudobulk mixtures with known cell fractions...")
    sim_cfg = PseudobulkConfig(
        n_samples=10,
        cells_per_sample=200,
        seed=Some(42),
    )
    sim_bulk, truth_df = simulate_pseudobulk(
        norm_adata,
        sim_cfg,
        cell_type_key="cell_type",
    ).unwrap()

    print(f"Simulated bulk AnnData: {sim_bulk.shape} (samples x genes)")
    print("Ground-truth cell-type proportions:")
    print(truth_df.head(3))

    print("\n" + "=" * 70)
    print("GSE120575 downstream pipeline completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
