"""End-to-End Tutorial Script: Exploring All Features of tme_datasets.

Run with:
    uv run python packages/tme_datasets/examples/tutorial_explore_all_features.py
"""

from __future__ import annotations

import tempfile
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from returns.maybe import Some
from returns.result import Failure, Success

# Core imports from tme_datasets
from tme_datasets import (
    # Core types & configs
    HarmonizeConfig,
    HarmonizeMode,
    NegativeBinomialConfig,
    PseudobulkConfig,
    SubsampleSpec,
    GeneReconcileConfig,
    GeneIDType,
    ChecksumSpec,
    StorageBackend,
    # Query & Registry
    list_registered_datasets,
    get_dataset_spec,
    load_dataset,
    query_datasets,
    load_geneset_collection,
    # Transforms & Perturbations
    subsample_cells,
    supersample_cells,
    randomize_negative_binomial,
    simulate_dropout,
    add_expression_jitter,
    in_silico_knockout,
    in_silico_overexpression,
    ComposeTransforms,
    # Gene sets & Scoring
    get_tme_major_lineage_collection,
    score_geneset_zscore,
    score_geneset_auc,
    compute_geneset_overlap,
    export_gmt,
    parse_gmt,
    # Simulation & Metrics & Harmonization
    simulate_pseudobulk,
    evaluate_integration_metrics,
    align_and_concatenate,
    # PyTorch
    TmeTorchDataset,
    create_tme_dataloader,
    # Storage & Genes & Verification
    convert_to_zarr,
    load_backed,
    reconcile_genes,
    compute_file_hash,
    verify_checksum,
)


def print_section(title: str) -> None:
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)


def create_demo_adata(n_cells: int = 120, n_genes: int = 30) -> ad.AnnData:
    """Helper creating a synthetic AnnData containing authentic TME marker genes."""
    rng = np.random.default_rng(42)
    genes = [f"GENE_{i}" for i in range(n_genes)]
    # Embed canonical TME marker genes
    tme_genes = ["CD8A", "CD8B", "CD4", "CD68", "CD14", "MS4A1", "PDCD1", "CD274", "CTLA4", "EPCAM"]
    for i, g in enumerate(tme_genes):
        genes[i] = g

    counts = rng.negative_binomial(n=4, p=0.25, size=(n_cells, n_genes)).astype(np.float32)
    cell_types = rng.choice(["CD8_T", "CD4_T", "Macrophage", "B_cell", "Tumor"], size=n_cells)
    response = rng.choice([1.0, 0.0], size=n_cells, p=[0.45, 0.55])
    batches = rng.choice(["Study_Alpha", "Study_Beta"], size=n_cells)

    obs = pd.DataFrame(
        {
            "cell_type": cell_types,
            "response_binary": response,
            "batch": batches,
            "sample_id": [f"cell_{i:04d}" for i in range(n_cells)],
        },
        index=[f"cell_{i:04d}" for i in range(n_cells)],
    )

    var = pd.DataFrame(index=genes)
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.obsm["X_pca"] = rng.normal(size=(n_cells, 10))
    return adata


def main() -> None:
    print_section("1. REGISTRY & METADATA INSPECTION")
    specs = list_registered_datasets()
    print(f"Total registered datasets in tme_datasets: {len(specs)}")
    print("\nFirst 6 registered datasets:")
    for s in specs[:6]:
        organ = s.organ.value_or("N/A") if isinstance(s.organ, Some) else "N/A"
        cells = s.n_samples_or_cells.value_or("N/A") if isinstance(s.n_samples_or_cells, Some) else "N/A"
        print(f"  - [{s.modality.value.upper()}] {s.id:<18} | {s.cancer_type:<22} | Organ: {organ:<10} | N: {cells}")

    print_section("2. GENE SET COLLECTIONS & SIGNATURE SCORING")
    coll = get_tme_major_lineage_collection()
    print(f"Loaded GeneSetCollection: '{coll.name}' ({len(coll.gene_sets)} signatures)")
    for name, gs in coll.gene_sets.items():
        print(f"  * {name:<18}: {', '.join(gs.genes[:6])}")

    demo_adata = create_demo_adata()

    # Score signatures using Z-Score and AUCell
    z_res = score_geneset_zscore(demo_adata, coll)
    auc_res = score_geneset_auc(demo_adata, coll)

    if isinstance(z_res, Success) and isinstance(auc_res, Success):
        df_z = z_res.unwrap()
        df_auc = auc_res.unwrap()
        print("\nSignature Mean Z-scores (Top 3 cells):")
        print(df_z.head(3))
        print("\nSignature Rank AUCell Scores (Top 3 cells):")
        print(df_auc.head(3))

    # Overlap analysis between signatures
    overlap_res = compute_geneset_overlap(coll, coll)
    if isinstance(overlap_res, Success):
        df_over = overlap_res.unwrap()
        cross = df_over.filter(df_over["signature_a"] != df_over["signature_b"])
        print("\nCross-Signature Overlap (First 3 pairs):")
        print(cross.select(["signature_a", "signature_b", "n_shared_genes", "jaccard_similarity"]).head(3))

    print_section("3. TRANSFORMS, SUBSAMPLING & PERTURBATIONS")
    print(f"Original AnnData shape: {demo_adata.shape}")

    # A. Stratified Subsampling
    sub_spec = SubsampleSpec(n_or_fraction=50, stratify_by=Some("cell_type"), balanced=True, seed=Some(42))
    sub_adata = subsample_cells(demo_adata, sub_spec).unwrap()
    print(f"Balanced Stratified Subsample (50 cells): {sub_adata.shape}")
    print(f"Cell type distribution:\n{sub_adata.obs['cell_type'].value_counts().to_dict()}")

    # B. Negative Binomial Count Randomization & Parameter Inference
    # Mode 1: Simple entry-wise Gamma-Poisson noise (zeros stay 0)
    nb_cfg = NegativeBinomialConfig(dispersion=0.20, seed=Some(123))
    nb_adata = randomize_negative_binomial(demo_adata, nb_cfg).unwrap()
    print(f"\nSimple Negative Binomial randomized counts preserved in .layers:")
    print(f"  Available layers: {list(nb_adata.layers.keys())}")
    print(f"  Mean count (original): {np.mean(demo_adata.X):.2f} | Mean count (NB): {np.mean(nb_adata.X):.2f}")

    # Mode 2: Parameter Inference Mode (Method of Moments / Empirical Bayes)
    from tme_datasets import NBEstimationMethod
    nb_mom_cfg = NegativeBinomialConfig(
        estimation_method=Some(NBEstimationMethod.MOMENTS),
        cluster_key=Some("cell_type"),
        seed=Some(42),
    )
    mom_adata = randomize_negative_binomial(demo_adata, nb_mom_cfg).unwrap()
    print(f"Parameter-Inferred NB (MoM stratified by cell_type): shape = {mom_adata.shape}")

    # C. In-silico Targeted Knockout & Overexpression
    ko_adata = in_silico_knockout(demo_adata, genes=("PDCD1", "CD274"), efficiency=1.0).unwrap()
    pd1_idx = list(ko_adata.var_names).index("PDCD1")
    print(f"\nIn-silico Knockout of PDCD1 (100% efficiency): Max expression = {np.max(ko_adata.X[:, pd1_idx])}")

    # D. Monadic Transform Pipeline Composition
    pipeline = ComposeTransforms([
        lambda a: subsample_cells(a, SubsampleSpec(n_or_fraction=80, seed=Some(1))),
        lambda a: simulate_dropout(a, rate=0.05, seed=2),
        lambda a: add_expression_jitter(a, sigma=0.05, seed=3),
    ])
    piped_adata = pipeline(demo_adata).unwrap()
    print(f"Chained Pipeline (Subsample -> Dropout -> Jitter): Final shape = {piped_adata.shape}")

    print_section("4. IN-SILICO PSEUDOBULK DECONVOLUTION BENCHMARK GENERATOR")
    pb_cfg = PseudobulkConfig(n_samples=5, cells_per_sample=500, noise_dispersion=Some(0.05), seed=Some(99))
    bulk_adata, truth_df = simulate_pseudobulk(demo_adata, pb_cfg, cell_type_key="cell_type").unwrap()
    print(f"Simulated Bulk RNA-seq Mixtures: {bulk_adata.shape} (Samples x Genes)")
    print("\nExact Known Dirichlet Ground-Truth Cell Proportions (pl.DataFrame):")
    print(truth_df)

    print_section("5. MULTI-DATASET HARMONIZATION & INTEGRATION METRICS")
    # Simulate two independent cohorts with partially overlapping genes
    cohort_a = demo_adata[:60, :25].copy()
    cohort_b = demo_adata[60:, 5:].copy()
    print(f"Cohort A: {cohort_a.shape} | Cohort B: {cohort_b.shape}")

    # Mode 1: Intersection
    cfg_inter = HarmonizeConfig(mode=HarmonizeMode.INTERSECTION, min_shared_genes=10)
    inter_adata = align_and_concatenate([cohort_a, cohort_b], ["Cohort_A", "Cohort_B"], cfg_inter).unwrap()
    print(f"Harmonized (INTERSECTION): {inter_adata.shape} (Retains {inter_adata.n_vars} common genes)")

    # Mode 2: Union Zero-Filled
    cfg_union = HarmonizeConfig(mode=HarmonizeMode.UNION_ZERO_FILLED)
    union_adata = align_and_concatenate([cohort_a, cohort_b], ["Cohort_A", "Cohort_B"], cfg_union).unwrap()
    print(f"Harmonized (UNION ZERO-FILLED): {union_adata.shape} (Pads unmeasured genes with 0s)")

    # Integration Quality Metrics (iLISI, cLISI, Silhouette, kBET)
    metrics = evaluate_integration_metrics(inter_adata, batch_key="dataset_id", label_key="cell_type").unwrap()
    print("\nIntegration Quality & Mixing Metrics:")
    print(f"  * Mean iLISI (Cohort Mixing, ideal -> 2.0):    {metrics.mean_ilisi}")
    print(f"  * Mean cLISI (Cell-Type Purity, ideal -> 1.0): {metrics.mean_clisi}")
    print(f"  * Silhouette Ratio (Bio / Batch):              {metrics.silhouette_ratio}")
    print(f"  * kBET Acceptance Rate:                        {metrics.kbet_acceptance_rate}")

    print_section("6. PYTORCH DATASET & DATALOADER BRIDGE")
    torch_ds = TmeTorchDataset(
        demo_adata,
        label_keys=("response_binary", "cell_type"),
        on_the_fly_nb_dispersion=Some(0.15),
    )
    print(f"TmeTorchDataset initialized with {len(torch_ds)} cells and dynamic NB augmentation.")
    loader = create_tme_dataloader(torch_ds, batch_size=32, shuffle=True)

    batch_x, batch_y = next(iter(loader))
    print(f"Batched Expression Tensor: shape={batch_x.shape}, dtype={batch_x.dtype}")
    print(f"Batched Response Labels:   shape={batch_y['response_binary'].shape}")
    print(f"Batched Cell-Type Labels:  shape={batch_y['cell_type'].shape}")

    print_section("7. OUT-OF-CORE STORAGE & CRYPTOGRAPHIC VERIFICATION")
    with tempfile.TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        h5ad_path = tmp_dir / "demo.h5ad"
        demo_adata.write_h5ad(h5ad_path)

        # Cryptographic Hash & Verification
        hash_res = compute_file_hash(h5ad_path)
        sha256 = hash_res.unwrap()
        print(f"Computed SHA-256: {sha256[:16]}...{sha256[-16:]}")
        chk_res = verify_checksum(h5ad_path, ChecksumSpec(expected_hash=sha256))
        print(f"Checksum Verification: {'PASSED' if chk_res.unwrap() else 'FAILED'}")

        # Out-of-core Backed AnnData Loading
        backed_adata = load_backed(h5ad_path, backend=StorageBackend.BACKED_H5AD).unwrap()
        print(f"Loaded Backed AnnData without RAM loading: is_backed = {backed_adata.isbacked}")

        # Chunked Zarr Conversion
        zarr_path = tmp_dir / "demo.zarr"
        convert_to_zarr(demo_adata, zarr_path).unwrap()
        print(f"Converted to Chunked Zarr directory: {zarr_path.name}")

    print_section("TUTORIAL COMPLETE - ALL FEATURES VERIFIED SUCCESSFULLY!")


if __name__ == "__main__":
    main()
