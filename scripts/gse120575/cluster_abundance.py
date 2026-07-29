import argparse
import sys
from pathlib import Path
import anndata as ad  # type: ignore
import pertpy as pt  # type: ignore
import pandas as pd  # type: ignore
import scanpy as sc  # type: ignore
from sklearn.cluster import KMeans  # type: ignore
from returns.result import Result, Success, Failure  # type: ignore
from returns.pipeline import flow
from returns.pointfree import bind

def load_anndata(adata_path: Path) -> Result[ad.AnnData, str]:
    try:
        return Success(ad.read_h5ad(adata_path))
    except Exception as e:
        return Failure(f"Failed to load AnnData: {e}")

def compute_clusters(adata: ad.AnnData) -> Result[tuple[ad.AnnData, list[str]], str]:
    try:
        cluster_keys = []
        # Leiden resolutions
        for res in [0.5, 1.0, 1.5, 2.0]:
            key = f"leiden_{res}"
            sc.tl.leiden(adata, resolution=res, key_added=key)
            cluster_keys.append(key)
        
        # KMeans
        for k in [8, 10, 15]:
            kmeans = KMeans(n_clusters=k, random_state=42)
            adata.obs[f"kmeans_{k}"] = kmeans.fit_predict(adata.obsm["X_pca"])
            adata.obs[f"kmeans_{k}"] = adata.obs[f"kmeans_{k}"].astype(str).astype("category")
            cluster_keys.append(f"kmeans_{k}")
        
        return Success((adata, cluster_keys))
    except Exception as e:
        return Failure(f"Failed to compute clusters: {e}")

def export_proportions(adata: ad.AnnData, cluster_key: str, condition_name: str, out_dir: Path) -> None:
    df = adata.obs[["melanoma-sample", "response", cluster_key]].copy()
    counts = df.groupby(["melanoma-sample", "response", cluster_key], observed=False).size().unstack(fill_value=0)
    props = counts.div(counts.sum(axis=1), axis=0)
    props = props.reset_index()
    props.to_csv(out_dir / f"sccoda_proportions_{condition_name}_{cluster_key}.csv", index=False)

def run_sccoda_for_condition(adata: ad.AnnData, cluster_key: str, condition_name: str, out_dir: Path) -> Result[bool, str]:
    try:
        export_proportions(adata, cluster_key, condition_name, out_dir)
        
        sccoda = pt.tl.Sccoda()
        mdata = sccoda.load(
            adata,
            type="cell_level",
            generate_sample_level=True,
            cell_type_identifier=cluster_key,
            sample_identifier="melanoma-sample",
            covariate_obs=["response"]
        )
        
        # Choose the most abundant cluster as the reference
        counts = adata.obs[cluster_key].value_counts()
        ref_cluster = str(counts.idxmax())
        
        mdata = sccoda.prepare(
            mdata,
            formula="response",
            reference_cell_type=ref_cluster,
            automatic_reference_absence_threshold=0.5
        )
        sccoda.run_nuts(mdata, num_warmup=200, num_samples=1000, rng_key=42)
        
        # Lower threshold on inclusion probability by increasing FDR
        sccoda.set_fdr(mdata, est_fdr=0.2)
        
        varm_keys = mdata["coda"].varm.keys()
        effect_keys = [k for k in varm_keys if k.startswith("effect_df_")]
        if not effect_keys:
            raise ValueError(f"No effect_df found in varm.")
        
        effect_df = mdata["coda"].varm[effect_keys[0]]
        effect_df.to_csv(out_dir / f"sccoda_results_{condition_name}_{cluster_key}.csv")
        
        with open(out_dir / f"sccoda_summary_{condition_name}_{cluster_key}.txt", "w") as f:
            f.write(str(sccoda.summary(mdata)))
        
        return Success(True)
    except Exception as e:
        return Failure(f"Failed scCODA on {condition_name} ({cluster_key}): {e}")

def process_condition(adata: ad.AnnData, cluster_key: str, timepoint: str, condition_name: str, out_dir: Path) -> Result[bool, str]:
    try:
        if timepoint == "All":
            subset = adata.copy()
        else:
            subset = adata[adata.obs["treatment_status"] == timepoint].copy()
            
        if subset.n_obs == 0:
            return Failure(f"No cells found for {condition_name}")
            
        return run_sccoda_for_condition(subset, cluster_key, condition_name, out_dir)
    except Exception as e:
        return Failure(f"Failed to subset for {condition_name}: {e}")

def run_all_conditions(adata_and_keys: tuple[ad.AnnData, list[str]], out_dir: Path) -> Result[bool, str]:
    adata, cluster_keys = adata_and_keys
    adata.obs["response"] = adata.obs["characteristics: response"]
    
    for key in cluster_keys:
        for tp, cname in [("Pre", "Pre"), ("Post", "Post"), ("All", "Combined")]:
            print(f"Running scCODA for {cname} using {key}...")
            res = process_condition(adata, key, tp, cname, out_dir)
            if not isinstance(res, Success):
                print(f"Warning: {res.failure()}")
                
    return Success(True)

def run_pipeline(adata_path: Path, out_dir: Path) -> Result[bool, str]:
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        return flow(
            load_anndata(adata_path),
            bind(compute_clusters),
            bind(lambda ak: run_all_conditions(ak, out_dir))
        )
    except Exception as e:
        return Failure(str(e))

def main() -> None:
    parser = argparse.ArgumentParser(description="Run scCODA composition analysis")
    parser.add_argument("--adata", required=True, help="Input h5ad path")
    parser.add_argument("--out-dir", required=True, help="Output directory for CSVs")
    args = parser.parse_args()
    
    match run_pipeline(Path(args.adata), Path(args.out_dir)):
        case Success(_):
            print("Successfully completed scCODA analysis.")
            sys.exit(0)
        case Failure(err):
            print(f"Error: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
