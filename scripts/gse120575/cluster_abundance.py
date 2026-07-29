import argparse
import sys
from pathlib import Path
import anndata as ad  # type: ignore
import pertpy as pt  # type: ignore
import pandas as pd  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
from returns.result import Result, Success, Failure  # type: ignore
from returns.pipeline import flow
from returns.pointfree import bind

def load_anndata(adata_path: Path) -> Result[ad.AnnData, str]:
    try:
        return Success(ad.read_h5ad(adata_path))
    except Exception as e:
        return Failure(f"Failed to load AnnData: {e}")

def run_sccoda_for_condition(adata: ad.AnnData, condition_name: str, out_dir: Path) -> Result[bool, str]:
    try:
        # Rename column to avoid patsy parsing errors with colons and spaces
        adata.obs["response"] = adata.obs["characteristics: response"]
        
        sccoda = pt.tl.Sccoda()
        mdata = sccoda.load(
            adata,
            type="cell_level",
            generate_sample_level=True,
            cell_type_identifier="leiden",
            sample_identifier="melanoma-sample",
            covariate_obs=["response"]
        )
        mdata = sccoda.prepare(
            mdata,
            formula="response",
            reference_cell_type="automatic",
            automatic_reference_absence_threshold=0.5
        )
        sccoda.run_nuts(mdata, num_warmup=1000, num_samples=5000, rng_key=42)
        
        # Extract effect DataFrame
        # Patsy creates variables based on categories, e.g. 'response[T.Responder]'
        # pertpy stores this in varm as 'effect_df_...'
        varm_keys = mdata["coda"].varm.keys()
        effect_keys = [k for k in varm_keys if k.startswith("effect_df_")]
        if not effect_keys:
            raise ValueError(f"No effect_df found in varm. Keys: {list(varm_keys)}")
        
        effect_df = mdata["coda"].varm[effect_keys[0]]
        effect_df.to_csv(out_dir / f"sccoda_results_{condition_name}.csv")
        
        # Save a textual summary as well (optional, but helpful for human reading)
        with open(out_dir / f"sccoda_summary_{condition_name}.txt", "w") as f:
            f.write(str(sccoda.summary(mdata)))
        
        # Plot stacked barplot manually
        df = adata.obs[["melanoma-sample", "response", "leiden"]].copy()
        counts = df.groupby(["melanoma-sample", "response", "leiden"], observed=False).size().unstack(fill_value=0)
        props = counts.div(counts.sum(axis=1), axis=0)
        
        props = props.reset_index()
        props = props.sort_values(["response", "melanoma-sample"])
        
        ax = props.set_index("melanoma-sample").drop(columns=["response"]).plot(
            kind="bar", stacked=True, figsize=(12, 6), cmap="tab20"
        )
        labels = [f"{row['melanoma-sample']} ({row['response']})" for _, row in props.iterrows()]
        ax.set_xticklabels(labels, rotation=90)
        
        plt.title(f"Cluster Proportions ({condition_name})")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Cluster")
        plt.ylabel("Proportion")
        plt.tight_layout()
        plt.savefig(out_dir / f"sccoda_abundance_{condition_name}.png", dpi=300)
        plt.close()
        
        return Success(True)
    except Exception as e:
        return Failure(f"Failed scCODA on {condition_name}: {e}")

def process_condition(adata: ad.AnnData, timepoint: str, condition_name: str, out_dir: Path) -> Result[bool, str]:
    try:
        if timepoint == "All":
            subset = adata.copy()
        else:
            subset = adata[adata.obs["treatment_status"] == timepoint].copy()
            
        if subset.n_obs == 0:
            return Failure(f"No cells found for {condition_name}")
            
        return run_sccoda_for_condition(subset, condition_name, out_dir)
    except Exception as e:
        return Failure(f"Failed to subset for {condition_name}: {e}")

def run_all_conditions(adata: ad.AnnData, out_dir: Path) -> Result[bool, str]:
    res_pre = process_condition(adata, "Pre", "Pre", out_dir)
    if not isinstance(res_pre, Success): return res_pre
    
    res_post = process_condition(adata, "Post", "Post", out_dir)
    if not isinstance(res_post, Success): return res_post
    
    res_all = process_condition(adata, "All", "Combined", out_dir)
    if not isinstance(res_all, Success): return res_all
    
    return Success(True)

def run_pipeline(adata_path: Path, out_dir: Path) -> Result[bool, str]:
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        return flow(
            load_anndata(adata_path),
            bind(lambda adata: run_all_conditions(adata, out_dir))
        )
    except Exception as e:
        return Failure(str(e))

def main() -> None:
    parser = argparse.ArgumentParser(description="Run scCODA composition analysis")
    parser.add_argument("--adata", required=True, help="Input h5ad path")
    parser.add_argument("--out-dir", required=True, help="Output directory for plots and CSVs")
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
