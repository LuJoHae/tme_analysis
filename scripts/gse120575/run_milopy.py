import argparse
import sys
from pathlib import Path
import anndata as ad  # type: ignore
import milopy  # type: ignore
from returns.result import Result, Success, Failure  # type: ignore
from returns.pipeline import flow  # type: ignore
from returns.pointfree import bind  # type: ignore

def load_anndata(adata_path: Path) -> Result[ad.AnnData, str]:
    try:
        return Success(ad.read_h5ad(adata_path))
    except Exception as e:
        return Failure(f"Failed to load AnnData: {e}")

def run_milopy_for_condition(adata: ad.AnnData, condition_name: str, out_dir: Path) -> Result[bool, str]:
    try:
        subset = adata.copy()
        
        milopy.core.make_nhoods(subset, prop=0.1, k=50, d=40)
        milopy.core.count_cells(subset, sample_col="melanoma-sample")
        design_df = subset.obs[["melanoma-sample", "response"]].drop_duplicates().set_index("melanoma-sample")
        milopy.core.test_nhoods(subset, design="~response", design_df=design_df)
        
        if "nhood_test_results" not in subset.uns:
            return Failure("nhood_test_results not found in subset.uns")
            
        res_df = subset.uns["nhood_test_results"]
        
        # Save results CSV
        res_df.to_csv(out_dir / f"milopy_results_{condition_name}.csv")
        
        # Save the neighborhood assignment sparse matrix for downstream plotting
        if "nhoods" in subset.obsm:
            import scipy.sparse as sp
            sp.save_npz(out_dir / f"milopy_nhoods_{condition_name}.npz", subset.obsm["nhoods"])
            
        return Success(True)
    except Exception as e:
        return Failure(f"Failed milopy on {condition_name}: {e}")

def process_condition(adata: ad.AnnData, timepoint: str, condition_name: str, out_dir: Path) -> Result[bool, str]:
    try:
        if timepoint == "All":
            subset = adata.copy()
        else:
            subset = adata[adata.obs["treatment_status"] == timepoint].copy()
            
        if subset.n_obs == 0:
            return Failure(f"No cells found for {condition_name}")
            
        return run_milopy_for_condition(subset, condition_name, out_dir)
    except Exception as e:
        return Failure(f"Failed to subset for {condition_name}: {e}")

def run_all_conditions(adata: ad.AnnData, out_dir: Path) -> Result[bool, str]:
    try:
        adata.obs["response"] = adata.obs["characteristics: response"]
    except Exception as e:
        return Failure(f"Failed to find characteristics: response in adata.obs: {e}")

    for tp, cname in [("Pre", "Pre"), ("Post", "Post"), ("All", "Combined")]:
        print(f"Running milopy for {cname}...")
        res = process_condition(adata, tp, cname, out_dir)
        match res:
            case Failure(err):
                print(f"Warning: {err}")
            case Success(_):
                pass
                
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
    parser = argparse.ArgumentParser(description="Run milopy DA analysis")
    parser.add_argument("--adata", required=True, help="Input h5ad path")
    parser.add_argument("--out-dir", required=True, help="Output directory for CSVs")
    args = parser.parse_args()
    
    match run_pipeline(Path(args.adata), Path(args.out_dir)):
        case Success(_):
            print("Successfully completed milopy analysis.")
            sys.exit(0)
        case Failure(err):
            print(f"Error: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
