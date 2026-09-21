"""Functional validator for milopy differential abundance readiness across preprocessed AnnData datasets."""

from pathlib import Path
from pydantic import BaseModel, ConfigDict
import anndata as ad
import scanpy as sc
from scipy.sparse import issparse
from returns.result import Result, Success, Failure
import milopy  # type: ignore


class MilopyCheckResult(BaseModel):
    """Immutable metadata result of milopy compatibility verification on a single .h5ad file."""
    model_config = ConfigDict(frozen=True)

    accession: str
    n_obs: int
    n_vars: int
    has_sparse_x: bool
    has_pca: bool
    has_knn: bool
    sample_col: str
    n_samples: int
    design_col: str
    n_conditions: int
    cell_type_col: str
    milopy_dry_run_pass: bool
    status_summary: str
    error_detail: str = ""


def find_column(obs_cols: list[str], candidates: list[str]) -> str:
    """Finds first matching column name from candidates in obs_cols, or returns 'Missing'."""
    for c in candidates:
        if c in obs_cols:
            return c
    return "Missing"


def verify_milopy_compatibility(adata_path: Path) -> Result[MilopyCheckResult, str]:
    """Evaluates whether a preprocessed .h5ad AnnData object meets all requirements for milopy DA analysis."""
    try:
        accession = adata_path.name.replace("_processed.h5ad", "").replace(".h5ad", "")
        adata = ad.read_h5ad(adata_path)

        n_obs = adata.n_obs
        n_vars = adata.n_vars
        has_sparse_x = issparse(adata.X)
        has_pca = "X_pca" in adata.obsm
        has_knn = "connectivities" in adata.obsp

        obs_cols = list(adata.obs.columns)
        sample_col = find_column(obs_cols, ["patient", "patient_id", "Patient_ID", "sample", "Sample", "original.barcode"])
        n_samples = adata.obs[sample_col].nunique() if sample_col != "Missing" else 0

        design_col = find_column(obs_cols, ["response", "therapy", "treatment_status", "cancer_code", "cohort"])
        n_conditions = adata.obs[design_col].nunique() if design_col != "Missing" else 0

        cell_type_col = find_column(obs_cols, ["cell_type", "cell_type_main", "celltype", "CellType", "major_cell_type"])

        # Check basic prerequisite eligibility
        if n_obs < 50 or n_vars < 100:
            return Success(MilopyCheckResult(
                accession=accession, n_obs=n_obs, n_vars=n_vars,
                has_sparse_x=has_sparse_x, has_pca=has_pca, has_knn=has_knn,
                sample_col=sample_col, n_samples=n_samples,
                design_col=design_col, n_conditions=n_conditions,
                cell_type_col=cell_type_col, milopy_dry_run_pass=False,
                status_summary="FAIL_INSUFFICIENT_CELLS",
                error_detail=f"Matrix dimensions too small ({n_obs} cells, {n_vars} genes)",
            ))

        if sample_col == "Missing" or n_samples < 2:
            return Success(MilopyCheckResult(
                accession=accession, n_obs=n_obs, n_vars=n_vars,
                has_sparse_x=has_sparse_x, has_pca=has_pca, has_knn=has_knn,
                sample_col=sample_col, n_samples=n_samples,
                design_col=design_col, n_conditions=n_conditions,
                cell_type_col=cell_type_col, milopy_dry_run_pass=False,
                status_summary="FAIL_SINGLE_SAMPLE",
                error_detail=f"Sample column '{sample_col}' has {n_samples} samples (milopy requires >= 2 samples)",
            ))

        if design_col == "Missing" or n_conditions < 2:
            return Success(MilopyCheckResult(
                accession=accession, n_obs=n_obs, n_vars=n_vars,
                has_sparse_x=has_sparse_x, has_pca=has_pca, has_knn=has_knn,
                sample_col=sample_col, n_samples=n_samples,
                design_col=design_col, n_conditions=n_conditions,
                cell_type_col=cell_type_col, milopy_dry_run_pass=False,
                status_summary="FAIL_NO_DESIGN_CONTRAST",
                error_detail=f"Design column '{design_col}' has {n_conditions} conditions (milopy requires >= 2 levels)",
            ))

        # Perform lightweight dry-run milopy test on a subset
        dry_run_pass = False
        err_msg = ""
        try:
            # Filter out 'Unknown' response cells — they cannot be assigned to a design contrast
            if design_col != "Missing":
                adata = adata[adata.obs[design_col] != "Unknown"].copy()
                n_conditions = adata.obs[design_col].nunique()
                n_samples = adata.obs[sample_col].nunique() if sample_col != "Missing" else 0
                if n_conditions < 2:
                    return Success(MilopyCheckResult(
                        accession=accession, n_obs=adata.n_obs, n_vars=adata.n_vars,
                        has_sparse_x=issparse(adata.X), has_pca="X_pca" in adata.obsm,
                        has_knn="connectivities" in adata.obsp,
                        sample_col=sample_col, n_samples=n_samples,
                        design_col=design_col, n_conditions=n_conditions,
                        cell_type_col=cell_type_col, milopy_dry_run_pass=False,
                        status_summary="FAIL_NO_DESIGN_CONTRAST",
                        error_detail=f"After filtering Unknown cells, design column has {n_conditions} conditions",
                    ))
                if adata.n_obs < 50:
                    return Success(MilopyCheckResult(
                        accession=accession, n_obs=adata.n_obs, n_vars=adata.n_vars,
                        has_sparse_x=issparse(adata.X), has_pca=False, has_knn=False,
                        sample_col=sample_col, n_samples=n_samples,
                        design_col=design_col, n_conditions=n_conditions,
                        cell_type_col=cell_type_col, milopy_dry_run_pass=False,
                        status_summary="FAIL_INSUFFICIENT_CELLS",
                        error_detail="Too few cells after filtering Unknown response",
                    ))

            # Check minimum patients per condition (milopy GLM needs >= 2 replicates per group)
            if sample_col != "Missing" and design_col != "Missing":
                per_condition = adata.obs.groupby(design_col, observed=True)[sample_col].nunique()
                if (per_condition < 2).any():
                    return Success(MilopyCheckResult(
                        accession=accession, n_obs=adata.n_obs, n_vars=adata.n_vars,
                        has_sparse_x=issparse(adata.X), has_pca="X_pca" in adata.obsm,
                        has_knn="connectivities" in adata.obsp,
                        sample_col=sample_col, n_samples=n_samples,
                        design_col=design_col, n_conditions=n_conditions,
                        cell_type_col=cell_type_col, milopy_dry_run_pass=False,
                        status_summary="FAIL_INSUFFICIENT_REPLICATES",
                        error_detail=f"Each condition needs >= 2 patients; got: {per_condition.to_dict()}",
                    ))

            sub_size = min(2000, adata.n_obs)
            sub_idx = adata.obs.sample(sub_size, random_state=42).index if adata.n_obs > 2000 else adata.obs.index
            sub = adata[sub_idx, :].copy()

            if "X_pca" not in sub.obsm:
                sc.pp.pca(sub, n_comps=min(30, sub.n_vars - 1))
            if "connectivities" not in sub.obsp:
                sc.pp.neighbors(sub, n_neighbors=min(15, sub_size - 1))

            k_val = min(30, max(5, sub_size // 10))
            milopy.core.make_nhoods(sub, prop=0.1, k=k_val)
            milopy.core.count_cells(sub, sample_col=sample_col)
            # Only include samples present in the subsample — absent samples have 0 library size
            # which causes edgepython TMM normalization to fail
            present_samples = set(sub.obs[sample_col].unique())
            design_df = (
                sub.obs[[sample_col, design_col]]
                .drop_duplicates()
                .set_index(sample_col)
                .loc[lambda df: df.index.isin(present_samples)]
            )
            # Re-check minimum patients per condition after subsample filtering
            per_cond_sub = design_df.groupby(design_col, observed=True).size()
            if (per_cond_sub < 2).any():
                err_msg = f"After subsampling, condition has <2 patients: {per_cond_sub.to_dict()}"
                raise ValueError(err_msg)

            milopy.core.test_nhoods(sub, design=f"~{design_col}", design_df=design_df)

            if "nhood_test_results" in sub.uns:
                dry_run_pass = True
            else:
                err_msg = "nhood_test_results missing in uns after test_nhoods"
        except Exception as ex:
            err_msg = f"Milopy execution error: {str(ex)}"


        status = "PASS_MILOPY_READY" if dry_run_pass else "FAIL_DRY_RUN_ERROR"

        return Success(MilopyCheckResult(
            accession=accession,
            n_obs=n_obs,
            n_vars=n_vars,
            has_sparse_x=has_sparse_x,
            has_pca=has_pca,
            has_knn=has_knn,
            sample_col=sample_col,
            n_samples=n_samples,
            design_col=design_col,
            n_conditions=n_conditions,
            cell_type_col=cell_type_col,
            milopy_dry_run_pass=dry_run_pass,
            status_summary=status,
            error_detail=err_msg,
        ))
    except Exception as e:
        return Failure(f"Failed to check AnnData file {adata_path.name}: {str(e)}")
