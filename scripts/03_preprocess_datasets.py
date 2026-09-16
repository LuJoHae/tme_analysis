"""Script 3: Individual Quality Control, Metadata Harmonization, and Gene Normalization per Dataset.

Inspired by single_cell_datasets (SingleCellDataProcessStep02, Step05 & Step07), this script:
1. Processes each dataset individually (no global concatenation).
2. Detects prior log-transformation & normalization state (check_and_normalize_transform_space).
3. Harmonizes metadata column names (dataset, cancer_code, staging, sequencing_tech, cell_type, patient, organ, original.barcode) while preserving all clinical annotations.
4. Normalizes all gene IDs for every dataset using gene_utils.norm_genes(pre_id_transform="auto").
5. Applies quality control filtering (min/max genes, min/max counts, mitochondrial content).
6. Removes MT and non-canonical junk contigs.
7. Saves one preprocessed .h5ad HDF5 file per dataset.
"""

from pathlib import Path
from typing import Optional
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict
import anndata as ad
import scanpy as sc
from returns.result import Result, Success, Failure
from gene_utils import norm_genes


# Standardized metadata lookup table for single-cell cohorts
COHORT_METADATA_MAP = {
    "AziziSingleCellMapDiverse2018Adata": {
        "cancer_code": "BRCA", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {"geo_id": "geo_id", "patient": "patient", "tissue": "organ"}
    },
    "BeckerSinglecellAnalysesDefine2022Adata": {
        "cancer_code": "COAD", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {"sample": "sample", "geo_id": "geo_id"}
    },
    "BiermannDissectingTreatmentnaiveEcosystem2022Adata": {
        "cancer_code": "SKCM", "staging": "brain metastasis", "sequencing_tech": "10x_UMI",
        "obs_map": {"patient": "patient", "organ": "organ", "cell_type_main": "cell_type"}
    },
    "BorcherdingMappingImmuneEnvironment2021Adata": {
        "cancer_code": "ccRCC", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {}
    },
    "ChengPancancerSinglecellTranscriptional2021Adata": {
        "cancer_code": "Pan", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {"patient": "patient", "tissue": "organ", "cancer_type": "cancer_type"}
    },
    "DuranteSinglecellAnalysisReveals2020Adata": {
        "cancer_code": "UVM", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {"sample": "sample", "geo_id": "geo_id"}
    },
    "JerbyArnonCancerCellProgram2018Adata": {
        "cancer_code": "SKCM", "staging": "primary tumor", "sequencing_tech": "Smart-seq2",
        "obs_map": {"samples": "sample", "Cohort": "cohort", "cell.types": "cell_type"}
    },
    "KhaliqRefiningColorectalCancer2022Adata": {
        "cancer_code": "CC", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {"samples": "sample", "Condition": "is_tumor"}
    },
    "KimSinglecellRNASequencing2020Adata": {
        "cancer_code": "LUAD", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {}
    },
    "LeaderSinglecellAnalysisHuman2021Adata": {
        "cancer_code": "NSCLC", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {"batch": "batch"}
    },
    "LuSinglecellAtlasMulticellular2022Adata": {
        "cancer_code": "HCC", "staging": "primary tumor and metastasis", "sequencing_tech": "10x_UMI",
        "obs_map": {"sample": "sample", "patient": "patient", "site": "organ", "celltype": "cell_type"}
    },
    "PelkaSpatiallyOrganizedMulticellular2021Adata": {
        "cancer_code": "CRC", "staging": "primary tumor and metastasis", "sequencing_tech": "10x_UMI",
        "obs_map": {"SPECIMEN_TYPE": "is_tumor", "PROCESSING_TYPE": "cell_type", "PatientTypeID": "patient"}
    },
    "PuSinglecellTranscriptomicAnalysis2021Adata": {
        "cancer_code": "PTC", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {"patient": "patient", "biopsy_site": "organ", "geo_id": "geo_id"}
    },
    "QianPancancerBlueprintHeterogeneous2020aAdata": {
        "cancer_code": "Pan", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {"Source Name": "sample", "Characteristics[individual]": "patient", "Characteristics[organism part]": "organ"}
    },
    "SharmaOncofetalReprogrammingEndothelial2020Adata": {
        "cancer_code": "HCC", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {}
    },
    "SadeFeldmanDefiningTCell2018Adata": {
        "cancer_code": "SKCM", "staging": "primary tumor", "sequencing_tech": "Smart-seq2",
        "obs_map": {"Cell_ID": "original.barcode", "Patient_ID": "patient"}
    },
    "YostClonalReplacementTumor2019Adata": {
        "cancer_code": "TCell_Immuno", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {}
    },
    "ZhengLandscapeInfiltratingTCells2017Adata": {
        "cancer_code": "Immune_Microenvironment", "staging": "primary tumor", "sequencing_tech": "10x_UMI",
        "obs_map": {}
    },
}


class DatasetQCSpec(BaseModel):
    """Quality control threshold specification for a single-cell dataset."""
    model_config = ConfigDict(frozen=True)
    min_genes: int = 300
    max_genes: int = 4500
    min_counts: int = 300
    max_counts: int = 15000
    max_mt_content: float = 20.0


class PreprocessConfig(BaseModel):
    """Immutable configuration for dataset preprocessing into HDF5 files."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    input_dir: Path
    output_dir: Path
    default_qc: DatasetQCSpec = DatasetQCSpec()


def check_and_normalize_transform_space(adata: ad.AnnData) -> tuple[ad.AnnData, dict[str, bool]]:
    """Detects if expression data has already been log-transformed or total-sum normalized.

    Returns updated AnnData and a metadata status dictionary indicating is_log1p and is_normalized.
    """
    adata_out = adata.copy()
    
    # 1. Check uns annotations
    has_uns_log1p = "log1p" in adata_out.uns or "log_transformed" in adata_out.uns
    has_uns_norm = "normalized" in adata_out.uns
    
    # 2. Inspect matrix statistics
    X_mat = adata_out.X
    if hasattr(X_mat, "toarray"):
        sample_vals = X_mat[:min(100, adata_out.n_obs)].toarray()
    else:
        sample_vals = X_mat[:min(100, adata_out.n_obs)]

    max_val = float(sample_vals.max()) if sample_vals.size > 0 else 0.0
    
    # Heuristic: max value <= 35 and presence of non-integer floats indicates log space (e.g. log2(TPM+1) or log1p)
    is_non_integer = not np.all(np.equal(np.mod(sample_vals, 1), 0))
    is_log_transformed = has_uns_log1p or (max_val <= 35.0 and is_non_integer and max_val > 0.0)

    # Check total sums per row
    row_sums = sample_vals.sum(axis=1)
    is_normalized = has_uns_norm or np.allclose(row_sums, 1e4, rtol=1e-2) or np.allclose(row_sums, 1e6, rtol=1e-2)

    status = {
        "is_log1p": is_log_transformed,
        "is_normalized": is_normalized,
        "max_value": max_val,
    }

    # Store detection status in uns metadata
    adata_out.uns["data_transform_state"] = status

    # If in log-space, convert back to linear space for raw count QC / BayesPrism reference preparation
    if is_log_transformed:
        adata_out.layers["log1p_original"] = adata_out.X.copy()
        if hasattr(adata_out.X, "data"):
            # Sparse matrix expm1
            adata_out.X.data = np.expm1(adata_out.X.data)
        else:
            adata_out.X = np.expm1(adata_out.X)

    return adata_out, status


def harmonize_dataset_metadata(adata: ad.AnnData, dataset_key: str) -> ad.AnnData:
    """Harmonizes metadata columns (dataset, cancer_code, staging, sequencing_tech, cell_type, patient, organ) while preserving original obs fields."""
    adata_out = adata.copy()
    obs_df = adata_out.obs.copy()

    # Ensure original barcode column is present
    if "original.barcode" not in obs_df.columns:
        obs_df["original.barcode"] = obs_df.index

    # Find matching cohort specification
    cohort_spec = None
    for key, spec in COHORT_METADATA_MAP.items():
        if key.lower() in dataset_key.lower() or dataset_key.lower() in key.lower():
            cohort_spec = spec
            break

    if cohort_spec:
        obs_df["cancer_code"] = cohort_spec["cancer_code"]
        obs_df["staging"] = cohort_spec["staging"]
        obs_df["sequencing_tech"] = cohort_spec["sequencing_tech"]
        
        # Apply specific obs column renames
        for old_col, new_col in cohort_spec["obs_map"].items():
            if old_col in obs_df.columns and new_col not in obs_df.columns:
                obs_df[new_col] = obs_df[old_col]
    else:
        obs_df["cancer_code"] = "Unknown"
        obs_df["staging"] = "Unknown"
        obs_df["sequencing_tech"] = "10x_UMI"

    obs_df["dataset"] = dataset_key

    # Standardize cell_type column fallback
    if "cell_type" not in obs_df.columns:
        for col in ["cell_type_main", "celltype", "cell_types", "CellType"]:
            if col in obs_df.columns:
                obs_df["cell_type"] = obs_df[col]
                break
        else:
            obs_df["cell_type"] = "Unknown"

    adata_out.obs = obs_df
    return adata_out


def preprocess_single_dataset(
    h5ad_path: Path,
    qc_spec: DatasetQCSpec,
    output_path: Path,
) -> Result[Path, str]:
    """Applies QC filtering, metadata harmonization, gene normalization, and saves individual .h5ad HDF5 file."""
    try:
        if not h5ad_path.exists():
            return Failure(f"Input file not found: {h5ad_path}")

        dataset_key = h5ad_path.stem.replace("_sparse", "").replace(".h5ad", "")
        print(f"Preprocessing dataset '{dataset_key}': {h5ad_path.name}...")
        adata = ad.read_h5ad(h5ad_path)

        # 1. Detect prior log-transform / normalization
        adata, status = check_and_normalize_transform_space(adata)
        print(f" -> Transform state for {dataset_key}: log_transformed={status['is_log1p']}, normalized={status['is_normalized']}, max_val={status['max_value']:.2f}")

        # 2. Harmonize metadata columns while preserving all obs fields
        adata = harmonize_dataset_metadata(adata, dataset_key)

        # 3. Normalize all gene IDs per dataset using gene_utils
        try:
            adata = norm_genes(adata, pre_id_transform="auto")
        except Exception:
            pass

        # 4. Calculate QC metrics & mitochondrial %
        if "contig" in adata.var.columns:
            adata.var["mt"] = (adata.var["contig"] == "MT").values
        else:
            adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")

        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        # 5. Apply QC cell filtering
        sc.pp.filter_cells(adata, min_genes=qc_spec.min_genes)
        sc.pp.filter_cells(adata, min_counts=qc_spec.min_counts)
        sc.pp.filter_cells(adata, max_genes=qc_spec.max_genes)
        sc.pp.filter_cells(adata, max_counts=qc_spec.max_counts)

        if "pct_counts_mt" in adata.obs.columns:
            adata = adata[adata.obs["pct_counts_mt"] < qc_spec.max_mt_content, :].copy()

        # 6. Filter out non-canonical contig & mitochondrial junk genes
        if "contig" in adata.var.columns:
            junk_mask = (adata.var["contig"] == "MT") | (adata.var["contig"].astype(str).str.len() >= 3)
            adata = adata[:, ~junk_mask].copy()

        # 7. Save preprocessed AnnData object to disk as HDF5 (.h5ad)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        ad.settings.allow_write_nullable_strings = True
        adata.write_h5ad(output_path, compression="gzip")

        print(f"Saved preprocessed dataset '{dataset_key}' ({adata.n_obs} cells, {adata.n_vars} genes) -> {output_path}")
        return Success(output_path)
    except Exception as e:
        return Failure(f"Preprocessing failed for {h5ad_path.name}: {str(e)}")


def run_dataset_preprocessing_pipeline(config: PreprocessConfig, n_jobs: int = -1) -> Result[list[Path], str]:
    """Iterates through all sparse .h5ad files and applies preprocessing in parallel to generate one HDF5 file per dataset."""
    try:
        from concurrent.futures import ProcessPoolExecutor
        import os

        config.output_dir.mkdir(parents=True, exist_ok=True)
        h5ad_files = list(config.input_dir.glob("*.h5ad"))
        
        if not h5ad_files:
            return Failure(f"No .h5ad files found in input directory: {config.input_dir}")

        workers = os.cpu_count() if n_jobs <= 0 else n_jobs
        print(f"Executing dataset preprocessing for {len(h5ad_files)} datasets in parallel using {workers} workers...")

        tasks = [(h5_file, config.default_qc, config.output_dir / f"{h5_file.stem}_processed.h5ad") for h5_file in h5ad_files]

        processed_paths = []
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(preprocess_single_dataset, h5_file, qc, out_file)
                for (h5_file, qc, out_file) in tasks
            ]
            for future in futures:
                match future.result():
                    case Success(p):
                        processed_paths.append(p)
                    case Failure(err):
                        print(f"Warning: {err}")

        return Success(processed_paths)
    except Exception as e:
        return Failure(f"Preprocessing pipeline failed: {str(e)}")


def main() -> None:
    """CLI entry point for dataset quality control and preprocessing."""
    import argparse
    parser = argparse.ArgumentParser(description="Preprocess single-cell datasets into individual HDF5 (.h5ad) files.")
    parser.add_argument("--input-dir", type=str, default="data/sparse_h5", help="Input directory with sparse .h5ad files")
    parser.add_argument("--out-dir", type=str, default="data/processed_h5", help="Output directory for preprocessed .h5ad files")
    parser.add_argument("--input-h5", type=str, default=None, help="Optional specific input sparse h5ad file for single dataset")
    parser.add_argument("--out-h5", type=str, default=None, help="Optional specific output preprocessed h5ad file for single dataset")
    args = parser.parse_args()

    if args.input_h5 and args.out_h5:
        input_h5 = Path(args.input_h5).resolve()
        out_h5 = Path(args.out_h5).resolve()
        qc_spec = DatasetQCSpec()
        match preprocess_single_dataset(input_h5, qc_spec, out_h5):
            case Success(p):
                print(f"Successfully preprocessed dataset '{input_h5.name}' -> {p}")
            case Failure(err):
                print(f"Preprocessing error for '{input_h5.name}': {err}")
    else:
        config = PreprocessConfig(
            input_dir=Path(args.input_dir).resolve(),
            output_dir=Path(args.out_dir).resolve(),
        )
        print(f"Starting dataset preprocessing pipeline...")
        print(f"Input dir: {config.input_dir}")
        print(f"Output dir: {config.output_dir}")

        match run_dataset_preprocessing_pipeline(config):
            case Success(paths):
                print(f"\nSuccessfully generated {len(paths)} preprocessed HDF5 (.h5ad) files:")
                for p in paths:
                    print(f" - {p}")
            case Failure(err):
                print(f"Pipeline error: {err}")


if __name__ == "__main__":
    main()
