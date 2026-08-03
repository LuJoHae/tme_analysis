"""Processing script to parse GSE123139 scRNA-seq & TCR data into AnnData HDF5 (.h5ad) files."""

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Final, assert_never

import anndata as ad  # type: ignore
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import polars as pl  # type: ignore
import scipy.sparse as sp  # type: ignore
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success


class ProcessConfig(BaseModel):
    """Immutable data processing configuration for GSE123139."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    input_dir: Path
    out_h5ad: Path
    tcr_file: Maybe[Path] = Nothing


def parse_series_matrix_metadata(matrix_path: Path) -> Result[pl.DataFrame, str]:
    """Parse GEO series matrix file to extract sample annotations."""
    try:
        sample_ids: list[str] = []
        characteristics_dict: dict[str, list[str]] = {}

        open_fn = gzip.open if matrix_path.suffix == ".gz" else open
        with open_fn(matrix_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("!Sample_geo_accession"):
                    sample_ids = [s.strip('"') for s in line.strip().split("\t")[1:]]
                elif line.startswith("!Sample_characteristics_ch1"):
                    parts = [p.strip('"') for p in line.strip().split("\t")[1:]]
                    if parts and ":" in parts[0]:
                        key = parts[0].split(":")[0].strip().lower().replace(" ", "_")
                        vals = [p.split(":", 1)[1].strip() if ":" in p else p for p in parts]
                        characteristics_dict[key] = vals

        if not sample_ids:
            return Failure("Could not locate !Sample_geo_accession in series matrix file.")

        data: dict[str, list[str]] = {"sample_id": sample_ids}
        for k, v in characteristics_dict.items():
            if len(v) == len(sample_ids):
                data[k] = v

        meta_df = pl.DataFrame(data)
        return Success(meta_df)
    except Exception as err:
        return Failure(f"Failed to parse series matrix file {matrix_path}: {err}")


def load_tcr_metadata(tcr_path: Path) -> Result[pl.DataFrame, str]:
    """Load single-cell TCR-seq metadata table using Polars."""
    try:
        tcr_df = pl.read_csv(
            tcr_path,
            separator="\t",
            truncate_ragged_lines=True,
            encoding="utf8-lossy",
        )
        return Success(tcr_df)
    except Exception as err:
        return Failure(f"Failed to load TCR file {tcr_path}: {err}")


def process_raw_count_files(raw_dir: Path) -> Result[tuple[pl.DataFrame, list[str], list[str]], str]:
    """Load individual sample expression text files and combine into a single expression matrix."""
    try:
        count_files = sorted(list(raw_dir.glob("*.txt.gz")) + list(raw_dir.glob("*.txt")))
        if not count_files:
            return Failure(f"No .txt or .txt.gz raw count files found in {raw_dir}")

        all_cell_names: list[str] = []
        combined_df: pl.DataFrame | None = None

        for idx, file_path in enumerate(count_files):
            sample_match = re.search(r"(GSM\d+)", file_path.name)
            sample_prefix = sample_match.group(1) if sample_match else f"Sample{idx}"

            df = pl.read_csv(
                file_path,
                separator="\t",
                has_header=True,
                truncate_ragged_lines=True,
            )

            gene_col = df.columns[0]
            expr_cols = df.columns[1:]

            renamed_cols = {col: f"{sample_prefix}_{col}" for col in expr_cols}
            df_renamed = df.select([pl.col(gene_col).alias("gene")] + [pl.col(c).alias(renamed_cols[c]) for c in expr_cols])

            all_cell_names.extend(list(renamed_cols.values()))

            if combined_df is None:
                combined_df = df_renamed
            else:
                combined_df = combined_df.join(df_renamed, on="gene", how="full", coalesce=True)

        if combined_df is None:
            return Failure("No count data loaded.")

        combined_df = combined_df.fill_null(0.0)
        genes = combined_df["gene"].to_list()
        expr_df = combined_df.drop("gene")

        return Success((expr_df, genes, all_cell_names))
    except Exception as err:
        return Failure(f"Failed processing raw count files: {err}")


def build_anndata_h5ad(
    expr_df: pl.DataFrame,
    genes: list[str],
    cell_names: list[str],
    series_meta: pl.DataFrame,
    tcr_meta: Maybe[pl.DataFrame],
    out_h5ad: Path,
) -> Result[Path, str]:
    """Construct AnnData object and serialize to HDF5 (.h5ad) format."""
    try:
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)

        # Expression values (cells x genes)
        X_matrix = expr_df.to_numpy().T.astype(np.float32)
        X_sparse = sp.csr_matrix(X_matrix)

        str_genes = [str(g) for g in genes]
        str_cells = [str(c) for c in cell_names]

        var_df = pd.DataFrame(index=str_genes)
        var_df.index.name = "gene_symbol"
        var_df.index = var_df.index.astype(str)

        obs_df = pd.DataFrame(index=str_cells)
        obs_df.index.name = "cell_id"
        obs_df.index = obs_df.index.astype(str)

        # Assign sample ID and well ID suffix matching
        obs_df["sample_id"] = [cell_id.split("_")[0] for cell_id in str_cells]
        obs_df["well_id"] = [cell_id.split("_")[-1] for cell_id in str_cells]

        # Merge series matrix annotations
        meta_pd = series_meta.to_pandas()
        meta_pd["sample_id"] = meta_pd["sample_id"].astype(str)
        meta_pd = meta_pd.set_index("sample_id")
        obs_df = obs_df.join(meta_pd, on="sample_id", how="left")

        # Optionally join TCR annotations if present matching well_id / Well_ID
        match tcr_meta:
            case Some(tcr_df):
                tcr_pd = tcr_df.to_pandas()
                key_col = "Well_ID" if "Well_ID" in tcr_pd.columns else ("well_id" if "well_id" in tcr_pd.columns else ("cell_id" if "cell_id" in tcr_pd.columns else None))
                if key_col:
                    tcr_pd[key_col] = tcr_pd[key_col].astype(str)
                    tcr_pd = tcr_pd.drop_duplicates(subset=[key_col]).set_index(key_col)
                    obs_df = obs_df.join(tcr_pd, on="well_id", how="left")
            case _:
                pass

        # Convert all object and categorical columns in obs to string to ensure clean HDF5 serialization
        for col in obs_df.columns:
            if obs_df[col].dtype == object or str(obs_df[col].dtype) == "category":
                obs_df[col] = obs_df[col].fillna("").astype(str)

        adata = ad.AnnData(X=X_sparse, obs=obs_df, var=var_df)

        # Write to HDF5 (.h5ad) format
        print(f"Writing AnnData object ({adata.n_obs} cells x {adata.n_vars} genes) to {out_h5ad}...")
        adata.write_h5ad(out_h5ad, compression="gzip")

        return Success(out_h5ad)
    except Exception as err:
        return Failure(f"Failed creating AnnData HDF5 file: {err}")


def run_pipeline(config: ProcessConfig) -> Result[Path, str]:
    """Execute complete functional pipeline to construct GSE123139 HDF5 file."""
    raw_dir = config.input_dir / "raw_counts"
    matrix_file = config.input_dir / "GSE123139_series_matrix.txt.gz"
    tcr_file = config.input_dir / "GSE123139_T_cells_tcrb_v2.txt.gz"

    # 1. Parse Series Matrix Metadata
    match parse_series_matrix_metadata(matrix_file):
        case Failure(err):
            return Failure(err)
        case Success(series_meta):
            # 2. Parse optional TCR Metadata
            tcr_meta: Maybe[pl.DataFrame] = Nothing
            if tcr_file.exists():
                match load_tcr_metadata(tcr_file):
                    case Success(df):
                        tcr_meta = Some(df)
                    case Failure(_):
                        tcr_meta = Nothing
                    case _ as unreachable:
                        assert_never(unreachable)

            # 3. Process raw count matrices
            target_raw_dir = raw_dir if raw_dir.exists() else config.input_dir
            match process_raw_count_files(target_raw_dir):
                case Failure(err):
                    return Failure(err)
                case Success((expr_df, genes, cell_names)):
                    # 4. Build and output AnnData HDF5 file
                    return build_anndata_h5ad(
                        expr_df=expr_df,
                        genes=genes,
                        cell_names=cell_names,
                        series_meta=series_meta,
                        tcr_meta=tcr_meta,
                        out_h5ad=config.out_h5ad,
                    )
                case _ as unreachable:
                    assert_never(unreachable)
        case _ as unreachable:
            assert_never(unreachable)


def main() -> None:
    """CLI entry point for processing GSE123139 into HDF5 format."""
    parser = argparse.ArgumentParser(description="Process GSE123139 raw counts and metadata into AnnData HDF5 (.h5ad) file.")
    _ = parser.add_argument("--input-dir", required=True, help="Directory containing downloaded GSE123139 files.")
    _ = parser.add_argument("--out-h5ad", required=True, help="Target output filepath for AnnData HDF5 (.h5ad) file.")

    args = parser.parse_args()
    config = ProcessConfig(
        input_dir=Path(args.input_dir),
        out_h5ad=Path(args.out_h5ad),
        tcr_file=Some(Path(args.input_dir) / "GSE123139_T_cells_tcrb_v2.txt.gz"),
    )

    match run_pipeline(config):
        case Success(out_path):
            print(f"Successfully processed GSE123139 dataset to HDF5: {out_path}")
            sys.exit(0)
        case Failure(err):
            print(f"Processing error: {err}", file=sys.stderr)
            sys.exit(1)
        case _ as unreachable:
            assert_never(unreachable)


if __name__ == "__main__":
    main()
