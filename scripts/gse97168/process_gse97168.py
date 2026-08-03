"""Raw processing script to convert GSE97168 UMI count table and metadata into AnnData HDF5 (.h5ad)."""

import argparse
import gzip
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
    """Immutable data processing configuration for GSE97168."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    input_dir: Path
    out_h5ad: Path


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


def load_umi_table(umitab_path: Path) -> Result[tuple[pl.DataFrame, list[str], list[str]], str]:
    """Load GSE97168 UMI count table using Polars."""
    try:
        print(f"Loading UMI table from {umitab_path}...")
        open_fn = gzip.open if umitab_path.suffix == ".gz" else open
        with open_fn(umitab_path, "rt", encoding="utf-8", errors="replace") as f:
            line1 = f.readline().strip("\r\n")

        cell_names = [c.strip('"') for c in line1.split("\t") if c.strip('"')]
        cols = ["gene"] + cell_names

        df = pl.read_csv(
            umitab_path,
            separator="\t",
            has_header=False,
            skip_rows=1,
            new_columns=cols,
            truncate_ragged_lines=True,
            encoding="utf8-lossy",
        )

        df = df.with_columns(pl.col("gene").str.strip_chars('"'))
        genes = df["gene"].to_list()
        expr_df = df.drop("gene")

        return Success((expr_df, genes, cell_names))
    except Exception as err:
        return Failure(f"Failed to load UMI count table {umitab_path}: {err}")


def load_metadata_table(meta_path: Path) -> Result[pl.DataFrame, str]:
    """Load cell metadata table using Polars, skipping header description comments."""
    try:
        skip_n = 0
        open_fn = gzip.open if meta_path.suffix == ".gz" else open
        with open_fn(meta_path, "rt", encoding="utf-8", errors="replace") as f:
            for idx, line in enumerate(f):
                parts = [p.strip('"').strip() for p in line.split("\t")]
                if parts and parts[0] == "well" and len(parts) >= 2 and parts[1] == "Amp_batch_ID":
                    skip_n = idx
                    break

        meta_df = pl.read_csv(
            meta_path,
            separator="\t",
            has_header=True,
            skip_rows=skip_n,
            truncate_ragged_lines=True,
            encoding="utf8-lossy",
        )
        return Success(meta_df)
    except Exception as err:
        return Failure(f"Failed loading metadata table {meta_path}: {err}")


def build_anndata_h5ad(
    expr_df: pl.DataFrame,
    genes: list[str],
    cell_names: list[str],
    suppl_meta: pl.DataFrame,
    series_meta: Maybe[pl.DataFrame],
    out_h5ad: Path,
) -> Result[Path, str]:
    """Construct AnnData object and serialize to HDF5 (.h5ad) format."""
    try:
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)

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

        # Merge supplementary metadata
        suppl_pd = suppl_meta.to_pandas()
        suppl_pd = suppl_pd.loc[:, [c for c in suppl_pd.columns if c and not c.startswith("_duplicated")]]

        id_col = suppl_pd.columns[0]
        for col in suppl_pd.columns:
            if col.lower() in ("cell_id", "cell_name", "cell", "sample", "well"):
                id_col = col
                break

        suppl_pd[id_col] = suppl_pd[id_col].astype(str)
        suppl_pd = suppl_pd.drop_duplicates(subset=[id_col]).set_index(id_col)
        obs_df = obs_df.join(suppl_pd, how="left")

        # Deduplicate columns if any duplicate names remain
        obs_df = obs_df.loc[:, ~obs_df.columns.duplicated()]

        # Optionally join series matrix metadata if available
        match series_meta:
            case Some(s_df):
                s_pd = s_df.to_pandas()
                if "sample_id" in s_pd.columns:
                    s_pd["sample_id"] = s_pd["sample_id"].astype(str)
                    s_pd = s_pd.set_index("sample_id")
                    if "sample_id" in obs_df.columns:
                        obs_df = obs_df.join(s_pd, on="sample_id", how="left")
            case _:
                pass

        # Sanitize metadata columns to string type
        for col in obs_df.columns:
            obs_df[col] = obs_df[col].astype(str).replace("nan", "")

        adata = ad.AnnData(X=X_sparse, obs=obs_df, var=var_df)

        print(f"Saving raw AnnData object ({adata.n_obs} cells x {adata.n_vars} genes) to {out_h5ad}...")
        adata.write_h5ad(str(out_h5ad), compression="gzip")

        return Success(out_h5ad)
    except Exception as err:
        return Failure(f"Failed creating AnnData HDF5 file: {err}")


def run_pipeline(config: ProcessConfig) -> Result[Path, str]:
    """Execute functional processing pipeline for GSE97168 raw AnnData HDF5."""
    umitab_file = config.input_dir / "GSE97168_umitab.txt.gz"
    meta_file = config.input_dir / "GSE97168_metadata.txt.gz"
    series_file = config.input_dir / "GSE97168_series_matrix.txt.gz"

    match load_umi_table(umitab_file):
        case Failure(err):
            return Failure(err)
        case Success((expr_df, genes, cell_names)):
            match load_metadata_table(meta_file):
                case Failure(err):
                    return Failure(err)
                case Success(suppl_meta):
                    series_meta: Maybe[pl.DataFrame] = Nothing
                    if series_file.exists():
                        match parse_series_matrix_metadata(series_file):
                            case Success(df):
                                series_meta = Some(df)
                            case Failure(_):
                                series_meta = Nothing
                            case _ as unreachable:
                                assert_never(unreachable)

                    return build_anndata_h5ad(
                        expr_df=expr_df,
                        genes=genes,
                        cell_names=cell_names,
                        suppl_meta=suppl_meta,
                        series_meta=series_meta,
                        out_h5ad=config.out_h5ad,
                    )
                case _ as unreachable:
                    assert_never(unreachable)
        case _ as unreachable:
            assert_never(unreachable)


def main() -> None:
    """CLI entry point for processing GSE97168 into raw AnnData HDF5 file."""
    parser = argparse.ArgumentParser(description="Process GSE97168 UMI count table and metadata into AnnData HDF5 (.h5ad).")
    _ = parser.add_argument("--input-dir", required=True, help="Directory containing downloaded GSE97168 files.")
    _ = parser.add_argument("--out-h5ad", required=True, help="Target output filepath for AnnData HDF5 file.")

    args = parser.parse_args()
    config = ProcessConfig(input_dir=Path(args.input_dir), out_h5ad=Path(args.out_h5ad))

    match run_pipeline(config):
        case Success(out_path):
            print(f"Successfully processed GSE97168 dataset to HDF5: {out_path}")
            sys.exit(0)
        case Failure(err):
            print(f"Processing error: {err}", file=sys.stderr)
            sys.exit(1)
        case _ as unreachable:
            assert_never(unreachable)


if __name__ == "__main__":
    main()
