import urllib.request
import os
from pathlib import Path
import polars as pl
from returns.result import Result, Success, Failure
import gzip

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl"
FILES = [
    "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
    "GSE120575_patient_ID_single_cells.txt.gz"
]

def _download_file(url: str, dest: Path) -> Result[Path, str]:
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            urllib.request.urlretrieve(url, dest)
        return Success(dest)
    except Exception as e:
        return Failure(f"Failed to download {url}: {str(e)}")

def download_gse120575(dest_dir: str) -> Result[tuple[Path, Path], str]:
    dest_path = Path(dest_dir)
    
    tpm_url = f"{BASE_URL}/{FILES[0]}"
    meta_url = f"{BASE_URL}/{FILES[1]}"
    
    return _download_file(tpm_url, dest_path / FILES[0]).bind(
        lambda tpm_file: _download_file(meta_url, dest_path / FILES[1]).map(
            lambda meta_file: (tpm_file, meta_file)
        )
    )

def _load_and_join_data(tpm_path: Path, meta_path: Path) -> Result[pl.DataFrame, str]:
    try:
        # Load TPM data (rows are genes, columns are cells)
        # Assuming the first column is gene names, we might need to transpose it.
        # But first let's just load it. The prompt says "save it as hdf files that can be read with polars" (amended to parquet).
        # We will save the expression matrix and metadata separately, or joined.
        # Given it's scRNA-seq (16000 cells x 20000 genes), doing a full join in polars might be huge if transposed.
        # Let's just save them as two parquet files for efficiency.
        df_tpm = pl.read_csv(tpm_path, separator='\t', null_values=["NA"])
        df_meta = pl.read_csv(meta_path, separator='\t', skip_rows=19) # GEO matrices often have a header, we need to be careful. Wait, patient ID might just be standard tab separated.
        
        # Actually, let's just convert the raw downloaded files directly to parquet without complex joining here.
        # We can do joining in the analysis script if needed.
        return Success((df_tpm, df_meta))
    except Exception as e:
        return Failure(f"Failed to load data: {str(e)}")

def _save_parquet(df: pl.DataFrame, out_path: Path) -> Result[Path, str]:
    try:
        df.write_parquet(out_path)
        return Success(out_path)
    except Exception as e:
        return Failure(f"Failed to save parquet to {out_path}: {str(e)}")

def process_to_parquet(tpm_path: Path, meta_path: Path, out_dir: str) -> Result[tuple[Path, Path, Path], str]:
    out_dir_path = Path(out_dir)
    out_dir_path.mkdir(parents=True, exist_ok=True)
    
    tpm_out = out_dir_path / "gse120575_tpm.parquet"
    tpm_meta_out = out_dir_path / "gse120575_tpm_cell_metadata.parquet"
    meta_out = out_dir_path / "gse120575_meta.parquet"
    
    # 1. Read the first two lines of the TPM file to separate column names from the metadata row
    try:
        try:
            with gzip.open(tpm_path, 'rt', encoding='utf-8', errors='replace') as f:
                line1 = f.readline().strip('\r\n')
                line2 = f.readline().strip('\r\n')
        except (gzip.BadGzipFile, UnicodeDecodeError):
            with open(tpm_path, 'rt', encoding='utf-8', errors='replace') as f:
                line1 = f.readline().strip('\r\n')
                line2 = f.readline().strip('\r\n')
    except Exception as e:
        return Failure(f"Failed to read TPM header: {str(e)}")

    cols = line1.split('\t')
    cols[0] = "gene" # Replace the empty first column with 'gene'
    
    meta_row = line2.split('\t')
    
    # Handle potential trailing tab (empty column name) in TPM file
    if len(cols) > 1 and not cols[-1].strip():
        cols[-1] = "__drop_me__"
        meta_row[-1] = ""

    # 2. Save the TPM metadata (Cell_ID -> Patient_ID)
    try:
        cell_ids = cols[1:]
        patient_ids = meta_row[1:]
        if cols[-1] == "__drop_me__":
            cell_ids = cell_ids[:-1]
            patient_ids = patient_ids[:-1]
            
        df_tpm_meta = pl.DataFrame({"Cell_ID": cell_ids, "Patient_ID": patient_ids})
        df_tpm_meta.write_parquet(tpm_meta_out)
    except Exception as e:
        return Failure(f"Failed to save TPM metadata: {str(e)}")

    # 3. Process the rest of the TPM matrix and patient metadata
    try:
        # scan TPM data, skipping the 2 header rows
        (
            pl.scan_csv(
                tpm_path, 
                separator='\t', 
                has_header=False,
                skip_rows=2, 
                new_columns=cols,
                truncate_ragged_lines=True,
                encoding="utf8-lossy"
            )
            .select([
                pl.col("gene").cast(pl.String),
                pl.all().exclude(["gene", "__drop_me__"]).cast(pl.Float32, strict=True)
            ])
            .sink_parquet(tpm_out)
        )
        
        # scan patient metadata, skipping 19 rows, take 7 cols, ignore footer
        (
            pl.scan_csv(
                meta_path, 
                separator='\t', 
                skip_rows=19,
                truncate_ragged_lines=True,
                encoding="utf8-lossy"
            )
            .select(pl.col("*").head(7)) # Only take the first 7 columns if there are trailing empty ones
            .filter(pl.col("Sample name").is_not_null() & pl.col("Sample name").str.starts_with("Sample"))
            .sink_parquet(meta_out)
        )
        
        # Assert no nulls in the TPM data
        null_counts = pl.scan_parquet(tpm_out).select(pl.all().is_null().sum()).collect()
        if null_counts.sum_horizontal().item(0) > 0:
            return Failure("Assertion failed: Null values found in TPM dataset!")
            
    except Exception as e:
        try:
            # Fallback to read_csv
            (
                pl.read_csv(
                    tpm_path, 
                    separator='\t', 
                    has_header=False,
                    skip_rows=2,
                    new_columns=cols,
                    truncate_ragged_lines=True, 
                    encoding="utf8-lossy"
                )
                .select([
                    pl.col("gene").cast(pl.String),
                    pl.all().exclude(["gene", "__drop_me__"]).cast(pl.Float32, strict=True)
                ])
                .write_parquet(tpm_out)
            )
            
            (
                pl.read_csv(
                    meta_path, 
                    separator='\t', 
                    skip_rows=19,
                    truncate_ragged_lines=True, 
                    encoding="utf8-lossy"
                )
                .select(pl.col("*").head(7))
                .filter(pl.col("Sample name").is_not_null() & pl.col("Sample name").str.starts_with("Sample"))
                .write_parquet(meta_out)
            )
            
            # Assert no nulls
            null_counts = pl.scan_parquet(tpm_out).select(pl.all().is_null().sum()).collect()
            if null_counts.sum_horizontal().item(0) > 0:
                return Failure("Assertion failed: Null values found in TPM dataset!")
                
        except Exception as e2:
            return Failure(f"Failed to convert to parquet: {str(e2)}")
            
    return Success((tpm_out, tpm_meta_out, meta_out))

def fetch_and_format_gse120575(dest_dir: str) -> Result[tuple[Path, Path, Path], str]:
    return download_gse120575(dest_dir).bind(
        lambda paths: process_to_parquet(paths[0], paths[1], dest_dir)
    )
