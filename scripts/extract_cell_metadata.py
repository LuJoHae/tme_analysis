import argparse
import polars as pl
import re
import sys
from pathlib import Path
from returns.result import Result, Success, Failure

# Regex definition with named capture groups (patient renamed to melanoma_sample)
REGEX = r"^(?P<row>[A-Za-z])(?P<col>[1-9]|1[0-9]|2[0-4])_P(?P<plate>\d{1,2})_(?P<melanoma_sample>[Mm]\d{1,2}|MMD\d+(?:-\d+[A-Z])?|M\d{2}-\d{1,2}-\d{1,2}-\d{2})(?:-B(?P<biopsy>\d+))?(?:_L(?P<lane>\d{3}))?(?:_(?P<enrichment>T|myeloid)_enriched)?$"
PATTERN = re.compile(REGEX)

SCHEMA = {
    "cell_id": pl.String,
    "plate-row": pl.String,
    "plate-col": pl.String,
    "plate": pl.String,
    "melanoma-sample": pl.String,
    "sample_id_category": pl.String,
    "biopsy": pl.String,
    "sequencing_lane": pl.String,
    "enrichment": pl.String,
}

def determine_sample_category(sample: str) -> str:
    if sample.startswith("MMD"):
        return "MMD"
    elif "-" in sample:
        return "M_with_date"
    else:
        return "normal_M"

def parse_cell_id(cell_id: str) -> Result[dict[str, str | None], str]:
    match = PATTERN.match(cell_id)
    if not match:
        return Failure(f"Cell ID does not match pattern: {cell_id}")
    
    d = match.groupdict()
    
    row = d["row"].upper()
    col = str(int(d["col"])).zfill(2)
    plate = str(int(d["plate"])).zfill(2)
    sample = d["melanoma_sample"]
    sample_category = determine_sample_category(sample)
    
    biopsy = d["biopsy"]
    if biopsy is not None:
        biopsy = str(int(biopsy)).zfill(2)
        
    lane = d["lane"]
    if lane is not None:
        lane = str(int(lane)).zfill(3)
        
    enrichment = d["enrichment"]
    
    return Success({
        "cell_id": cell_id,
        "plate-row": row,
        "plate-col": col,
        "plate": plate,
        "melanoma-sample": sample,
        "sample_id_category": sample_category,
        "biopsy": biopsy,
        "sequencing_lane": lane,
        "enrichment": enrichment
    })

def process_metadata(meta_path: Path) -> Result[pl.DataFrame, str]:
    try:
        # Read the original metadata
        meta_df = pl.read_parquet(meta_path)
        cell_ids = meta_df["title"].to_list()
        
        parsed_records = []
        for cid in cell_ids:
            res = parse_cell_id(cid)
            match res:
                case Success(record):
                    parsed_records.append(record)
                case Failure(err):
                    return Failure(err)
                    
        parsed_df = pl.DataFrame(parsed_records, schema=SCHEMA)
        
        # Join with metadata
        pat_col = "characteristics: patinet ID (Pre=baseline; Post= on treatment)"
        joined_df = parsed_df.join(meta_df, left_on="cell_id", right_on="title")
        
        # Extract and split patient info
        final_df = joined_df.with_columns(
            pl.col(pat_col).str.split_exact("_", 1).struct.rename_fields(["treatment_status", "patient_id"]).alias("split_pat")
        ).unnest("split_pat").select([
            "cell_id", "plate-row", "plate-col", "plate", "melanoma-sample", "sample_id_category",
            "treatment_status", "patient_id", "biopsy", "sequencing_lane", "enrichment"
        ])
        
        initial_count = final_df.height
        
        # 1. Exclude the cell with a small 'm' in melanoma-sample
        final_df = final_df.filter(~pl.col("melanoma-sample").str.starts_with("m"))
        
        # 2. Exclude cells where the sample ID maps to a minority patient ID (typos)
        # Find the dominant patient_id for each melanoma-sample
        dominant_mapping = (
            final_df
            .group_by(["melanoma-sample", "patient_id"])
            .len()
            .sort("len", descending=True)
            .group_by("melanoma-sample")
            .first()
            .select(["melanoma-sample", "patient_id"])
        )
        
        # Keep only the rows that match the dominant mapping
        final_df = final_df.join(dominant_mapping, on=["melanoma-sample", "patient_id"], how="inner")
        
        final_count = final_df.height
        dropped_count = initial_count - final_count
        
        # Assert that these rules only excluded a very small number of cells (e.g., < 10)
        if dropped_count > 10:
            return Failure(f"Assertion failed: Too many cells excluded ({dropped_count} cells). Expected < 10.")
            
        print(f"Data cleaning successful. Excluded {dropped_count} mislabeled/typo cells.")
        
        return Success(final_df)
    except Exception as e:
        return Failure(f"Failed to process metadata: {str(e)}")

def save_dataframe(df: pl.DataFrame, out_path: Path) -> Result[Path, str]:
    try:
        # Save explicitly with string types inferred naturally from the dicts
        df.write_parquet(out_path)
        return Success(out_path)
    except Exception as e:
        return Failure(f"Failed to save parquet: {str(e)}")

def run_extraction(meta_path: Path, out_path: Path) -> Result[Path, str]:
    return process_metadata(meta_path).bind(
        lambda df: save_dataframe(df, out_path)
    )

def main() -> None:
    parser = argparse.ArgumentParser(description="Extract metadata from GSE120575 cell IDs")
    parser.add_argument("--meta", required=True, help="Input metadata parquet file")
    parser.add_argument("--out", required=True, help="Output parsed cell IDs parquet file")
    args = parser.parse_args()

    meta_path = Path(args.meta)
    out_path = Path(args.out)
    
    print(f"Extracting cell ID metadata from {meta_path}...")
    match run_extraction(meta_path, out_path):
        case Success(path):
            print(f"Successfully saved parsed metadata to {path}")
            sys.exit(0)
        case Failure(err):
            print(f"Error: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
