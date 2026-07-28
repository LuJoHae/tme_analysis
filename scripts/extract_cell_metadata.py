from typing import Literal, assert_never
import polars as pl
import re
import sys
import argparse
from pathlib import Path
from returns.result import Result, Success, Failure

# Regex definition with named capture groups
REGEX = r"^(?P<row>[A-Za-z])(?P<col>[1-9]|1[0-9]|2[0-4])_P(?P<plate>\d{1,2})_(?P<patient>[Mm]\d{1,2}|MMD\d+(?:-\d+[A-Z])?|M\d{2}-\d{1,2}-\d{1,2}-\d{2})(?:-B(?P<biopsy>\d+))?(?:_L(?P<lane>\d{3}))?(?:_(?P<enrichment>T|myeloid)_enriched)?$"
PATTERN = re.compile(REGEX)

SCHEMA = {
    "cell_id": pl.String,
    "plate-row": pl.String,
    "plate-col": pl.String,
    "plate": pl.String,
    "patient": pl.String,
    "patient_id_category": pl.String,
    "biopsy": pl.String,
    "sequencing_lane": pl.String,
    "enrichment": pl.String,
}

def determine_patient_category(patient: str) -> str:
    if patient.startswith("MMD"):
        return "MMD"
    elif "-" in patient:
        return "M_with_date"
    else:
        return "normal_M"

def parse_cell_id(cell_id: str) -> Result[dict[str, str | None], str]:
    match = PATTERN.match(cell_id)
    if not match:
        return Failure(f"Cell ID does not match pattern: {cell_id}")
    
    d = match.groupdict()
    
    # Format with leading zeros
    row = d["row"].upper()
    col = str(int(d["col"])).zfill(2)
    plate = str(int(d["plate"])).zfill(2)
    patient = d["patient"]
    patient_category = determine_patient_category(patient)
    
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
        "patient": patient,
        "patient_id_category": patient_category,
        "biopsy": biopsy,
        "sequencing_lane": lane,
        "enrichment": enrichment
    })

def process_metadata(meta_path: Path) -> Result[pl.DataFrame, str]:
    try:
        df = pl.read_parquet(meta_path)
        cell_ids = df["title"].to_list()
        
        parsed_records = []
        for cid in cell_ids:
            res = parse_cell_id(cid)
            match res:
                case Success(record):
                    parsed_records.append(record)
                case Failure(err):
                    return Failure(err)
                    
        parsed_df = pl.DataFrame(parsed_records, schema=SCHEMA)
        return Success(parsed_df)
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
