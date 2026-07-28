import polars as pl
import re

def main():
    meta_path = "/storage/halu/data/GSE120575/gse120575_meta.parquet"
    print(f"Loading metadata from {meta_path}...")
    df = pl.read_parquet(meta_path)
    
    cell_ids = df["title"].to_list()
    
    regex = r"^[A-P](?:[1-9]|1[0-9]|2[0-4])_P\d_M\d{2}(?:-B\d+)?(?:_L\d{3})?(?:_T_enriched)?$"
    pattern = re.compile(regex)
    
    valid_count = 0
    invalid_count = 0
    invalid_examples = []
    
    for cid in cell_ids:
        if pattern.match(cid):
            valid_count += 1
        else:
            invalid_count += 1
            if len(invalid_examples) < 10:
                invalid_examples.append(cid)
                
    print(f"Total cells checked: {len(cell_ids)}")
    print(f"Valid IDs (matched regex): {valid_count}")
    print(f"Invalid IDs (did not match): {invalid_count}")
    
    if invalid_count > 0:
        print("\nExamples of invalid IDs:")
        for ex in invalid_examples:
            print(ex)

if __name__ == "__main__":
    main()
