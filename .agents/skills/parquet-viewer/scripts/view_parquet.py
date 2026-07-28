import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(description="Inspect a Parquet file.")
    parser.add_argument("file_path", help="Path to the .parquet file")
    parser.add_argument("-n", "--rows", type=int, default=5, help="Number of rows to show for head and tail")
    parser.add_argument("-c", "--columns", type=str, help="Comma-separated list of columns to view")
    args = parser.parse_args()

    if not os.path.exists(args.file_path):
        print(f"Error: File not found: {args.file_path}")
        sys.exit(1)

    try:
        # Import polars locally so the script fails faster if the file doesn't exist
        import polars as pl
    except ImportError:
        print("Error: The 'polars' library is not installed in the current environment.")
        sys.exit(1)

    try:
        df = pl.read_parquet(args.file_path)
    except Exception as e:
        print(f"Error reading Parquet file: {e}")
        sys.exit(1)
        
    if args.columns:
        cols = [c.strip() for c in args.columns.split(",")]
        # Ensure selected columns exist
        missing = [c for c in cols if c not in df.columns]
        if missing:
            print(f"Error: Columns not found: {', '.join(missing)}")
            sys.exit(1)
        df = df.select(cols)

    print(f"========================================")
    print(f"FILE: {args.file_path}")
    print(f"SHAPE: {df.height} rows x {df.width} columns")
    print(f"========================================")
    
    print("\n[ SCHEMA ]")
    schema_dict = df.schema
    for col_name, dtype in schema_dict.items():
        print(f"  {col_name}: {dtype}")
        
    print("\n[ NULL COUNTS ]")
    null_counts = df.select(pl.all().is_null().sum())
    has_nulls = False
    for col in null_counts.columns:
        count = null_counts[col].item(0)
        if count > 0:
            print(f"  {col}: {count} null(s)")
            has_nulls = True
    if not has_nulls:
        print("  No null values found in any column.")

    print(f"\n[ HEAD ({args.rows} rows) ]")
    print(df.head(args.rows))
    
    print(f"\n[ TAIL ({args.rows} rows) ]")
    print(df.tail(args.rows))

if __name__ == "__main__":
    main()
