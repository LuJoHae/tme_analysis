---
name: parquet-viewer
description: Skill for cleanly viewing the schema, data, and null-counts of Parquet files using a robust helper script.
---

# Parquet Viewer Skill

When the user asks you to view, inspect, or summarize a `.parquet` file, you should use the provided helper script rather than writing ad-hoc python code in the terminal.

## Usage

You can use the `run_command` tool to execute the helper script:

```bash
python .agents/skills/parquet-viewer/scripts/view_parquet.py /path/to/data.parquet
```

### Options
- `-n`, `--rows`: Number of rows to show from both the **head** and **tail** of the file (default: 5)
- `-c`, `--columns`: Comma-separated list of columns to view (useful if the file is extremely wide and you only want to see a specific subset)

## Behavior
The script uses `polars` to load the dataset and will output:
1. The **shape** of the dataset (total rows x total columns).
2. The **schema** (a list of all column names and their data types).
3. **Null counts** for each column (if any nulls exist, it lists them; otherwise it confirms there are zero nulls).
4. The first N rows (`[ HEAD ]`).
5. The last N rows (`[ TAIL ]`).
