# Split Computation and Plotting

Always separate heavy computations (like model fitting, data aggregation, and statistical testing) from data visualization.

1. **Computation Scripts**: Scripts performing calculations or model fitting must save their final results (e.g., DataFrames, model summaries) to disk in a standard format (like `.csv`, `.parquet`, `.h5ad`, or `.json`). They should *not* generate plots.
2. **Plotting Scripts**: Scripts responsible for visualization must load the pre-computed results from disk and generate the plots. They should *not* perform heavy calculations.

This pattern ensures that visualizations can be rapidly iterated on and tweaked (e.g., fixing label rotations or colors) without having to re-run expensive computational pipelines.

3. **Output Directories**: Unless explicitly specified otherwise, always route raw data files (`.csv`, `.parquet`, `.h5ad`) to the `output/` directory, and strictly route generated plots (`.svg`, `.png`) to the `results/` directory.
