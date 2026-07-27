# Snakefile
# Orchestrates data processing on the remote server

# TODO: Update this path to the absolute path of your data directory on 'olm'
DATA_DIR = "/TODO/UPDATE/THIS/PATH/TO/DATA"
RESULTS_DIR = f"{DATA_DIR}/results"

# The 'all' rule defines what files should ultimately be generated.
# Snakemake will figure out the dependencies to build these.
rule all:
    input:
        f"{RESULTS_DIR}/example_output.csv"

# Example rule: This rule will be triggered to build 'example_output.csv'
# Replace this with your actual scripts and data files.
rule run_example_script:
    input:
        # Example input data file (you don't strictly need to define one if it's generated dynamically, but it's good practice)
        # raw_data = f"{DATA_DIR}/raw_data.csv"
    output:
        # The script must generate this file
        summary = f"{RESULTS_DIR}/example_output.csv"
    script:
        # Snakemake will execute this script and inject the `snakemake` object
        # so the script can access `snakemake.output.summary`
        "scripts/example_script.py"
