import os
from pathlib import Path
import polars as pl

DISCOVERED_DATA_DIR = os.environ.get("DISCOVERED_DATA_DIR", config.get("discovered_data_dir", "/storage/halu/data-test"))
REGISTRY_PATH = os.environ.get("REGISTRY_PATH", "data/registry/discovered_solid_tumor_sc_datasets.parquet")
RAW_DIR = f"{DISCOVERED_DATA_DIR}/raw"
PREPROCESSED_DIR = f"{DISCOVERED_DATA_DIR}/preprocessed"


def get_tier0_benchmark_cohorts():
    """Retrieve all Tier 0 Benchmark Core cohort accessions from registry."""
    if not os.path.exists(REGISTRY_PATH):
        return []
    df = pl.read_parquet(REGISTRY_PATH)
    return df.filter(pl.col("tier") == "Tier 0 (Benchmark Core)")["accession"].to_list()


def get_tier1_response_cohorts():
    """Retrieve all Tier 0 and Tier 1 ICB Response cohort accessions from registry."""
    if not os.path.exists(REGISTRY_PATH):
        return []
    df = pl.read_parquet(REGISTRY_PATH)
    return df.filter(pl.col("tier").is_in(["Tier 0 (Benchmark Core)", "Tier 1 (ICB Response)"]))["accession"].to_list()


def get_all_cohorts():
    """Retrieve all 352 discovered cohort accessions from registry."""
    if not os.path.exists(REGISTRY_PATH):
        return []
    df = pl.read_parquet(REGISTRY_PATH)
    return df["accession"].to_list()


TIER0_COHORTS = get_tier0_benchmark_cohorts()
TIER1_COHORTS = get_tier1_response_cohorts()
ALL_COHORTS = get_all_cohorts()


rule process_all_tier0_benchmark_cohorts:
    """Target rule to download and preprocess Tier 0 Premier Benchmark Core cohorts (22 cohorts)."""
    input:
        expand(f"{PREPROCESSED_DIR}/{{cohort}}.h5ad", cohort=TIER0_COHORTS)


rule process_all_discovered_cohorts:
    """Target rule to download and preprocess ALL 352 discovered cohorts across all tiers."""
    input:
        expand(f"{PREPROCESSED_DIR}/{{cohort}}.h5ad", cohort=ALL_COHORTS)


rule process_all_tier1_response_cohorts:
    """Target rule to download and preprocess all Tier 0 + Tier 1 ICB Response cohorts (79 cohorts)."""
    input:
        expand(f"{PREPROCESSED_DIR}/{{cohort}}.h5ad", cohort=TIER1_COHORTS)


rule download_discovered_cohort:
    """Download raw matrix supplementary files or stream H5AD from CZ CELLxGENE in parallel."""
    output:
        directory(f"{RAW_DIR}/{{cohort}}")
    params:
        registry=REGISTRY_PATH
    threads: 2
    resources:
        mem_mb=2000
    shell:
        """
        python scripts/download_cohorts.py \
            --registry {params.registry} \
            --data-raw-dir {RAW_DIR} \
            --cohorts {wildcards.cohort} \
            --max-workers {threads}
        """


rule preprocess_discovered_cohort:
    """Convert raw matrices into standardized dual-layer AnnData H5AD with strict RAM limit."""
    input:
        f"{RAW_DIR}/{{cohort}}"
    output:
        f"{PREPROCESSED_DIR}/{{cohort}}.h5ad"
    params:
        registry=REGISTRY_PATH
    threads: 4
    resources:
        mem_mb=12000
    shell:
        """
        python scripts/batch_preprocess_cohorts.py \
            --registry {params.registry} \
            --data-raw-dir {RAW_DIR} \
            --output-dir {PREPROCESSED_DIR} \
            --cohorts {wildcards.cohort} \
            --max-ram-gb 150.0
        """
