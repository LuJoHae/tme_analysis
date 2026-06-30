import datalair
from gene_utils import download_from_cbioportal, calculate_maf_tmb
import pandas as pd
from pyensembl import EnsemblRelease
from typing import Callable
from pathlib import Path


class CBioPortalDataset(datalair.Dataset):
    _registry = {
        "VanAllen": "skcm_dfci_2015.tar.gz",
        "Snyder": "skcm_mskcc_2014.tar.gz",
        "Liu": "mel_dfci_2019.tar.gz",
        "Hugo": "mel_ucla_2016.tar.gz",
        "Catalanotti": "skcm_vanderbilt_mskcc_2015.tar.gz",
        "Cloughesy-iAtlas": "gbm_iatlas_prins_2019.tar.gz",
        "Rosenberg-iAtlas": "blca_iatlas_imvigor210_2017.tar.gz",
        "Riaz-iAtlas": "mel_iatlas_riaz_nivolumab_2017.tar.gz",
        "Liu-iAtlas": "mel_iatlas_liu_2019.tar.gz",
        "Gide-iAtlas": "mel_iatlas_gide_2019.tar.gz",
        "Hugo-iAtlas": "mel_iatlas_hugo_ucla_2016.tar.gz",
        "Padron-iAtlas": "paad_iatlas_prince_2022.tar.gz",
        "Choueiri-iAtlas": "ccrcc_iatlas_choueiri_2016.tar.gz",
        "Anders-iAtlas": "brca_iatlas_anders_2022.tar.gz",
        "McDermott-iAtlas": "rcc_iatlas_immotion150_2018.tar.gz"
    }

    def __init__(self, name: str) -> None:
        if name not in self._registry:
            raise ValueError(f"Unmapped dataset identifier: {name}")
        super().__init__(namespace="CBioPortalDataset", dataset_name=name)

    def derive(self, lair: datalair.Lair) -> None:
        output_dir = lair.get_path(self)
        download_from_cbioportal(self._registry[self._dataset_name], output_dir)

    @classmethod
    def get_dataset_name(cls):
        return list(cls._registry.keys())



def transform_expression_to_tpm(counts_df: pd.DataFrame, release: int = 110) -> pd.DataFrame:
    """
    Computes TPM from raw counts using pyensembl for local gene metrics lookup.

    Parameters:
    counts_df: DataFrame of raw counts (rows = Hugo symbols, cols = samples).
    release: Ensembl release version number.
    """
    # 1. Initialize the local Ensembl database manager
    data_source = EnsemblRelease(release)

    def extract_kb_span(symbol: str) -> float | None:
        """Maps a Hugo symbol to its genomic span in kilobases."""
        try:
            genes = data_source.genes_by_name(symbol)
            # Map to the first transcript variant if multiple exist, else return None
            return (genes[0].end - genes[0].start + 1) / 1e3 if genes else None
        except ValueError:
            return None

    # 2. Filter-map operation over the index keyspace using assignment expression
    gene_lengths = pd.Series({
        symbol: length
        for symbol in counts_df.index
        if (length := extract_kb_span(symbol)) is not None
    })

    # 3. Restrict matrix operation to the intersection of existing annotations
    common_genes = counts_df.index.intersection(gene_lengths.index)
    counts = counts_df.loc[common_genes]
    lengths = gene_lengths.loc[common_genes]

    # 4. Matrix Transformations: RPK followed by TPM normalization
    rpk = counts.div(lengths, axis=0)
    return rpk.div(rpk.sum(axis=0), axis=1) * 1e6

def transform_rpkm_to_tpm(rpkm_df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts an RPKM/FPKM normalized DataFrame to TPM.
    """
    # Normalize each column such that its sum is strictly 10^6
    tpm = rpkm_df.div(rpkm_df.sum(axis=0), axis=1) * 1e6

    return tpm

def transform_tpm_to_tpm(df: pd.DataFrame) -> pd.DataFrame: return df


def load_and_transform_data_mrna(data_dir: str | Path) -> pd.DataFrame:
    p_dir = Path(data_dir)
    dispatch_map: dict[str, tuple[Callable[[pd.DataFrame], pd.DataFrame], str]] = {
        "data_mrna_seq_expression.txt": (transform_expression_to_tpm, "none"),
        "data_mrna_seq_tpm.txt": (transform_tpm_to_tpm, "tpm"),
        "data_mrna_seq_rpkm.txt": (transform_rpkm_to_tpm, "rpkm"),
    }

    existing_files = tuple(filter(lambda f: (p_dir / f).is_file(), dispatch_map.keys()))
    if len(existing_files) != 1:
        raise ValueError(
            f"Invariant violation: Expected exactly 1 file, but found {len(existing_files)} in {p_dir}."
        )
    target = existing_files[0]
    transform_fn, normalization = dispatch_map[target]

    data_mrna = transform_fn(pd.read_csv(p_dir / target, sep='\t', index_col=0))
    data_mrna.attrs["original_normalization"] = normalization
    if data_mrna.isna().any(axis=None):
        raise ValueError("Data integrity violation: DataFrame contains NA values.")
    return data_mrna


def load_and_process_data(data_dir):
    data_clinical = load_data_clinical(data_dir)
    data_mrna_seq_expression = load_and_transform_data_mrna(data_dir)
    return data_clinical, data_mrna_seq_expression


def load_data_clinical(data_dir) -> pd.DataFrame:
    # Known aliases for the response column across cBioPortal datasets
    RESPONSE_COL_CANDIDATES = [
    "RESPONSE", "Response", "response",
    "BEST_RESPONSE", "Best_Response", "best_response",
    "RECIST", "recist",
]

    RESPONSE_VALUE_MAP = {
        "Complete Response": "R",
        "Partial Response": "R",
        "Stable Disease": "NR",
        "Progressive Disease": "NR",
        # Add other common encodings as needed:
        "CR": "R",
        "PR": "R",
        "SD": "NR",
        "PD": "NR",
    }

    def _resolve_column(df, candidates):
        """Return the actual column name matching any candidate (case/space-insensitive)."""
        normalized = {c.strip().lower(): c for c in df.columns}
        for cand in candidates:
            key = cand.strip().lower()
            if key in normalized:
                return normalized[key]
        raise KeyError(
            f"None of {candidates} found. Available columns: {list(df.columns)}"
        )

    filepath_clinical_sample = data_dir / "data_clinical_sample.txt"
    data_clinical_sample = pd.read_csv(filepath_clinical_sample, sep="\t", index_col=0, skiprows=4)

    response_col = _resolve_column(data_clinical_sample, RESPONSE_COL_CANDIDATES)
    if response_col is None:
        raise KeyError(
            "No response column found. Available columns: "
            f"{list(data_clinical_sample.columns)}"
        )

    data_clinical_sample = data_clinical_sample.set_index("SAMPLE_ID")
    data_clinical_sample = data_clinical_sample[response_col]


    filepath_mutations = data_dir / "data_mutations.txt"
    if filepath_mutations.exists():
        data_mutations = pd.read_csv(data_dir / "data_mutations.txt", sep="\t", low_memory=False)
        tmb_results = calculate_maf_tmb(data_mutations)
        tmb_results = tmb_results.set_index("Tumor_Sample_Barcode")
        data_clinical_sample = pd.concat([tmb_results, data_clinical_sample], join="inner", axis=1)

    data_clinical_sample["response"] = data_clinical_sample[response_col].map(RESPONSE_VALUE_MAP)
    return data_clinical_sample


def get_dataset_dir(lair, dataset_class, ds_name):
    ds = dataset_class(name=ds_name)
    filepaths = lair.get_dataset_filepaths(ds)
    unpacked_file_key = next(f for f in filepaths.keys() if not f.endswith('.tar.gz'))
    filepath = filepaths[unpacked_file_key]
    return filepath / filepath.name