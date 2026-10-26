"""Pure functional mapping of internal machine identifiers to publication-grade display labels.

Guarantees:
- Predictor names contain zero underscores ('_') and append publication citations.
- Cohort labels contain zero hyphens ('-') and include author et al., year, n samples, and cancer type.
"""

from __future__ import annotations

from typing import Final, Mapping, Sequence
import polars as pl

PREDICTOR_DISPLAY_MAP: Final[Mapping[str, str]] = {
    "Fehrenbacher_Teff": "Teff (Fehrenbacher et al., 2016)",
    "Gene_CD8_Duo": "CD8 Duo (Chen et al., 2016)",
    "MCP_CD8_T_cells": "CD8 T cells (Becht et al., 2016)",
    "Gene_CD8A": "CD8A (Single Gene)",
    "Davoli_CIS": "CIS (Davoli et al., 2017)",
    "Jiang_CTLs": "CTLs (Jiang et al., 2018)",
    "Ayers_IFNG_6": "IFNG 6 (Ayers et al., 2017)",
    "Gene_CXCL9": "CXCL9 (Single Gene)",
    "Composite_Gene_CXCL9_TMB": "CXCL9 + TMB (Composite)",
    "CYT": "CYT (Rooney et al., 2015)",
    "Cristescu_GEP": "Merck GEP (Cristescu et al., 2018)",
    "TIDE_Inverted": "Inverted TIDE (Jiang et al., 2018)",
    "Ayers_GEP": "Ayers GEP (Ayers et al., 2017)",
    "Wu_MIAS": "MIAS (Wu et al., 2020)",
    "Composite_Ayers_GEP_TMB": "GEP + TMB (Composite)",
    "GeneBio_TargetComposite": "Target Composite (Kong et al., 2022)",
    "TMB_Genomic": "TMB (Genomic)",
    "Roh_IS": "Immune Score (Roh et al., 2017)",
    "Gene_PDCD1": "PDCD1 (Single Gene)",
    "Messina_CKS": "12 Chemokine Score (Messina et al., 2012)",
    "Gene_CTLA4": "CTLA4 (Single Gene)",
    "Gene_CXCL9_AP_Gated": "CXCL9 AP Gated (Composite)",
    "Gene_CD274": "CD274 (Single Gene)",
    "Kong_NetBio_PD1": "NetBio PD1 (Kong et al., 2022)",
    "Freeman_PGM": "PGM (Freeman et al., 2022)",
    "Huang_NRS": "NRS (Huang et al., 2019)",
    "IMPRES": "IMPRES (Auslander et al., 2018)",
    "IPRES_Inverted": "Inverted IPRES (Hugo et al., 2016)",
    "Nurmik_CAFs_Inverted": "Inverted CAFs (Nurmik et al., 2020)",
    "Jiang_TAMs_Inverted": "Inverted TAMs (Jiang et al., 2018)",
    "Jiang_Texh": "T cell Exhaustion (Jiang et al., 2018)",
    "Jiang_Texh_Inverted": "Inverted T cell Exhaustion (Jiang et al., 2018)",
    "Deconvolution_Inference_Our": "Deconvolution Inference (Our)",
    "Cohort_Optimized_Deconv": "Cohort-Optimized Deconvolution (Our Model, 2026)",
}

PREDICTOR_CATEGORY_MAP: Final[Mapping[str, str]] = {
    "MCP_CD8_T_cells": "Signature",
    "TIDE_Inverted": "Signature",
    "Deconvolution_Inference_Our": "Cell Deconvolution",
}

CANONICAL_CANCER_TYPE_MAP: Final[Mapping[str, str]] = {
    "Triple-Negative Breast": "Triple Negative Breast",
    "Triple-Negative Breast Cancer": "Triple Negative Breast",
    "Triple Negative Breast Cancer": "Triple Negative Breast",
    "Triple Negative Breast": "Triple Negative Breast",
    "Breast": "Triple Negative Breast",
    "Clear Cell Renal Cell": "Clear Cell Renal Cell Carcinoma",
    "Clear Cell Renal Cell Carcinoma": "Clear Cell Renal Cell Carcinoma",
    "Renal Cell": "Clear Cell Renal Cell Carcinoma",
    "Pancreatic Ductal": "Pancreatic Ductal",
    "Pancreatic Ductal Adenocarcinoma": "Pancreatic Ductal",
    "Pancreatic": "Pancreatic Ductal",
    "Urothelial Bladder": "Urothelial Bladder",
    "Urothelial Bladder Cancer": "Urothelial Bladder",
    "Bladder": "Urothelial Bladder",
    "Melanoma": "Melanoma",
    "Renal Cell Carcinoma": "Renal Cell Carcinoma",
    "Pan-Cancer": "Pan Cancer",
    "Pan Cancer": "Pan Cancer",
}

COHORT_METADATA_MAP: Final[Mapping[str, Mapping[str, object]]] = {
    "Hugo-iAtlas": {"author": "Hugo et al.", "year": "2016", "cancer": "Melanoma", "n": 27},
    "Riaz-iAtlas": {"author": "Riaz et al.", "year": "2017", "cancer": "Melanoma", "n": 98},
    "Liu-iAtlas": {"author": "Liu et al.", "year": "2019", "cancer": "Melanoma", "n": 122},
    "Gide-iAtlas": {"author": "Gide et al.", "year": "2019", "cancer": "Melanoma", "n": 91},
    "Rosenberg-iAtlas": {"author": "Rosenberg et al.", "year": "2016", "cancer": "Urothelial Bladder", "n": 298},
    "McDermott-iAtlas": {"author": "McDermott et al.", "year": "2018", "cancer": "Clear Cell Renal Cell Carcinoma", "n": 247},
    "Padron-iAtlas": {"author": "Padron et al.", "year": "2022", "cancer": "Pancreatic Ductal", "n": 85},
    "Anders-iAtlas": {"author": "Anders et al.", "year": "2022", "cancer": "Triple Negative Breast", "n": 31},
    "Choueiri-iAtlas": {"author": "Choueiri et al.", "year": "2016", "cancer": "Clear Cell Renal Cell Carcinoma", "n": 16},
    "Melanoma-Combined": {"prefix": "Combined", "cancer": "Melanoma", "n": 338},
    "RenalCell-Combined": {"prefix": "Combined", "cancer": "Renal Cell Carcinoma", "n": 263},
    "PanCancer-Combined": {"prefix": "Combined", "cancer": "Pan Cancer", "n": 1015},
}

CANONICAL_CANCER_ORDER: Final[tuple[str, ...]] = (
    "Clear Cell Renal Cell Carcinoma",
    "Melanoma",
    "Pancreatic Ductal",
    "Triple Negative Breast",
    "Urothelial Bladder",
)


def clean_string_remove_hyphens(text: str) -> str:
    """Replace hyphens and dashes with spaces to ensure zero hyphens in output."""
    return text.replace("-", " ").replace("—", " ").replace("–", " ").strip()


def normalize_cancer_type(raw_cancer: str | None) -> str:
    """Normalize cancer type string to standard canonical representation without hyphens."""
    if not raw_cancer:
        return "Unknown"
    cleaned = clean_string_remove_hyphens(raw_cancer)
    if raw_cancer in CANONICAL_CANCER_TYPE_MAP:
        return CANONICAL_CANCER_TYPE_MAP[raw_cancer]
    if cleaned in CANONICAL_CANCER_TYPE_MAP:
        return CANONICAL_CANCER_TYPE_MAP[cleaned]
    return cleaned


def get_predictor_display_label(raw_name: str) -> str:
    """Return a human-readable predictor name with publication citation and zero underscores."""
    if raw_name in PREDICTOR_DISPLAY_MAP:
        return PREDICTOR_DISPLAY_MAP[raw_name]
    # Fallback: clean underscores
    cleaned = raw_name.replace("_", " ").strip()
    return cleaned


def get_cohort_display_label(
    cohort_id: str,
    cancer_type: str | None = None,
    n_samples: int | None = None,
) -> str:
    """Format a cohort display label with zero hyphens, including author et al., year, n, and cancer type."""
    meta = COHORT_METADATA_MAP.get(cohort_id)
    if meta is not None:
        c_type = (
            normalize_cancer_type(cancer_type)
            if cancer_type
            else str(meta.get("cancer", "Pan Cancer"))
        )
        if "author" in meta and "year" in meta:
            base_name = f"{meta['author']}, {meta['year']}"
        else:
            base_name = f"{meta.get('prefix', 'Combined')} {c_type}"

        canonical_n = n_samples if n_samples is not None and n_samples > 0 else meta.get("n")
        if canonical_n is not None and int(str(canonical_n)) > 0:
            return f"{base_name} (n={canonical_n}, {c_type})"
        elif c_type and c_type not in base_name:
            return f"{base_name} ({c_type})"
        return base_name

    # Generic fallback: remove hyphens
    clean_id = clean_string_remove_hyphens(cohort_id.replace("iAtlas", "").strip())
    c_type = normalize_cancer_type(cancer_type)
    base_name = clean_id

    if n_samples is not None and n_samples > 0:
        return f"{base_name} (n={n_samples}, {c_type})"
    elif c_type and c_type not in base_name:
        return f"{base_name} ({c_type})"
    return base_name


def with_display_labels(
    df: pl.DataFrame,
    cohort_col: str = "cohort_id",
    predictor_col: str = "predictor_name",
    cancer_type_col: str = "cancer_type",
    n_samples_col: str = "n_samples",
) -> pl.DataFrame:
    """Pure transformation enriching a DataFrame with 'cohort_display' and 'predictor_display' columns."""
    res = df

    # 1. Predictor display mapping
    if predictor_col in res.columns:
        res = res.with_columns(
            pl.col(predictor_col)
            .map_elements(get_predictor_display_label, return_dtype=pl.String)
            .alias("predictor_display")
        )

    # 2. Cohort display mapping
    if cohort_col in res.columns:
        has_cancer = cancer_type_col in res.columns
        has_n = n_samples_col in res.columns

        cols_to_select = [cohort_col]
        if has_cancer:
            cols_to_select.append(cancer_type_col)
        if has_n:
            cols_to_select.append(n_samples_col)

        def _map_cohort_cell(s: dict[str, object]) -> str:
            cid = str(s[cohort_col])
            if cid in COHORT_METADATA_MAP:
                return get_cohort_display_label(cid)
            ctype = str(s[cancer_type_col]) if has_cancer and s.get(cancer_type_col) else None
            n_s = int(str(s[n_samples_col])) if has_n and s.get(n_samples_col) is not None else None
            return get_cohort_display_label(cid, ctype, n_s)

        if len(cols_to_select) == 1:
            res = res.with_columns(
                pl.col(cohort_col)
                .map_elements(lambda c: get_cohort_display_label(str(c)), return_dtype=pl.String)
                .alias("cohort_display")
            )
        else:
            res = res.with_columns(
                pl.struct(cols_to_select)
                .map_elements(_map_cohort_cell, return_dtype=pl.String)
                .alias("cohort_display")
            )

    # 3. Category display cleaning (remove underscores, title case, apply overrides)
    if "category" in res.columns:
        if predictor_col in res.columns:
            res = res.with_columns(
                pl.when(pl.col(predictor_col).is_in(list(PREDICTOR_CATEGORY_MAP.keys())))
                .then(pl.col(predictor_col).replace(PREDICTOR_CATEGORY_MAP))
                .otherwise(
                    pl.col("category").cast(pl.String).str.replace_all("_", " ").str.to_titlecase()
                )
                .alias("category")
            )
        else:
            res = res.with_columns(
                pl.col("category").cast(pl.String).str.replace_all("_", " ").str.to_titlecase().alias("category")
            )

    return res


def get_canonical_cohort_order(
    df: pl.DataFrame,
    x_col: str = "cohort_display",
    cancer_order: Sequence[str] | None = None,
) -> list[str]:
    """Return an explicitly ordered list of cohort labels for x-axis heatmap presentation.

    Ordering guarantee:
    1. Single clinical cohorts appear first, ordered by cancer type (according to cancer_order,
       defaulting to CANONICAL_CANCER_ORDER: Clear Cell Renal Cell Carcinoma -> Melanoma -> Pancreatic Ductal ->
       Triple Negative Breast -> Urothelial Bladder).
    2. Within each cancer type, cohorts are sorted chronologically by publication year and author.
    3. Combined standardized pools appear at the end, with disease-specific pools first,
       and the pan-cancer omnibus pool strictly last.
    """
    c_order = tuple(cancer_order) if cancer_order is not None else CANONICAL_CANCER_ORDER
    cancer_rank_map = {c: i for i, c in enumerate(c_order)}

    work_df = (
        with_display_labels(df)
        if x_col == "cohort_display" and "cohort_display" not in df.columns
        else df
    )

    if x_col not in work_df.columns:
        if "cohort_id" in work_df.columns:
            work_df = with_display_labels(work_df)
        else:
            return []

    cols_to_keep = [x_col]
    for c in ["cohort_id", "cancer_type", "pooling_strategy"]:
        if c in work_df.columns and c not in cols_to_keep:
            cols_to_keep.append(c)

    unique_df = work_df.select(cols_to_keep).unique()

    single_cohorts: list[tuple[tuple[int, str, int, str], str]] = []
    combined_cohorts: list[tuple[tuple[int, str, str], str]] = []

    for row in unique_df.iter_rows(named=True):
        label = str(row[x_col])
        cid = str(row.get("cohort_id", ""))
        meta = COHORT_METADATA_MAP.get(cid, {})

        is_combined = (
            row.get("pooling_strategy") == "standardized"
            or cid.endswith("-Combined")
            or "prefix" in meta
            or "Combined" in label
        )

        raw_cancer = (
            str(row.get("cancer_type"))
            if row.get("cancer_type")
            else str(meta.get("cancer", ""))
        )
        norm_cancer = normalize_cancer_type(raw_cancer) if raw_cancer else "Unknown"

        if is_combined:
            is_pancan = (
                1
                if (
                    "PanCancer" in cid
                    or "Pan Cancer" in norm_cancer
                    or "Pan Cancer" in label
                )
                else 0
            )
            combined_cohorts.append(((is_pancan, norm_cancer, label), label))
        else:
            c_rank = cancer_rank_map.get(norm_cancer, 999)
            year = int(str(meta.get("year", 9999))) if "year" in meta else 9999
            author = str(meta.get("author", label))
            single_cohorts.append(((c_rank, norm_cancer, year, author), label))

    single_cohorts.sort(key=lambda x: x[0])
    combined_cohorts.sort(key=lambda x: x[0])

    ordered_labels: list[str] = []
    seen: set[str] = set()
    for _, lbl in single_cohorts + combined_cohorts:
        if lbl not in seen:
            seen.add(lbl)
            ordered_labels.append(lbl)

    return ordered_labels
