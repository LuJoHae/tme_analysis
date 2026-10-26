"""Unit tests for display labels mapping module.

Guarantees:
- Zero underscores in all mapped predictor names.
- Zero hyphens in all formatted cohort labels.
- Exhaustive coverage of all 31 baseline predictors and 12 cohorts/pools.
"""

from __future__ import annotations

import polars as pl
from tme_response.visualization.display_labels import (
    CANONICAL_CANCER_ORDER,
    COHORT_METADATA_MAP,
    PREDICTOR_DISPLAY_MAP,
    get_canonical_cohort_order,
    get_cohort_display_label,
    get_predictor_display_label,
    with_display_labels,
)

ALL_31_PREDICTORS = [
    "Ayers_GEP",
    "Ayers_IFNG_6",
    "CYT",
    "Composite_Ayers_GEP_TMB",
    "Composite_Gene_CXCL9_TMB",
    "Cristescu_GEP",
    "Davoli_CIS",
    "Fehrenbacher_Teff",
    "Freeman_PGM",
    "GeneBio_TargetComposite",
    "Gene_CD274",
    "Gene_CD8A",
    "Gene_CD8_Duo",
    "Gene_CTLA4",
    "Gene_CXCL9",
    "Gene_CXCL9_AP_Gated",
    "Gene_PDCD1",
    "Huang_NRS",
    "IMPRES",
    "IPRES_Inverted",
    "Jiang_CTLs",
    "Jiang_TAMs_Inverted",
    "Jiang_Texh",
    "Jiang_Texh_Inverted",
    "Kong_NetBio_PD1",
    "MCP_CD8_T_cells",
    "Messina_CKS",
    "Nurmik_CAFs_Inverted",
    "Roh_IS",
    "TIDE_Inverted",
    "TMB_Genomic",
    "Wu_MIAS",
    "Deconvolution_Inference_Our",
]

ALL_12_COHORTS = [
    ("Hugo-iAtlas", "Melanoma", 27),
    ("Riaz-iAtlas", "Melanoma", 98),
    ("Liu-iAtlas", "Melanoma", 122),
    ("Gide-iAtlas", "Melanoma", 91),
    ("Rosenberg-iAtlas", "Bladder", 298),
    ("McDermott-iAtlas", "Clear Cell Renal Cell Carcinoma", 247),
    ("Padron-iAtlas", "Pancreatic", 85),
    ("Anders-iAtlas", "Breast", 31),
    ("Choueiri-iAtlas", "Clear Cell Renal Cell Carcinoma", 16),
    ("Melanoma-Combined", "Melanoma", 338),
    ("RenalCell-Combined", "Clear Cell Renal Cell Carcinoma", 263),
    ("PanCancer-Combined", "Pan Cancer", 1015),
]


def test_all_31_predictors_have_display_labels_without_underscores() -> None:
    """Verify that all 31 predictors are explicitly mapped and contain NO underscores."""
    for pred_id in ALL_31_PREDICTORS:
        assert pred_id in PREDICTOR_DISPLAY_MAP, f"Missing mapping for predictor: {pred_id}"
        label = get_predictor_display_label(pred_id)
        assert "_" not in label, f"Predictor label contains underscore: {label}"
        assert "(" in label and ")" in label, f"Predictor label missing citation / indicator: {label}"


def test_all_cohorts_have_display_labels_without_hyphens() -> None:
    """Verify that all 12 cohorts/pools format to clean strings with NO hyphens."""
    for cid, ctype, n in ALL_12_COHORTS:
        assert cid in COHORT_METADATA_MAP, f"Missing metadata for cohort: {cid}"
        label = get_cohort_display_label(cid, ctype, n)
        assert "-" not in label, f"Cohort label contains hyphen: {label}"
        assert f"n={n}" in label, f"Cohort label missing sample size: {label}"
        assert ctype in label, f"Cohort label missing cancer type: {label}"


def test_with_display_labels_dataframe_transformation() -> None:
    """Verify that with_display_labels enriches a Polars DataFrame correctly."""
    df = pl.DataFrame({
        "cohort_id": ["Gide-iAtlas", "Rosenberg-iAtlas", "PanCancer-Combined"],
        "cancer_type": ["Melanoma", "Urothelial-Bladder", "Pan-Cancer"],
        "predictor_name": ["Fehrenbacher_Teff", "Gene_CD8_Duo", "MCP_CD8_T_cells"],
        "n_samples": [91, 298, 1015],
        "roc_auc": [0.72, 0.65, 0.62],
    })

    enriched = with_display_labels(df)
    assert "cohort_display" in enriched.columns
    assert "predictor_display" in enriched.columns

    for c_disp in enriched["cohort_display"].to_list():
        assert "-" not in c_disp, f"Transformed cohort contains hyphen: {c_disp}"

    for p_disp in enriched["predictor_display"].to_list():
        assert "_" not in p_disp, f"Transformed predictor contains underscore: {p_disp}"

    # Check exact values
    assert enriched["predictor_display"][0] == "Teff (Fehrenbacher et al., 2016)"
    assert enriched["predictor_display"][1] == "CD8 Duo (Chen et al., 2016)"
    assert enriched["predictor_display"][2] == "CD8 T cells (Becht et al., 2016)"
    assert enriched["cohort_display"][0] == "Gide et al., 2019 (n=91, Melanoma)"


def test_cancer_type_synonyms_produce_identical_cohort_labels() -> None:
    """Assert that minor variations in cancer type strings resolve to the identical canonical label."""
    synonym_pairs = [
        ("Anders-iAtlas", "Triple-Negative Breast", "Triple-Negative Breast Cancer"),
        ("Choueiri-iAtlas", "Clear Cell Renal Cell", "Clear Cell Renal Cell Carcinoma"),
        ("Padron-iAtlas", "Pancreatic Ductal", "Pancreatic Ductal Adenocarcinoma"),
        ("Rosenberg-iAtlas", "Urothelial Bladder", "Urothelial Bladder Cancer"),
        ("PanCancer-Combined", "Pan-Cancer", "Pan Cancer"),
    ]
    for cid, name1, name2 in synonym_pairs:
        label1 = get_cohort_display_label(cid, cancer_type=name1)
        label2 = get_cohort_display_label(cid, cancer_type=name2)
        assert label1 == label2, f"Synonym mismatch for {cid}: '{label1}' != '{label2}'"


def test_with_display_labels_prevents_split_columns() -> None:
    """Assert that mixed cancer type strings and sample sizes do not split a cohort into multiple columns."""
    df = pl.DataFrame({
        "cohort_id": [
            "Anders-iAtlas",
            "Anders-iAtlas",
            "Anders-iAtlas",
            "Rosenberg-iAtlas",
            "Rosenberg-iAtlas",
        ],
        "cancer_type": [
            "Triple-Negative Breast",
            "Triple-Negative Breast Cancer",
            "Triple Negative Breast",
            "Urothelial Bladder",
            "Urothelial Bladder Cancer",
        ],
        "n_samples": [31, 31, 25, 298, 214],
        "predictor_name": ["Fehrenbacher_Teff", "Deconvolution_Inference_Our", "TMB_Genomic", "Gene_CD8A", "TMB_Genomic"],
    })

    enriched = with_display_labels(df)
    unique_anders = enriched.filter(pl.col("cohort_id") == "Anders-iAtlas")["cohort_display"].unique().to_list()
    assert len(unique_anders) == 1, f"Anders split into multiple columns: {unique_anders}"

    unique_rosenberg = enriched.filter(pl.col("cohort_id") == "Rosenberg-iAtlas")["cohort_display"].unique().to_list()
    assert len(unique_rosenberg) == 1, f"Rosenberg split into multiple columns: {unique_rosenberg}"


def test_predictor_category_overrides() -> None:
    """Verify that MCP CD8 T cells and TIDE are in Signature, and Deconvolution Inference is in Cell Deconvolution."""
    df = pl.DataFrame({
        "predictor_name": ["MCP_CD8_T_cells", "TIDE_Inverted", "Deconvolution_Inference_Our", "Gene_CD8A"],
        "category": ["cell_deconvolution", "systems_model", "signature", "single_gene"],
    })
    enriched = with_display_labels(df)
    cats = dict(zip(enriched["predictor_name"].to_list(), enriched["category"].to_list()))
    assert cats["MCP_CD8_T_cells"] == "Signature"
    assert cats["TIDE_Inverted"] == "Signature"
    assert cats["Deconvolution_Inference_Our"] == "Cell Deconvolution"
    assert cats["Gene_CD8A"] == "Single Gene"


def test_canonical_cohort_order() -> None:
    """Verify that get_canonical_cohort_order groups single cohorts by cancer type then combined cohorts."""
    all_cids = [
        "Hugo-iAtlas",
        "Riaz-iAtlas",
        "Liu-iAtlas",
        "Gide-iAtlas",
        "Rosenberg-iAtlas",
        "McDermott-iAtlas",
        "Padron-iAtlas",
        "Anders-iAtlas",
        "Choueiri-iAtlas",
        "Melanoma-Combined",
        "RenalCell-Combined",
        "PanCancer-Combined",
    ]
    df = pl.DataFrame({"cohort_id": all_cids})

    ordered_labels = get_canonical_cohort_order(df, x_col="cohort_display")
    assert len(ordered_labels) == 12

    # Single cohorts (indices 0..8) must come before combined pools (indices 9..11)
    for lbl in ordered_labels[:9]:
        assert not lbl.startswith("Combined"), f"Combined pool found in single cohort range: {lbl}"
    for lbl in ordered_labels[9:]:
        assert lbl.startswith("Combined"), f"Single cohort found in combined pool range: {lbl}"

    # Clear Cell Renal Cell Carcinoma (Choueiri, McDermott) must be indices 0 and 1
    assert "Choueiri" in ordered_labels[0]
    assert "McDermott" in ordered_labels[1]
    assert "Clear Cell Renal Cell Carcinoma" in ordered_labels[0]
    assert "Clear Cell Renal Cell Carcinoma" in ordered_labels[1]

    # Melanoma (Hugo, Riaz, Gide, Liu) must be indices 2, 3, 4, 5 in chronological order
    assert "Hugo" in ordered_labels[2]
    assert "Riaz" in ordered_labels[3]
    assert "Gide" in ordered_labels[4]
    assert "Liu" in ordered_labels[5]

    # Pancreatic, Breast, Bladder
    assert "Padron" in ordered_labels[6]
    assert "Anders" in ordered_labels[7]
    assert "Rosenberg" in ordered_labels[8]

    # Combined pools: Melanoma, Renal Cell, Pan Cancer last
    assert "Melanoma" in ordered_labels[9]
    assert "Renal Cell" in ordered_labels[10]
    assert "Pan Cancer" in ordered_labels[11]

    # Also verify with cohort_id
    ordered_cids = get_canonical_cohort_order(df, x_col="cohort_id")
    assert ordered_cids == [
        "Choueiri-iAtlas",
        "McDermott-iAtlas",
        "Hugo-iAtlas",
        "Riaz-iAtlas",
        "Gide-iAtlas",
        "Liu-iAtlas",
        "Padron-iAtlas",
        "Anders-iAtlas",
        "Rosenberg-iAtlas",
        "Melanoma-Combined",
        "RenalCell-Combined",
        "PanCancer-Combined",
    ]


def test_create_cancer_group_track_and_heatmap_grouping() -> None:
    """Verify that create_cancer_group_track produces valid Altair chart and create_auc_heatmap succeeds."""
    from returns.result import Success
    from tme_response.visualization.heatmap import create_auc_heatmap, create_cancer_group_track

    df = pl.DataFrame({
        "cohort_id": [cid for cid, _, _ in ALL_12_COHORTS],
        "cancer_type": [ctype for _, ctype, _ in ALL_12_COHORTS],
        "predictor_name": ["Fehrenbacher_Teff"] * len(ALL_12_COHORTS),
        "n_samples": [n for _, _, n in ALL_12_COHORTS],
        "roc_auc": [0.70] * len(ALL_12_COHORTS),
        "p_value_vs_half": [0.01] * len(ALL_12_COHORTS),
    })

    enriched = with_display_labels(df)
    order = get_canonical_cohort_order(enriched, x_col="cohort_display")

    # Test track creation
    track = create_cancer_group_track(cohort_order=order, df=enriched, x_col="cohort_display")
    assert track is not None

    # Test heatmap with cancer groups
    res = create_auc_heatmap(enriched, metric="roc_auc", show_cancer_groups=True)
    assert isinstance(res, Success)


def test_renal_cell_carcinoma_full_string() -> None:
    """Verify that renal cohorts resolve to the full 'Clear Cell Renal Cell Carcinoma' string."""
    for cid in ["Choueiri-iAtlas", "McDermott-iAtlas"]:
        lbl = get_cohort_display_label(cid)
        assert "Clear Cell Renal Cell Carcinoma" in lbl
        assert not lbl.endswith("Clear Cell Renal Cell)")





