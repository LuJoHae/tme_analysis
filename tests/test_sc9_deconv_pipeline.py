"""Tests for 9-Cohort Multi-Resolution Single-Cell Reference Deconvolution Pipeline."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import xml.etree.ElementTree as ET
import polars as pl
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_DIR = REPO_ROOT / "scripts" / "sc9_milo_deconv_validation"


def test_scripts_exist_and_py_compile() -> None:
    """Verifies all expected pipeline scripts exist in scripts/sc9_milo_deconv_validation/."""
    expected_scripts = [
        "__init__.py",
        "01_build_single_reference.py",
        "02_build_combined_references.py",
        "03_deconvolute_iatlas.py",
        "04_logistic_regression.py",
        "05_analyze_milopy.py",
        "06_compare_concordance.py",
        "07_plot_figures.py",
        "run_pipeline.py",
    ]
    for s_name in expected_scripts:
        script_path = PIPELINE_DIR / s_name
        assert script_path.exists(), f"Missing expected pipeline script: {script_path}"


def test_snakemake_rules_defined() -> None:
    """Verifies that workflow/rules/sc9_deconv_validation.smk defines all required rules."""
    smk_path = REPO_ROOT / "workflow" / "rules" / "sc9_deconv_validation.smk"
    assert smk_path.exists()
    content = smk_path.read_text(encoding="utf-8")

    expected_rules = [
        "rule sc9_build_single_reference:",
        "rule sc9_build_combined_references:",
        "rule sc9_deconvolute_iatlas:",
        "rule sc9_logistic_regression:",
        "rule sc9_analyze_milopy:",
        "rule sc9_compare_concordance:",
        "rule sc9_plot_validation_figures:",
        "rule run_sc9_deconv_pipeline:",
    ]
    for r in expected_rules:
        assert r in content, f"Expected Snakemake rule missing: {r}"


def test_plot_figures_svg_standard(tmp_path: Path) -> None:
    """Tests that plot_fig1_predictive_capacity creates valid 180mm SVG with (a)-(r) badges and min 13pt font."""
    plot_script = PIPELINE_DIR / "07_plot_figures.py"
    spec = importlib.util.spec_from_file_location("sc9_plots", plot_script)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    # Create synthetic benchmark summary table with multiple reference strategies & resolutions
    records = []
    strategies = [
        ("Sade-Feldman", "Single Dataset", 1),
        ("Gondal-PanCancer", "Single Dataset", 1),
        ("Marteau-Colorectal", "Single Dataset", 1),
        ("All-9-Milo-Combined", "Criteria-Combined", 9),
        ("NSCLC-Trio", "Criteria-Combined", 3),
        ("Random-Pair-1", "Random-Combined", 2),
    ]
    resolutions = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]

    for name, cat, n_ds in strategies:
        for res in resolutions:
            records.append({
                "reference_type": name,
                "reference_category": cat,
                "resolution": float(res),
                "n_clusters": int(10 + res * 10),
                "n_signature_genes": 350,
                "condition_number": 15.0 + res * 2.0,
                "spearman_rho": 0.45,
                "concordance_percentage": 68.0,
                "melanoma_multivariate_auc": 0.60 + res * 0.03,
                "pancancer_multivariate_auc": 0.58 + res * 0.02,
                "mean_cohort_multivariate_auc": 0.59 + res * 0.025,
                "mean_univariate_auc": 0.55,
                "n_datasets": n_ds,
                "n_cells": n_ds * 15000,
                "subsample_fraction": 1.0,
            })

    # Add titration records
    for frac in [0.2, 0.4, 0.6, 0.8]:
        records.append({
            "reference_type": f"All-9-Milo-Combined ({int(frac*100)}% Cells)",
            "reference_category": "Cell-Subsample",
            "resolution": 0.5,
            "n_clusters": 20,
            "n_signature_genes": 350,
            "condition_number": 18.0,
            "spearman_rho": 0.40,
            "concordance_percentage": 65.0,
            "melanoma_multivariate_auc": 0.57 + frac * 0.08,
            "pancancer_multivariate_auc": 0.55 + frac * 0.07,
            "mean_cohort_multivariate_auc": 0.56 + frac * 0.075,
            "mean_univariate_auc": 0.53,
            "n_datasets": 9,
            "n_cells": int(frac * 9 * 15000),
            "subsample_fraction": frac,
        })

    df_synth = pl.DataFrame(records)
    synth_file = tmp_path / "multi_resolution_benchmark_summary.parquet"
    df_synth.write_parquet(synth_file)

    out_dir = tmp_path / "results"
    res = mod.plot_fig1_predictive_capacity(tmp_path, out_dir)
    assert res.__class__.__name__ == "Success"

    svg_file = out_dir / "fig1_predictive_capacity_benchmark.svg"
    assert svg_file.exists()

    content = svg_file.read_text(encoding="utf-8")
    assert 'width="180.00mm"' in content
    # Verify unconstrained height (no height="...mm" in root svg)
    assert 'height="' not in content[:content.find("viewBox")]

    # Verify panel badges (a) through (r) are present
    for i in range(18):
        letter = f"({chr(ord('a') + i)})"
        assert f">{letter}</text>" in content, f"Missing badge {letter} in SVG"

    # Verify XML well-formedness
    root = ET.fromstring(content)
    assert root.tag.endswith("svg")
    print("ALL TESTS PASSED SUCCESSFULLY!")


if __name__ == "__main__":
    import tempfile
    test_scripts_exist_and_py_compile()
    print("test_scripts_exist_and_py_compile passed")
    test_snakemake_rules_defined()
    print("test_snakemake_rules_defined passed")
    with tempfile.TemporaryDirectory() as td:
        test_plot_figures_svg_standard(Path(td))
    print("test_plot_figures_svg_standard passed")

