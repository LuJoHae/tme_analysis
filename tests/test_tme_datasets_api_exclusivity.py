"""Unit and AST invariant tests enforcing 100% tme_datasets API loading exclusivity.

Verifies that:
1. Ingestion, batch preprocessing, and single-cell QC pipelines live inside
   ``packages/tme_datasets/src/tme_datasets/pipeline/``, leaving 37 analytical and
   plotting scripts (``2a1`` through ``4e``) in ``scripts/manuscript_pipeline/``.
2. Zero scripts in ``scripts/manuscript_pipeline/`` or downstream consumer suites
   call ``read_h5ad``, ``read_zarr``, ``find_dataset_h5ad``, or ``datalair.Lair``
   directly, or glob ``preprocessed/*.h5ad`` outside the ``tme_datasets`` API.
3. ``tme_datasets.load_backed`` resolves dataset IDs via the ``tme_datasets``
   registry/path resolver without requiring manual ``.h5ad`` path construction.
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Final

import anndata as ad  # type: ignore[import-untyped]
from hypothesis import given, settings, strategies as st
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict
from returns.pipeline import is_successful

import tme_datasets
from tme_datasets import (
    load_backed,
    run_batch_preprocess_cohorts,
    run_parallel_download_and_preprocess,
    run_single_cell_qc_pipeline,
)

REPO_ROOT: Final[Path] = Path(__file__).resolve().parent.parent
MANUSCRIPT_PIPELINE_DIR: Final[Path] = REPO_ROOT / "scripts" / "manuscript_pipeline"

FORBIDDEN_LOADER_CALLS: Final[frozenset[str]] = frozenset(
    {
        "read_h5ad",
        "read_zarr",
        "find_dataset_h5ad",
        "Lair",
    }
)

ALLOWED_TME_API_LOADERS: Final[frozenset[str]] = frozenset(
    {
        "load_dataset",
        "load_empirical_pseudobulk",
        "aggregate_empirical_pseudobulk",
        "load_iatlas_cohort_or_combined",
        "load_backed",
        "resolve_milopy_cohorts",
        "list_preprocessed_datasets",
        "sample_single_cell_cohorts",
        "IMMUNOTHERAPY_GENE_PANEL",
        "load_coulton_tam_atlas",
        "load_gao_caf_atlas",
        "load_cords_caf_atlas",
        "load_buechler_caf_atlas",
        "load_zheng_tcell_atlas",
        "load_chu_tcell_atlas",
        "load_shi_pancancer_atlas",
    }
)

DOWNSTREAM_CONSUMER_SCRIPTS: Final[tuple[str, ...]] = (
    "scripts/stability_selection_iatlas/shared.py",
    "scripts/reference_sampling_hpo/run_all_cohorts_hpo_and_benchmarking.py",
    "scripts/reference_sampling_hpo/run_single_cohort_hpo.py",
    "scripts/concordance/01_sc_patient_response_association.py",
    "scripts/concordance/02_build_curated_reference.py",
    "scripts/concordance/03_insilico_pseudobulk_validation.py",
    "scripts/concordance_calibration/01_extract_sc_ground_truth.py",
    "scripts/ingestion/04_check_milopy_compatibility.py",
    "scripts/ingestion/04_plot_dataset_embeddings.py",
    "scripts/ingestion/05_run_milopy_da.py",
    "scripts/ingestion/06_investigate_subtle_datasets.py",
    "scripts/ingestion/07_integrate_combined_cohorts.py",
    "scripts/ingestion/plot_unified_embeddings.py",
    "scripts/ingestion/parallel_download_and_preprocess.py",
    "scripts/gse123139/plot_gse123139.py",
    "scripts/gse97168/plot_gse97168.py",
    "scripts/coulton_tam_comparison/compare_meta42_coulton_tam.py",
    "scripts/pan_cancer_atlas_comparison/compare_meta42_external_atlases.py",
)


class ScriptAstAudit(BaseModel):
    """Immutable AST audit result for a single Python script."""

    model_config = ConfigDict(frozen=True)

    rel_path: str
    forbidden_calls: tuple[str, ...]
    h5ad_globs: tuple[str, ...]
    tme_api_imports: tuple[str, ...]


def extract_call_name(node: ast.Call) -> str:
    """Extract function or attribute name from an AST Call node."""
    match node.func:
        case ast.Name(id=fn_id):
            return fn_id
        case ast.Attribute(attr=attr_name):
            return attr_name
        case _:
            return ""


def audit_python_script(script_path: Path) -> ScriptAstAudit:
    """Parse a Python file into an AST and collect loader calls and imports purely."""
    source = script_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(script_path))

    forbidden: list[str] = []
    globs: list[str] = []
    tme_imports: list[str] = []

    for node in ast.walk(tree):
        match node:
            case ast.Call() as call_node:
                fn_name = extract_call_name(call_node)
                if fn_name in FORBIDDEN_LOADER_CALLS:
                    forbidden.append(f"{fn_name}@L{call_node.lineno}")
                if fn_name in {"glob", "rglob"} and call_node.args:
                    first_arg = call_node.args[0]
                    if (
                        isinstance(first_arg, ast.Constant)
                        and isinstance(first_arg.value, str)
                        and ".h5ad" in first_arg.value
                    ):
                        globs.append(f"{fn_name}({first_arg.value!r})@L{call_node.lineno}")
            case ast.ImportFrom(module=mod_name, names=aliases) if (
                mod_name is not None and mod_name.startswith("tme_datasets")
            ):
                for alias in aliases:
                    if alias.name in ALLOWED_TME_API_LOADERS:
                        tme_imports.append(alias.name)
            case _:
                pass

    return ScriptAstAudit(
        rel_path=str(script_path.relative_to(REPO_ROOT)),
        forbidden_calls=tuple(forbidden),
        h5ad_globs=tuple(globs),
        tme_api_imports=tuple(sorted(set(tme_imports))),
    )


def test_pipeline_migration_and_37_manuscript_scripts() -> None:
    """Verify 1a-1c migrated to tme_datasets.pipeline and 37 scripts remain in manuscript_pipeline."""
    assert callable(run_parallel_download_and_preprocess)
    assert callable(run_batch_preprocess_cohorts)
    assert callable(run_single_cell_qc_pipeline)
    assert hasattr(tme_datasets, "pipeline")

    mp_scripts = sorted(
        p.name for p in MANUSCRIPT_PIPELINE_DIR.glob("*.py") if p.name != "__init__.py"
    )
    assert len(mp_scripts) >= 37, f"Expected >= 37 scripts in manuscript_pipeline, got {len(mp_scripts)}: {mp_scripts}"

    migrated = {
        "1a_download_and_preprocess_all.py",
        "1b_batch_preprocess_cohorts.py",
        "1c_run_qc_pipeline_all_single_cell.py",
        "3a_sade_feldman_pseudobulk_deconv_benchmark.py",
    }
    assert migrated.isdisjoint(set(mp_scripts))
    assert "3a_cross_cohort_pseudobulk_deconv_benchmark.py" in mp_scripts


def test_manuscript_pipeline_100_percent_tme_datasets_exclusivity() -> None:
    """Ensure zero manuscript_pipeline scripts bypass tme_datasets to read preprocessed .h5ad/.zarr/Lair."""
    mp_paths = sorted(p for p in MANUSCRIPT_PIPELINE_DIR.glob("*.py") if p.name != "__init__.py")
    audits = tuple(audit_python_script(p) for p in mp_paths)

    violations = tuple(
        a for a in audits if len(a.forbidden_calls) > 0 or len(a.h5ad_globs) > 0
    )
    assert len(violations) == 0, f"Forbidden non-API dataset loads detected: {violations}"

    # Verify primary cohort-loading scripts explicitly import tme_datasets loaders
    audit_by_name = {Path(a.rel_path).name: a for a in audits}
    expected_api_loaders: dict[str, str] = {
        "2a1_sample_and_analyze_cohorts.py": "resolve_milopy_cohorts",
        "2a2_combine_cohorts_and_analyze.py": "load_dataset",
        "2b_run_stability_selection_iatlas.py": "IMMUNOTHERAPY_GENE_PANEL",
        "3a_cross_cohort_pseudobulk_deconv_benchmark.py": "load_empirical_pseudobulk",
        "3b_pan_cohort_unified_deconv_benchmark.py": "load_empirical_pseudobulk",
        "3c_iatlas_instaprism_milopy_benchmark.py": "load_dataset",
        "3d_curate_expert_gene_list.py": "load_iatlas_cohort_or_combined",
        "3e1_benchmark_iatlas_icb_predictors.py": "load_dataset",
        "3e2_run_hpo_pipeline.py": "sample_single_cell_cohorts",
        "3f1_build_sf_reference.py": "load_dataset",
        "3f2_build_integrated_reference.py": "load_dataset",
        "3f3_build_dataset_reference.py": "load_dataset",
        "3f5_deconvolute_iatlas_multires.py": "load_dataset",
        "3f6_logistic_regression_multires.py": "load_dataset",
        "3f7_analyze_milopy_multires.py": "load_dataset",
        "4a1_visualize_cohort_umaps.py": "load_backed",
    }
    for script_name, expected_fn in expected_api_loaders.items():
        script_audit = audit_by_name[script_name]
        assert expected_fn in script_audit.tme_api_imports, (
            f"{script_name} must import {expected_fn} from tme_datasets; "
            f"found {script_audit.tme_api_imports}"
        )

    assert "load_dataset" in audit_by_name["3a_cross_cohort_pseudobulk_deconv_benchmark.py"].tme_api_imports
    stage3b_src = (MANUSCRIPT_PIPELINE_DIR / "3b_pan_cohort_unified_deconv_benchmark.py").read_text(
        encoding="utf-8"
    )
    assert "pseudobulk.parquet" not in stage3b_src, (
        "3b_pan_cohort_unified_deconv_benchmark.py must load pseudobulk strictly via "
        "tme_datasets.load_empirical_pseudobulk() rather than reading pseudobulk.parquet directly."
    )

    shared_stab_audit = audit_python_script(REPO_ROOT / "scripts/stability_selection_iatlas/shared.py")
    assert "load_iatlas_cohort_or_combined" in shared_stab_audit.tme_api_imports


def test_downstream_consumer_scripts_tme_datasets_exclusivity() -> None:
    """Verify refactored auxiliary consumer scripts also load preprocessed cohorts strictly via tme_datasets."""
    for rel_str in DOWNSTREAM_CONSUMER_SCRIPTS:
        script_path = REPO_ROOT / rel_str
        assert script_path.exists(), f"Missing consumer script: {rel_str}"
        audit = audit_python_script(script_path)
        assert len(audit.forbidden_calls) == 0, (
            f"{rel_str} has forbidden non-API dataset calls: {audit.forbidden_calls}"
        )
        assert len(audit.h5ad_globs) == 0, (
            f"{rel_str} globs .h5ad directly instead of calling list_preprocessed_datasets: {audit.h5ad_globs}"
        )


def test_load_backed_resolves_dataset_id_via_base_dir(tmp_path: Path) -> None:
    """Verify load_backed resolves a registered dataset ID from preprocessed/ without manual path building."""
    preproc_dir = tmp_path / "preprocessed"
    preproc_dir.mkdir(parents=True, exist_ok=True)
    obs = pd.DataFrame({"response": ["R", "NR"]}, index=["c1", "c2"])
    var = pd.DataFrame(index=["PDCD1", "CD8A"])
    adata = ad.AnnData(X=np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32), obs=obs, var=var)
    adata.write_h5ad(preproc_dir / "GSE120575.h5ad")

    res = load_backed("GSE120575", mode="r", base_dir=tmp_path)
    assert is_successful(res), f"Expected load_backed to resolve GSE120575, got {res}"
    backed_adata = res.unwrap()
    assert backed_adata.n_obs == 2
    assert backed_adata.n_vars == 2
    if backed_adata.file.is_open:
        backed_adata.file.close()


@given(
    fn_name=st.sampled_from(sorted(FORBIDDEN_LOADER_CALLS)),
    is_attr=st.booleans(),
)
@settings(max_examples=16)
def test_extract_call_name_detects_all_forbidden_patterns(fn_name: str, is_attr: bool) -> None:
    """Property test: AST call extractor detects both bare and module-qualified forbidden loader calls."""
    expr_src = f"mod.{fn_name}('x.h5ad')" if is_attr else f"{fn_name}('x.h5ad')"
    tree = ast.parse(expr_src, mode="eval")
    assert isinstance(tree.body, ast.Call)
    assert extract_call_name(tree.body) == fn_name
