"""Unit tests for the reference sampling HPO CLI script and diagnostic plotting."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import polars as pl
import pytest
from returns.result import Success

# Append script directory to import plotting
script_dir = Path(__file__).resolve().parents[3] / "scripts" / "reference_sampling_hpo"
sys.path.insert(0, str(script_dir))
from plotting import plot_held_out_validation, plot_pareto_frontier  # type: ignore[import-not-found]
from run_hpo_pipeline import (  # type: ignore[import-not-found]
    evaluate_held_out_generalization,
    generate_synthetic_benchmark_data,
    verify_execution_host,
)
from returns.result import Failure, Success
from tme_datasets.sampling.classifier_hpo import (  # type: ignore[import-untyped]
    ClassifierType,
    FeatureSelectorType,
    InnerClassifierConfig,
    TransformType,
)


def test_generate_synthetic_benchmark_data() -> None:
    """Verify synthetic benchmark data structure and class labels."""
    disc, held = generate_synthetic_benchmark_data()

    assert "Hugo-iAtlas" in disc
    assert "Riaz-iAtlas" in disc
    assert "Gide-iAtlas" in held
    assert "Rosenberg-iAtlas" in held

    X, genes, y = disc["Hugo-iAtlas"]
    assert X.shape[0] == 20
    assert len(genes) == 200
    assert len(np.unique(y)) == 2


def test_evaluate_held_out_generalization_synthetic() -> None:
    """Verify out-of-cohort generalization evaluation on synthetic data."""
    disc, held = generate_synthetic_benchmark_data()

    # Create synthetic reference matrix: 3 states x 200 genes
    rng = np.random.default_rng(42)
    genes = disc["Hugo-iAtlas"][1]
    n_states = 3
    ref_mat = rng.uniform(0.1, 1.0, size=(n_states, len(genes)))
    # Add strong signal in state 0 for marker genes
    ref_mat[0, :30] += 10.0
    state_labels = ("CD8_T", "Macrophage", "Tumor")

    cfg = InnerClassifierConfig(
        transform=TransformType.CLR,
        selector=FeatureSelectorType.NONE,
        classifier=ClassifierType.LOGISTIC_REGRESSION,
    )

    res = evaluate_held_out_generalization(
        best_ref_mat=ref_mat,
        best_ref_genes=genes,
        state_labels=state_labels,
        discovery_bulk=disc,
        held_out_bulk=held,
        best_inner_config=cfg,
    )

    assert isinstance(res, Success)
    held_aucs = res.unwrap()
    assert "Gide-iAtlas" in held_aucs
    assert "Rosenberg-iAtlas" in held_aucs
    assert all(0.0 <= score <= 1.0 for score in held_aucs.values())


def test_plot_pareto_frontier_svg(tmp_path: Path) -> None:
    """Verify that plot_pareto_frontier exports a valid vector SVG file."""
    df = pl.DataFrame({
        "trial_id": [1, 2, 3],
        "rung": ["rung_0_screening", "rung_1_refinement", "rung_2_full"],
        "mean_loco_auc": [0.55, 0.70, 0.82],
        "mean_loco_pr_auc": [0.50, 0.65, 0.75],
        "collinearity_max": [0.90, 0.82, 0.75],
        "condition_number": [500.0, 120.0, 35.0],
        "n_cell_states": [12, 10, 8],
        "n_shared_genes": [2000, 2200, 2400],
        "total_cells_sampled": [400, 1000, 2500],
        "elapsed_seconds": [2.5, 8.0, 22.0],
        "is_pruned": [False, False, False],
        "prune_reason": ["", "", ""],
    })

    out_file = tmp_path / "pareto_test.svg"
    res = plot_pareto_frontier(eval_df=df, pareto_trial_ids=[3], out_file=out_file)

    assert isinstance(res, Success)
    assert out_file.exists()
    assert out_file.stat().st_size > 500

    content = out_file.read_text(encoding="utf-8")
    assert "<svg" in content
    assert "</svg>" in content


def test_plot_held_out_validation_svg(tmp_path: Path) -> None:
    """Verify that plot_held_out_validation exports a valid vector SVG file."""
    held_aucs = {"Gide-iAtlas": 0.76, "Rosenberg-iAtlas": 0.68, "Choueiri-iAtlas": 0.72}
    out_file = tmp_path / "held_out_test.svg"

    res = plot_held_out_validation(held_out_aucs=held_aucs, discovery_mean_auc=0.74, out_file=out_file)

    assert isinstance(res, Success)
    assert out_file.exists()
    assert out_file.stat().st_size > 500

    content = out_file.read_text(encoding="utf-8")
    assert "<svg" in content
    assert "</svg>" in content


def test_verify_execution_host_blocks_local() -> None:
    """Verify that verify_execution_host blocks heavy runs on local workstation."""
    res = verify_execution_host(
        benchmark_mode=False,
        force_local=False,
        hostname="macbook-pro.local",
        allow_local_env="",
        storage_dir=Path("/nonexistent_storage_dir"),
    )
    assert isinstance(res, Failure)
    assert "Execution blocked" in res.failure()
    assert "macbook-pro.local" in res.failure()


def test_verify_execution_host_allows_benchmark() -> None:
    """Verify that benchmark mode is permitted on any host."""
    res = verify_execution_host(
        benchmark_mode=True,
        force_local=False,
        hostname="macbook-pro.local",
        allow_local_env="",
    )
    assert isinstance(res, Success)


def test_verify_execution_host_allows_force_local() -> None:
    """Verify that --force-local overrides host check."""
    res = verify_execution_host(
        benchmark_mode=False,
        force_local=True,
        hostname="macbook-pro.local",
        allow_local_env="",
    )
    assert isinstance(res, Success)


def test_verify_execution_host_allows_env_override() -> None:
    """Verify that ALLOW_LOCAL_HPO environment variable bypasses host restriction."""
    for val in ("1", "true", "True", "yes"):
        res = verify_execution_host(
            benchmark_mode=False,
            force_local=False,
            hostname="macbook-pro.local",
            allow_local_env=val,
        )
        assert isinstance(res, Success)


def test_verify_execution_host_allows_olm() -> None:
    """Verify that remote host 'olm' is accepted."""
    for host in ("olm", "olm.cluster.internal", "node-olm-01"):
        res = verify_execution_host(
            benchmark_mode=False,
            force_local=False,
            hostname=host,
            allow_local_env="",
        )
        assert isinstance(res, Success)

