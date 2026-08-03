import pytest
import numpy as np
import torch
import polars as pl
from returns.result import Success, Failure
from returns.maybe import Some, Nothing
from bayesprism.models import GibbsControl, OptControl
from bayesprism.pipeline import new_prism, run_prism, get_fraction, get_exp


def test_end_to_end_pipeline() -> None:
    # 4 cells, 3 genes
    ref = np.array([
        [100, 20, 10],
        [90, 15, 5],
        [10, 80, 200],
        [5, 95, 180]
    ], dtype=np.float32)

    cell_type_labels = ("T_cell", "T_cell", "B_cell", "B_cell")
    cell_state_labels = ("T_state1", "T_state2", "B_state1", "B_state2")

    mixture = np.array([
        [500, 400, 600],
        [200, 800, 1000]
    ], dtype=np.float32)

    gene_names = ("geneA", "geneB", "geneC")

    # Step 1: new_prism
    prism_res = new_prism(
        reference=ref,
        cell_type_labels=cell_type_labels,
        cell_state_labels=Some(cell_state_labels),
        mixture=mixture,
        gene_names=Some(gene_names),
        pseudo_min=1e-8,
    )

    match prism_res:
        case Failure(err):
            pytest.fail(f"new_prism failed: {err}")
        case Success(prism_obj):
            pass

    # Step 2: run_prism
    g_ctrl = GibbsControl(chain_length=50, burn_in=10, thinning=2, seed=123)
    o_ctrl = OptControl(maxit=50, optimizer="MAP")

    run_res = run_prism(
        prism=prism_obj,
        update_gibbs=True,
        gibbs_control=Some(g_ctrl),
        opt_control=Some(o_ctrl),
    )

    match run_res:
        case Failure(err):
            pytest.fail(f"run_prism failed: {err}")
        case Success(bp_obj):
            pass

    # Step 3: get_fraction
    frac_res = get_fraction(bp_obj, which_theta="final")
    match frac_res:
        case Failure(err):
            pytest.fail(f"get_fraction failed: {err}")
        case Success(df_frac):
            assert isinstance(df_frac, pl.DataFrame)
            assert df_frac.shape == (2, 3)  # bulk_id + 2 cell types

    # Step 4: get_exp
    exp_res = get_exp(bp_obj, state_or_type="type", cell_name="T_cell")
    match exp_res:
        case Failure(err):
            pytest.fail(f"get_exp failed: {err}")
        case Success(df_exp):
            assert isinstance(df_exp, pl.DataFrame)
            assert df_exp.shape == (2, 4)  # bulk_id + 3 genes
