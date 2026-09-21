"""Unit tests for harmonization, metrics, pseudobulk simulation, and PyTorch bridge."""

import anndata as ad
import numpy as np
import pandas as pd
from returns.result import Success
from tme_datasets.harmonization.align import align_and_concatenate
from tme_datasets.harmonization.metrics import evaluate_integration_metrics
from tme_datasets.models import HarmonizeConfig, PseudobulkConfig
from tme_datasets.simulation.simulator import simulate_pseudobulk
from tme_datasets.torch.dataloader import create_tme_dataloader
from tme_datasets.torch.dataset import TmeTorchDataset
from tme_datasets.types import HarmonizeMode


def test_align_and_concatenate_intersection(mock_single_cell_adata: ad.AnnData) -> None:
    adata1 = mock_single_cell_adata[:50].copy()
    adata2 = mock_single_cell_adata[50:].copy()

    config = HarmonizeConfig(mode=HarmonizeMode.INTERSECTION, min_shared_genes=5)
    res = align_and_concatenate([adata1, adata2], ["batch1", "batch2"], config)
    assert isinstance(res, Success)
    comb = res.unwrap()
    assert comb.n_obs == 100
    assert "dataset_id" in comb.obs.columns


def test_align_and_concatenate_union_zero_filled() -> None:
    # Adata 1 with genes A, B
    a1 = ad.AnnData(X=np.array([[1.0, 2.0]], dtype=np.float32), var=pd.DataFrame(index=["GENE_A", "GENE_B"]))
    # Adata 2 with genes B, C
    a2 = ad.AnnData(X=np.array([[3.0, 4.0]], dtype=np.float32), var=pd.DataFrame(index=["GENE_B", "GENE_C"]))

    config = HarmonizeConfig(mode=HarmonizeMode.UNION_ZERO_FILLED)
    res = align_and_concatenate([a1, a2], ["d1", "d2"], config)
    assert isinstance(res, Success)
    comb = res.unwrap()
    assert comb.n_obs == 2
    assert comb.n_vars == 3
    assert set(comb.var_names) == {"GENE_A", "GENE_B", "GENE_C"}


def test_evaluate_integration_metrics(mock_single_cell_adata: ad.AnnData) -> None:
    res = evaluate_integration_metrics(mock_single_cell_adata, batch_key="batch", label_key="cell_type")
    assert isinstance(res, Success)
    metrics = res.unwrap()
    assert metrics.mean_ilisi >= 1.0
    assert metrics.mean_clisi >= 1.0
    assert 0.0 <= metrics.kbet_acceptance_rate <= 1.0


def test_simulate_pseudobulk(mock_single_cell_adata: ad.AnnData) -> None:
    config = PseudobulkConfig(n_samples=10, cells_per_sample=100)
    res = simulate_pseudobulk(mock_single_cell_adata, config, cell_type_key="cell_type")
    assert isinstance(res, Success)
    bulk_adata, truth_df = res.unwrap()
    assert bulk_adata.n_obs == 10
    assert truth_df.height == 10
    assert "sample_id" in truth_df.columns
    assert "CD8_T" in truth_df.columns


def test_torch_dataset_and_loader(mock_single_cell_adata: ad.AnnData) -> None:
    torch_ds = TmeTorchDataset(mock_single_cell_adata, label_keys=("response_binary",))
    assert len(torch_ds) == mock_single_cell_adata.n_obs

    x_tensor, label_dict = torch_ds[0]
    assert x_tensor.shape[0] == mock_single_cell_adata.n_vars
    assert "response_binary" in label_dict

    loader = create_tme_dataloader(torch_ds, batch_size=16, shuffle=False)
    for batch_x, batch_y in loader:
        assert batch_x.shape == (16, mock_single_cell_adata.n_vars)
        assert "response_binary" in batch_y
        break
