"""Native PyTorch Dataset with on-the-fly Negative Binomial augmentations."""

from __future__ import annotations

from typing import Callable, Mapping
import anndata as ad
import numpy as np
import scipy.sparse as sp
import torch
from returns.maybe import Maybe, Nothing, Some
from returns.result import Result
from torch.utils.data import Dataset


class TmeTorchDataset(Dataset):
    """PyTorch Dataset delivering single-cell or sample expression tensors with dynamic augmentations."""

    def __init__(
        self,
        adata: ad.AnnData,
        label_keys: tuple[str, ...] = ("response_binary", "cell_type"),
        on_the_fly_nb_dispersion: Maybe[float] = Nothing,
        transform: Maybe[Callable[[ad.AnnData], Result[ad.AnnData, str]]] = Nothing,
        device: str = "cpu",
    ) -> None:
        self.adata = adata
        self.label_keys = label_keys
        self.nb_dispersion = (
            on_the_fly_nb_dispersion.value_or(None)
            if isinstance(on_the_fly_nb_dispersion, Some)
            else None
        )
        self.transform = transform
        self.device = device

        # Precompute dense or CSR index
        self.X = adata.X
        self.is_sparse = sp.issparse(self.X)
        self.n_obs = adata.n_obs

        # Pre-extract labels
        self.labels: dict[str, np.ndarray] = {}
        for k in label_keys:
            if k in adata.obs.columns:
                col = adata.obs[k].to_numpy()
                if np.issubdtype(col.dtype, np.number):
                    self.labels[k] = np.nan_to_num(col.astype(np.float32), nan=-1.0)
                else:
                    unique_vals = {v: idx for idx, v in enumerate(np.unique(col))}
                    self.labels[k] = np.array([unique_vals[v] for v in col], dtype=np.int64)

    def __len__(self) -> int:
        return self.n_obs

    def __getitem__(self, idx: int) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
        # Extract row
        row = (
            np.asarray(self.X[idx].toarray()).flatten()
            if self.is_sparse
            else np.asarray(self.X[idx]).flatten()
        )

        # On-the-fly Negative Binomial augmentation if enabled
        if self.nb_dispersion is not None and self.nb_dispersion > 0.0:
            alpha = self.nb_dispersion
            shape = 1.0 / alpha
            mu = np.clip(row, 0.0, None)
            scale = alpha * mu
            lam = np.random.gamma(shape=shape, scale=np.maximum(scale, 1e-8))
            lam = np.where(mu > 0, lam, 0.0)
            row = np.random.poisson(lam).astype(np.float32)

        x_tensor = torch.from_numpy(row.astype(np.float32)).to(self.device)

        labels_tensor = {
            k: torch.tensor(v[idx]).to(self.device)
            for k, v in self.labels.items()
        }

        return x_tensor, labels_tensor
