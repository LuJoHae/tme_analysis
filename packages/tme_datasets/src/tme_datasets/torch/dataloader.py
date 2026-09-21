"""DataLoader construction and custom collators for TME PyTorch datasets."""

from __future__ import annotations

from torch.utils.data import DataLoader

from .dataset import TmeTorchDataset


def create_tme_dataloader(
    dataset: TmeTorchDataset,
    batch_size: int = 64,
    shuffle: bool = True,
    num_workers: int = 0,
    pin_memory: bool = False,
) -> DataLoader:
    """Construct a standard PyTorch DataLoader delivering augmented TME batches."""
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        pin_memory=pin_memory,
    )
