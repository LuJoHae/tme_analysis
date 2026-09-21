"""Compatibility shim mirroring legacy ici_datasets interface."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
from returns.result import Success

from ..query import load_dataset


class CBioPortalDataset:
    """Legacy CBioPortalDataset shim delegating to tme_datasets."""

    _registry = {
        "Hugo-iAtlas": "mel_iatlas_hugo_ucla_2016.tar.gz",
        "Riaz-iAtlas": "mel_iatlas_riaz_nivolumab_2017.tar.gz",
        "Liu-iAtlas": "mel_iatlas_liu_2019.tar.gz",
        "Gide-iAtlas": "mel_iatlas_gide_2019.tar.gz",
        "Rosenberg-iAtlas": "blca_iatlas_imvigor210_2017.tar.gz",
        "Padron-iAtlas": "paad_iatlas_prince_2022.tar.gz",
        "Anders-iAtlas": "brca_iatlas_anders_2022.tar.gz",
        "McDermott-iAtlas": "rcc_iatlas_immotion150_2018.tar.gz",
        "Choueiri-iAtlas": "ccrcc_iatlas_choueiri_2016.tar.gz",
    }

    def __init__(self, name: str) -> None:
        self.name = name

    def load(self, base_dir: Path | None = None) -> ad.AnnData:
        res = load_dataset(self.name, base_dir=base_dir)
        if isinstance(res, Success):
            return res.unwrap()
        raise RuntimeError(f"Failed to load via compat: {res}")
