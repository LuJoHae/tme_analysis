import datalair
from gene_utils import download_from_cbioportal


class CBioPortalDataset(datalair.Dataset):
    _registry = {
        "VanAllen": "skcm_dfci_2015.tar.gz",
        "Snyder": "skcm_mskcc_2014.tar.gz",
        "Liu": "mel_dfci_2019.tar.gz",
        "Hugo": "mel_ucla_2016.tar.gz",
        "Catalanotti": "skcm_vanderbilt_mskcc_2015.tar.gz",
        "Cloughesy-iAtlas": "gbm_iatlas_prins_2019.tar.gz",
        "Rosenberg-iAtlas": "blca_iatlas_imvigor210_2017.tar.gz",
        "Riaz-iAtlas": "mel_iatlas_riaz_nivolumab_2017.tar.gz",
        "Liu-iAtlas": "mel_iatlas_liu_2019.tar.gz",
        "Gide-iAtlas": "mel_iatlas_gide_2019.tar.gz",
        "Hugo-iAtlas": "mel_iatlas_hugo_ucla_2016.tar.gz",
        "Padron-iAtlas": "paad_iatlas_prince_2022.tar.gz",
        "Choueiri-iAtlas": "ccrcc_iatlas_choueiri_2016.tar.gz",
        "Anders-iAtlas": "brca_iatlas_anders_2022.tar.gz",
        "McDermott-iAtlas": "rcc_iatlas_immotion150_2018.tar.gz"
    }

    def __init__(self, name: str) -> None:
        if name not in self._registry:
            raise ValueError(f"Unmapped dataset identifier: {name}")
        super().__init__(namespace="CBioPortalDataset", dataset_name=name)

    def derive(self, lair: datalair.Lair) -> None:
        output_dir = lair.get_path(self)
        download_from_cbioportal(self._registry[self._dataset_name], output_dir)

    @classmethod
    def get_dataset_name(cls):
        return list(cls._registry.keys())
