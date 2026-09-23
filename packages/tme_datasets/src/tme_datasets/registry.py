"""Declarative registry of all single-cell, spatial, bulk, and signature datasets."""

from __future__ import annotations

from typing import Mapping
from returns.maybe import Maybe, Nothing, Some

from .models import DatasetSpec, RegisteredDatasets
from .types import Modality

DATASET_REGISTRY: Mapping[str, DatasetSpec] = {
    # -------------------------------------------------------------------------
    # Single-Cell Reference Pool (docs/dataset_catalog/01_single_cell_reference_datasets.md)
    # -------------------------------------------------------------------------
    "GSE120575": DatasetSpec(
        id="GSE120575",
        title="Sade-Feldman et al. 2018 Melanoma T-cell ICB Response",
        modality=Modality.SINGLE_CELL,
        cancer_type="Melanoma",
        platform="Smart-seq2",
        organ=Some("Skin"),
        has_response_labels=True,
        n_samples_or_cells=Some(16288),
    ),
    "GSE115978": DatasetSpec(
        id="GSE115978",
        title="Jerby-Arnon et al. 2018 Melanoma T-cell Exclusion",
        modality=Modality.SINGLE_CELL,
        cancer_type="Melanoma",
        platform="Smart-seq2",
        organ=Some("Skin"),
        has_response_labels=False,
        n_samples_or_cells=Some(7186),
    ),
    "GSE125449": DatasetSpec(
        id="GSE125449",
        title="Ma et al. 2019 Liver Cancer (HCC / ICC) TME",
        modality=Modality.SINGLE_CELL,
        cancer_type="Liver Cancer",
        platform="10x_v2",
        organ=Some("Liver"),
        has_response_labels=False,
        n_samples_or_cells=Some(5115),
    ),
    "GSE123813": DatasetSpec(
        id="GSE123813",
        title="Yost et al. 2019 BCC Clonal Replacement",
        modality=Modality.SINGLE_CELL,
        cancer_type="Basal Cell Carcinoma",
        platform="10x_5prime",
        organ=Some("Skin"),
        has_response_labels=False,
        n_samples_or_cells=Some(3500),
    ),
    "Maynard_NSCLC": DatasetSpec(
        id="Maynard_NSCLC",
        title="Maynard et al. 2020 Therapy-Induced NSCLC Evolution",
        modality=Modality.SINGLE_CELL,
        cancer_type="Non-Small Cell Lung",
        platform="10x_v3",
        organ=Some("Lung"),
        has_response_labels=False,
        n_samples_or_cells=Some(3000),
    ),
    "GSE179994": DatasetSpec(
        id="GSE179994",
        title="Tietscher et al. Pan-Cancer T-cell Infiltrate Atlas",
        modality=Modality.SINGLE_CELL,
        cancer_type="Pan-Cancer",
        platform="10x_v3",
        organ=Some("Multi-organ"),
        has_response_labels=False,
        n_samples_or_cells=Some(150849),
    ),

    # -------------------------------------------------------------------------
    # 17 Pan-Cancer Single-Cell Reference Atlas Cohorts
    # -------------------------------------------------------------------------
    "GSE178341": DatasetSpec(
        id="GSE178341",
        title="Pelka et al. 2021 Colorectal Cancer Multicellular Hubs",
        modality=Modality.SINGLE_CELL,
        cancer_type="Colorectal Cancer",
        platform="10x_v3",
        organ=Some("Colon"),
        has_response_labels=False,
        n_samples_or_cells=Some(65000),
    ),
    "GSE114727": DatasetSpec(
        id="GSE114727",
        title="Azizi et al. 2018 Breast Cancer Immune Phenotypes",
        modality=Modality.SINGLE_CELL,
        cancer_type="Breast Cancer",
        platform="10x_v2",
        organ=Some("Breast"),
        has_response_labels=False,
        n_samples_or_cells=Some(45000),
    ),
    "E-MTAB-8107": DatasetSpec(
        id="E-MTAB-8107",
        title="Qian et al. 2020 Pan-Cancer Blueprint (Lung/CRC/OV/BRCA)",
        modality=Modality.SINGLE_CELL,
        cancer_type="Pan-Cancer",
        platform="10x_v2",
        organ=Some("Multi-organ"),
        has_response_labels=False,
        n_samples_or_cells=Some(200000),
    ),
    "GSE154763": DatasetSpec(
        id="GSE154763",
        title="Cheng et al. 2021 Pan-Cancer T-cell Atlas",
        modality=Modality.SINGLE_CELL,
        cancer_type="Pan-Cancer T-Cells",
        platform="10x_v3",
        organ=Some("Multi-organ"),
        has_response_labels=False,
        n_samples_or_cells=Some(390000),
    ),
    "GSE154826": DatasetSpec(
        id="GSE154826",
        title="Leader et al. 2021 NSCLC Microenvironment",
        modality=Modality.SINGLE_CELL,
        cancer_type="Non-Small Cell Lung",
        platform="10x_v3",
        organ=Some("Lung"),
        has_response_labels=False,
        n_samples_or_cells=Some(35000),
    ),
    "GSE131907": DatasetSpec(
        id="GSE131907",
        title="Kim et al. 2020 Lung Adenocarcinoma Heterogeneity",
        modality=Modality.SINGLE_CELL,
        cancer_type="Lung Adenocarcinoma",
        platform="10x_v2",
        organ=Some("Lung"),
        has_response_labels=False,
        n_samples_or_cells=Some(40000),
    ),
    "GSE201349": DatasetSpec(
        id="GSE201349",
        title="Becker et al. 2022 Colorectal Cancer Continuum",
        modality=Modality.SINGLE_CELL,
        cancer_type="Colorectal Cancer",
        platform="10x_v3",
        organ=Some("Colon"),
        has_response_labels=False,
        n_samples_or_cells=Some(30000),
    ),
    "GSE200997": DatasetSpec(
        id="GSE200997",
        title="Khaliq et al. 2022 Colorectal Cancer Classification",
        modality=Modality.SINGLE_CELL,
        cancer_type="Colorectal Cancer",
        platform="10x_v3",
        organ=Some("Colon"),
        has_response_labels=False,
        n_samples_or_cells=Some(25000),
    ),
    "GSE121638": DatasetSpec(
        id="GSE121638",
        title="Borcherding et al. 2021 ccRCC Immune Environment",
        modality=Modality.SINGLE_CELL,
        cancer_type="Clear Cell Renal Cell",
        platform="10x_v3",
        organ=Some("Kidney"),
        has_response_labels=False,
        n_samples_or_cells=Some(25000),
    ),
    "GSE156625": DatasetSpec(
        id="GSE156625",
        title="Sharma et al. 2020 HCC Onco-fetal Reprogramming",
        modality=Modality.SINGLE_CELL,
        cancer_type="Hepatocellular Carcinoma",
        platform="10x_v2",
        organ=Some("Liver"),
        has_response_labels=False,
        n_samples_or_cells=Some(15000),
    ),
    "GSE149614": DatasetSpec(
        id="GSE149614",
        title="Lu et al. 2022 HCC Multicellular Ecosystem",
        modality=Modality.SINGLE_CELL,
        cancer_type="Hepatocellular Carcinoma",
        platform="10x_v3",
        organ=Some("Liver"),
        has_response_labels=False,
        n_samples_or_cells=Some(18000),
    ),
    "GSE184362": DatasetSpec(
        id="GSE184362",
        title="Pu et al. 2021 Papillary Thyroid Carcinoma",
        modality=Modality.SINGLE_CELL,
        cancer_type="Papillary Thyroid",
        platform="10x_v3",
        organ=Some("Thyroid"),
        has_response_labels=False,
        n_samples_or_cells=Some(20000),
    ),
    "GSE139829": DatasetSpec(
        id="GSE139829",
        title="Durante et al. 2020 Uveal Melanoma Dynamics",
        modality=Modality.SINGLE_CELL,
        cancer_type="Uveal Melanoma",
        platform="10x_v2",
        organ=Some("Eye"),
        has_response_labels=False,
        n_samples_or_cells=Some(10000),
    ),
    "GSE200218": DatasetSpec(
        id="GSE200218",
        title="Biermann et al. 2022 Melanoma Brain Metastasis",
        modality=Modality.SINGLE_CELL,
        cancer_type="Melanoma Brain Met",
        platform="10x_v3",
        organ=Some("Brain"),
        has_response_labels=False,
        n_samples_or_cells=Some(12000),
    ),
    "GSE180661": DatasetSpec(
        id="GSE180661",
        title="Vazquez et al. 2022 Ovarian Cancer Ecosystem",
        modality=Modality.SINGLE_CELL,
        cancer_type="Ovarian Cancer",
        platform="10x_v3",
        organ=Some("Ovary"),
        has_response_labels=False,
        n_samples_or_cells=Some(15000),
    ),
    "GSE169246": DatasetSpec(
        id="GSE169246",
        title="Zhang et al. 2021 Triple-Negative Breast Cancer",
        modality=Modality.SINGLE_CELL,
        cancer_type="Triple-Negative Breast",
        platform="10x_v3",
        organ=Some("Breast"),
        has_response_labels=False,
        n_samples_or_cells=Some(18000),
    ),
    "GSE215120": DatasetSpec(
        id="GSE215120",
        title="Zhang et al. 2022 Pan-Cancer Myeloid Dynamics",
        modality=Modality.SINGLE_CELL,
        cancer_type="Pan-Cancer Myeloid",
        platform="10x_v3",
        organ=Some("Multi-organ"),
        has_response_labels=False,
        n_samples_or_cells=Some(50000),
    ),

    # -------------------------------------------------------------------------
    # Bulk iAtlas Validation Cohorts (docs/dataset_catalog/02_bulk_validation_cohorts.md)
    # -------------------------------------------------------------------------
    "Hugo-iAtlas": DatasetSpec(
        id="Hugo-iAtlas",
        title="Hugo et al. 2016 Melanoma anti-PD-1",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
        n_samples_or_cells=Some(27),
    ),
    "Riaz-iAtlas": DatasetSpec(
        id="Riaz-iAtlas",
        title="Riaz et al. 2017 Melanoma Nivolumab",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
        n_samples_or_cells=Some(107),
    ),
    "Liu-iAtlas": DatasetSpec(
        id="Liu-iAtlas",
        title="Liu et al. 2019 Melanoma anti-PD-1",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
        n_samples_or_cells=Some(122),
    ),
    "Gide-iAtlas": DatasetSpec(
        id="Gide-iAtlas",
        title="Gide et al. 2019 Melanoma anti-PD-1 +/- anti-CTLA-4",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
        n_samples_or_cells=Some(91),
    ),
    "Rosenberg-iAtlas": DatasetSpec(
        id="Rosenberg-iAtlas",
        title="Rosenberg et al. 2016 / IMvigor210 Bladder Atezolizumab",
        modality=Modality.BULK_RNA,
        cancer_type="Urothelial Bladder",
        platform="RNA-seq",
        organ=Some("Bladder"),
        has_response_labels=True,
        n_samples_or_cells=Some(347),
    ),
    "Padron-iAtlas": DatasetSpec(
        id="Padron-iAtlas",
        title="Padron et al. 2022 Pancreatic anti-PD-1 + CD40",
        modality=Modality.BULK_RNA,
        cancer_type="Pancreatic Ductal",
        platform="RNA-seq",
        organ=Some("Pancreas"),
        has_response_labels=True,
        n_samples_or_cells=Some(93),
    ),
    "Anders-iAtlas": DatasetSpec(
        id="Anders-iAtlas",
        title="Anders et al. 2021 Triple-Negative Breast Atezolizumab",
        modality=Modality.BULK_RNA,
        cancer_type="Triple-Negative Breast",
        platform="RNA-seq",
        organ=Some("Breast"),
        has_response_labels=True,
        n_samples_or_cells=Some(31),
    ),
    "McDermott-iAtlas": DatasetSpec(
        id="McDermott-iAtlas",
        title="McDermott et al. 2018 Renal Cell Atezolizumab +/- Bevacizumab",
        modality=Modality.BULK_RNA,
        cancer_type="Clear Cell Renal Cell",
        platform="RNA-seq",
        organ=Some("Kidney"),
        has_response_labels=True,
        n_samples_or_cells=Some(263),
    ),
    "Choueiri-iAtlas": DatasetSpec(
        id="Choueiri-iAtlas",
        title="Choueiri et al. 2016 Renal Cell Nivolumab",
        modality=Modality.BULK_RNA,
        cancer_type="Clear Cell Renal Cell",
        platform="RNA-seq",
        organ=Some("Kidney"),
        has_response_labels=True,
        n_samples_or_cells=Some(16),
    ),

    # -------------------------------------------------------------------------
    # Direct Paper ICI Response Cohorts (dataset_papers/README.md)
    # -------------------------------------------------------------------------
    "Auslander": DatasetSpec(
        id="Auslander",
        title="Auslander et al. 2018 Melanoma Response",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
    ),
    "Chen-CTLA4": DatasetSpec(
        id="Chen-CTLA4",
        title="Chen et al. 2016 Melanoma anti-CTLA-4",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
    ),
    "Chen-PD1": DatasetSpec(
        id="Chen-PD1",
        title="Chen et al. 2016 Melanoma anti-PD-1",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
    ),
    "Freeman": DatasetSpec(
        id="Freeman",
        title="Freeman et al. 2022 Melanoma ICB Response",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
    ),
    "VanAllen": DatasetSpec(
        id="VanAllen",
        title="Van Allen et al. 2015 Melanoma CTLA-4",
        modality=Modality.BULK_RNA,
        cancer_type="Melanoma",
        platform="RNA-seq",
        organ=Some("Skin"),
        has_response_labels=True,
    ),
    "EGAD00001006631": DatasetSpec(
        id="EGAD00001006631",
        title="Genentech IMvigor210 Align ReadsPerGene Counts",
        modality=Modality.BULK_RNA,
        cancer_type="Urothelial Bladder",
        platform="RNA-seq",
        organ=Some("Bladder"),
        has_response_labels=True,
    ),
}

# -------------------------------------------------------------------------
# Canonical Alias Mapping for Single-Cell and Atlas Datasets
# -------------------------------------------------------------------------
DATASET_ALIASES: Mapping[str, str] = {
    # Existing Single-Cell Benchmarks
    "Sade-Feldman": "GSE120575",
    "SadeFeldman": "GSE120575",
    "SadeFeldmanDefiningTCell2018Adata": "GSE120575",
    "Jerby-Arnon": "GSE115978",
    "JerbyArnon": "GSE115978",
    "JerbyArnonCancerCellProgram2018Adata": "GSE115978",
    "Ma": "GSE125449",
    "Ma_Liver": "GSE125449",
    "Yost": "GSE123813",
    "Yost_BCC": "GSE123813",
    "YostClonalReplacementTumor2019Adata": "GSE123813",
    "Maynard": "Maynard_NSCLC",
    "Tietscher": "GSE179994",

    # 17 Pan-Cancer Atlas Cohorts
    "Pelka": "GSE178341",
    "Pelka_CRC": "GSE178341",
    "Pelka2021": "GSE178341",
    "PelkaSpatiallyOrganizedMulticellular2021": "GSE178341",
    "PelkaSpatiallyOrganizedMulticellular2021Adata": "GSE178341",

    "Azizi": "GSE114727",
    "Azizi_BRCA": "GSE114727",
    "Azizi2018": "GSE114727",
    "AziziSingleCellMapDiverse2018": "GSE114727",
    "AziziSingleCellMapDiverse2018Adata": "GSE114727",

    "Qian": "E-MTAB-8107",
    "Qian_PanCancer": "E-MTAB-8107",
    "Qian2020": "E-MTAB-8107",
    "QianPancancerBlueprintHeterogeneous2020": "E-MTAB-8107",
    "QianPancancerBlueprintHeterogeneous2020aAdata": "E-MTAB-8107",

    "Cheng": "GSE154763",
    "Cheng_PanCancer": "GSE154763",
    "Cheng2021": "GSE154763",
    "ChengPancancerSinglecellTranscriptional2021": "GSE154763",
    "ChengPancancerSinglecellTranscriptional2021Adata": "GSE154763",

    "Leader": "GSE154826",
    "Leader_NSCLC": "GSE154826",
    "Leader2021": "GSE154826",
    "LeaderSinglecellAnalysisHuman2021": "GSE154826",
    "LeaderSinglecellAnalysisHuman2021Adata": "GSE154826",

    "Kim": "GSE131907",
    "Kim_LUAD": "GSE131907",
    "Kim2020": "GSE131907",
    "KimSinglecellRNASequencing2020": "GSE131907",
    "KimSinglecellRNASequencing2020Adata": "GSE131907",

    "Becker": "GSE201349",
    "Becker_COAD": "GSE201349",
    "Becker2022": "GSE201349",
    "BeckerSinglecellAnalysesDefine2022": "GSE201349",
    "BeckerSinglecellAnalysesDefine2022Adata": "GSE201349",

    "Khaliq": "GSE200997",
    "Khaliq_CC": "GSE200997",
    "Khaliq2022": "GSE200997",
    "KhaliqRefiningColorectalCancer2022": "GSE200997",
    "KhaliqRefiningColorectalCancer2022Adata": "GSE200997",

    "Borcherding": "GSE121638",
    "Borcherding_ccRCC": "GSE121638",
    "Borcherding2021": "GSE121638",
    "BorcherdingMappingImmuneEnvironment2021": "GSE121638",
    "BorcherdingMappingImmuneEnvironment2021Adata": "GSE121638",

    "Sharma": "GSE156625",
    "Sharma_HCC": "GSE156625",
    "Sharma2020": "GSE156625",
    "SharmaOncofetalReprogrammingEndothelial2020": "GSE156625",
    "SharmaOncofetalReprogrammingEndothelial2020Adata": "GSE156625",

    "Lu": "GSE149614",
    "Lu_HCC": "GSE149614",
    "Lu2022": "GSE149614",
    "LuSinglecellAtlasMulticellular2022": "GSE149614",
    "LuSinglecellAtlasMulticellular2022Adata": "GSE149614",

    "Pu": "GSE184362",
    "Pu_PTC": "GSE184362",
    "Pu2021": "GSE184362",
    "PuSinglecellTranscriptomicAnalysis2021": "GSE184362",
    "PuSinglecellTranscriptomicAnalysis2021Adata": "GSE184362",

    "Durante": "GSE139829",
    "Durante_UVM": "GSE139829",
    "Durante2020": "GSE139829",
    "DuranteSinglecellAnalysisReveals2020": "GSE139829",
    "DuranteSinglecellAnalysisReveals2020Adata": "GSE139829",

    "Biermann": "GSE200218",
    "Biermann_BrainMet": "GSE200218",
    "Biermann2022": "GSE200218",
    "BiermannDissectingTreatmentnaiveEcosystem2022": "GSE200218",
    "BiermannDissectingTreatmentnaiveEcosystem2022Adata": "GSE200218",

    "Vazquez": "GSE180661",
    "Vazquez_OV": "GSE180661",
    "Vazquez2022": "GSE180661",
    "VazquezOvarianCancerMutational2022": "GSE180661",
    "VazquezOvarianCancerMutational2022Adata": "GSE180661",

    "Zhang2021": "GSE169246",
    "Zhang_TNBC": "GSE169246",
    "ZhangSinglecellAnalysesReveal2021": "GSE169246",
    "ZhangSinglecellAnalysesReveal2021Adata": "GSE169246",

    "Zhang2022": "GSE215120",
    "Zhang_Myeloid": "GSE215120",
    "ZhangSinglecellAnalysisReveals2022": "GSE215120",
    "ZhangSinglecellAnalysisReveals2022Adata": "GSE215120",
}


def resolve_dataset_id(dataset_id: str) -> str:
    """Resolve an alias or study name to its canonical primary dataset ID."""
    return DATASET_ALIASES.get(dataset_id, dataset_id)


def get_dataset_spec(dataset_id: str) -> Maybe[DatasetSpec]:
    """Retrieve the immutable dataset metadata specification by ID or alias."""
    canonical_id = resolve_dataset_id(dataset_id)
    return Some(DATASET_REGISTRY[canonical_id]) if canonical_id in DATASET_REGISTRY else Nothing


def list_registered_datasets() -> RegisteredDatasets:
    """Return an enriched collection of all registered dataset specifications."""
    return RegisteredDatasets(DATASET_REGISTRY.values())

