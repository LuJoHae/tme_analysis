from typing import Any, List

import pandas as pd
from datalair import Lair
from matplotlib import pyplot as plt
from pandas import DataFrame, Series

from ici_datasets.bagaev_datasets import get_pan_melanoma_bagaev_dataset, Signature
from ici_datasets.other_datasets import get_vanallen_bagaev_signature
from ml_pipelines import rf_pipeline, plot_subtype_responser_fractions


def create_dataset(cohorts, lair) -> tuple[DataFrame, Any, Series, DataFrame]:
    annotations, signature = get_pan_melanoma_bagaev_dataset(lair)
    mask_has_response = annotations["Response"].isin(["R", "NR"])
    annotations = annotations.loc[mask_has_response]
    signature = signature.loc[mask_has_response]

    if cohorts is not None:
        annotations = annotations[
            annotations["Cohort"].isin(cohorts)]
    signature = signature.loc[annotations.index]
    data = annotations.loc[
        annotations["Response"].isin(["R", "NR"])
    ][["MFP", "Response", "Cohort", "Therapy", "Cohort_group"]]

    x = signature.loc[data.index]
    y = data["Response"]
    m = data[["Cohort"]]
    return data, m, x, y


def pipeline_individual_cohort(
    lair: Lair,
    cohorts: List[str] | None = None,
    n_components: int = 2,
    n_estimators: int = 50,
    use_cohort: bool = True,
    max_depth: int | None = None
):
    data, m, x, y = create_dataset(cohorts, lair)
    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    plot_subtype_responser_fractions(data, axes[0:2], title="ALL")
    if use_cohort:
        y_proba, axes, clf = rf_pipeline(x, y,
                                         categorical_metadata=m,
                                         axes=axes[2:5],
                                         n_components=n_components,
                                         n_estimators=n_estimators,
                                         max_depth=max_depth
                                         )
    else:
        y_proba, axes, clf = rf_pipeline(x, y,
                                         axes=axes[2:5],
                                         n_components=n_components,
                                         n_estimators=n_estimators,
                                         max_depth=max_depth
                                         )
    data["y_proba"] = y_proba
    return data, fig, axes, clf

def pipeline_all(lair):
    annotations, _ = get_pan_melanoma_bagaev_dataset(lair)
    cohorts = [
        'Lauss', 'Augustine', 'Ulloa-Montoya', 'Hugo', 'Nathanson', 'VanAllen', 'Riaz', 'Liang', 'Liu', 'Auslander', 'Gide'
    ]
    fig, axes_2d = plt.subplots(len(cohorts), 4, figsize=(20, 4*len(cohorts)))
    for axes, cohort_name in zip(axes_2d, cohorts, strict=True):
        print(cohort_name)
        data, m, x, y = create_dataset([cohort_name])
        assert len(data) > 0, cohort_name
        plot_subtype_responser_fractions(data, axes[0:2], title=cohort_name)
        _ = rf_pipeline(x, y, axes=axes[2:4])
    return fig

def pipelin_vanallen_via_freeman(lair):
    ds = Signature()
    lair.safe_derive(ds)
    filepaths = lair.get_dataset_filepaths(ds)
    geneset_order = pd.read_csv(filepaths["gene_signatures_order.tsv"], header=None)

    vanallen_signatures_scaled = get_vanallen_bagaev_signature(lair)

    annotations, _ = get_pan_melanoma_bagaev_dataset(lair)
    annotations_vanallen = annotations[annotations["Cohort"]=="VanAllen"]
    data = annotations_vanallen.loc[
        annotations_vanallen["Response"].isin(["R", "NR"])
    ][["MFP", "Response"]]

    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    data.index = data.index.str.capitalize()
    common_patients = sorted(list(set.intersection(*[set(data.index), set(vanallen_signatures_scaled.index)])))

    plot_subtype_responser_fractions(data, axes[0:2], title="VanAllenViaFreeman")

    x = vanallen_signatures_scaled.loc[common_patients, geneset_order[0]]
    y = data.loc[common_patients, "Response"]

    fig, _ = rf_pipeline(x, y, axes=axes[2:4])
    return fig

