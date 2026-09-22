import polars as pl
import pytest
from returns.maybe import Maybe, Nothing, Some

from tme_datasets.models import RegisteredDatasets
from tme_datasets.registry import get_dataset_spec, list_registered_datasets
from tme_datasets.types import Modality


def test_list_registered_datasets() -> None:
    specs = list_registered_datasets()
    assert isinstance(specs, tuple)
    assert isinstance(specs, RegisteredDatasets)
    assert len(specs) >= 20
    ids = [s.id for s in specs]
    assert "GSE120575" in ids
    assert "Hugo-iAtlas" in ids
    assert "Rosenberg-iAtlas" in ids
    assert "Auslander" in ids


def test_get_dataset_spec_found() -> None:
    spec_maybe = get_dataset_spec("GSE120575")
    assert isinstance(spec_maybe, Some)
    spec = spec_maybe.unwrap()
    assert spec.id == "GSE120575"
    assert spec.modality == Modality.SINGLE_CELL
    assert spec.has_response_labels is True


def test_get_dataset_spec_not_found() -> None:
    spec_maybe = get_dataset_spec("NON_EXISTENT_ID")
    assert not isinstance(spec_maybe, Some)


def test_registered_datasets_dual_indexing() -> None:
    specs = list_registered_datasets()
    # Positional indexing
    first_spec = specs[0]
    assert first_spec.id == specs.ids()[0]

    # Slicing returns RegisteredDatasets instance
    sub_specs = specs[0:3]
    assert isinstance(sub_specs, RegisteredDatasets)
    assert len(sub_specs) == 3

    # String ID indexing
    sade = specs["GSE120575"]
    assert sade.id == "GSE120575"
    assert sade.modality == Modality.SINGLE_CELL

    # String ID KeyError
    with pytest.raises(KeyError, match="not found"):
        _ = specs["NON_EXISTENT_DATASET_ID"]


def test_registered_datasets_contains_and_get() -> None:
    specs = list_registered_datasets()
    assert "GSE120575" in specs
    assert "NON_EXISTENT" not in specs
    assert specs[0] in specs

    # .get() method returns Maybe
    assert isinstance(specs.get("GSE120575"), Some)
    assert specs.get("GSE120575").unwrap().id == "GSE120575"
    assert specs.get("NON_EXISTENT") is Nothing


def test_registered_datasets_filter() -> None:
    specs = list_registered_datasets()
    sc_specs = specs.filter(modality=Modality.SINGLE_CELL)
    assert isinstance(sc_specs, RegisteredDatasets)
    assert len(sc_specs) > 0
    assert all(s.modality == Modality.SINGLE_CELL for s in sc_specs)

    # Filter with cancer type
    melanoma = specs.filter(cancer_type="Melanoma")
    assert all(s.cancer_type.lower() == "melanoma" for s in melanoma)


def test_registered_datasets_attribute_accessors() -> None:
    specs = list_registered_datasets()

    # Direct list accessors
    ids = specs.ids()
    titles = specs.titles()
    modalities = specs.modalities()
    cancer_types = specs.cancer_types()
    platforms = specs.platforms()
    has_resp = specs.has_response_labels()

    assert len(ids) == len(specs)
    assert len(titles) == len(specs)
    assert len(modalities) == len(specs)
    assert len(cancer_types) == len(specs)
    assert len(platforms) == len(specs)
    assert len(has_resp) == len(specs)

    # Optional attribute accessors: wrapped (default) vs unwrapped
    organs_wrapped = specs.organs()
    assert all(isinstance(o, Maybe) for o in organs_wrapped)

    organs_unwrapped = specs.organs(unwrapped=True)
    assert any(o is not None for o in organs_unwrapped)

    counts_wrapped = specs.n_samples_or_cells()
    counts_unwrapped = specs.n_samples_or_cells(unwrapped=True)
    assert len(counts_wrapped) == len(counts_unwrapped) == len(specs)

    urls_wrapped = specs.raw_source_urls()
    urls_unwrapped = specs.raw_source_urls(unwrapped=True)
    assert len(urls_wrapped) == len(urls_unwrapped) == len(specs)

    checksums = specs.checksums()
    assert len(checksums) == len(specs)

    qc_specs = specs.qc_specs()
    assert len(qc_specs) == len(specs)

    local_paths = specs.local_paths()
    assert len(local_paths) == len(specs)


def test_registered_datasets_to_polars() -> None:
    specs = list_registered_datasets()
    df = specs.to_polars()

    assert isinstance(df, pl.DataFrame)
    assert len(df) == len(specs)

    expected_cols = {
        "id",
        "title",
        "modality",
        "cancer_type",
        "platform",
        "organ",
        "has_response_labels",
        "n_samples_or_cells",
        "raw_source_url",
        "checksum",
        "has_qc_spec",
        "local_path",
    }
    assert expected_cols.issubset(set(df.columns))
    assert "GSE120575" in df["id"].to_list()

    # to_dataframe alias
    df_alias = specs.to_dataframe()
    assert df.shape == df_alias.shape

    # Check repr
    assert "<RegisteredDatasets len=" in repr(specs)
