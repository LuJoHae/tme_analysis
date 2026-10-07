import polars as pl
import pytest
from returns.maybe import Maybe, Nothing, Some

from tme_datasets.models import (
    PreprocessedDatasetSpec,
    PreprocessedDatasets,
    RegisteredDatasets,
)
from tme_datasets.registry import (
    get_dataset_spec,
    list_preprocessed_datasets,
    list_registered_datasets,
    query_preprocessed_datasets,
)
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
    is_dys = specs.is_dysfunctional()

    assert len(ids) == len(specs)
    assert len(titles) == len(specs)
    assert len(modalities) == len(specs)
    assert len(cancer_types) == len(specs)
    assert len(platforms) == len(specs)
    assert len(has_resp) == len(specs)
    assert len(is_dys) == len(specs)

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
        "is_dysfunctional",
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


def test_list_preprocessed_datasets() -> None:
    """Verify list_preprocessed_datasets scans and finds on-disk H5AD files."""
    prep = list_preprocessed_datasets()
    assert isinstance(prep, PreprocessedDatasets)
    assert len(prep) >= 5
    ids = prep.ids()
    assert "GSE120575" in ids
    # By default, dysfunctional datasets such as Hugo-iAtlas are excluded
    assert "Hugo-iAtlas" not in ids

    # When exclude_dysfunctional=False, dysfunctional datasets are included
    prep_with_dys = list_preprocessed_datasets(exclude_dysfunctional=False)
    assert len(prep_with_dys) >= 10
    assert "Hugo-iAtlas" in prep_with_dys.ids()

    # Validate model invariants
    for spec in prep:
        assert isinstance(spec, PreprocessedDatasetSpec)
        assert spec.h5ad_path.exists()
        assert spec.h5ad_path.is_file()
        assert spec.file_size_mb > 0
        assert spec.is_dysfunctional is False


def test_list_preprocessed_datasets_filters() -> None:
    """Verify filtering by modality, cancer type, and response labels."""
    # Modality filter with Enum and string
    sc_enum = list_preprocessed_datasets(modality=Modality.SINGLE_CELL)
    sc_str = list_preprocessed_datasets(modality="single_cell")
    assert len(sc_enum) == len(sc_str)
    assert len(sc_enum) >= 4
    assert all(s.modality == Modality.SINGLE_CELL for s in sc_enum)
    assert "GSE120575" in sc_enum.ids()

    bulk_prep = list_preprocessed_datasets(modality=Modality.BULK_RNA)
    assert all(s.modality == Modality.BULK_RNA for s in bulk_prep)
    assert "Anders-iAtlas" in bulk_prep.ids()
    assert "Hugo-iAtlas" not in bulk_prep.ids()

    bulk_with_dys = list_preprocessed_datasets(modality=Modality.BULK_RNA, exclude_dysfunctional=False)
    assert "Hugo-iAtlas" in bulk_with_dys.ids()

    # Cancer type filter (case-insensitive)
    mel_prep = list_preprocessed_datasets(cancer_type="melanoma")
    assert len(mel_prep) >= 2
    assert all(s.cancer_type.lower() == "melanoma" for s in mel_prep)
    mel_with_dys = list_preprocessed_datasets(cancer_type="melanoma", exclude_dysfunctional=False)
    assert len(mel_with_dys) >= 4

    # Response labels filter
    resp_prep = list_preprocessed_datasets(has_response=True)
    assert all(s.has_response_labels is True for s in resp_prep)
    assert "GSE120575" in resp_prep.ids()


def test_preprocessed_datasets_dual_indexing_and_accessors() -> None:
    """Verify indexing by int, slice, ID string, and list accessors."""
    prep = list_preprocessed_datasets()

    # Index by position
    first = prep[0]
    assert isinstance(first, PreprocessedDatasetSpec)
    assert first.id == prep.ids()[0]

    # Slicing returns PreprocessedDatasets
    sub = prep[0:2]
    assert isinstance(sub, PreprocessedDatasets)
    assert len(sub) == 2

    # Index by string ID
    sade = prep["GSE120575"]
    assert sade.id == "GSE120575"
    assert sade.file_size_mb > 0

    with pytest.raises(KeyError, match="not found"):
        _ = prep["NON_EXISTENT_ID"]

    # Contains and get
    assert "GSE120575" in prep
    assert "NON_EXISTENT_COHORT" not in prep

    spec_some = prep.get("GSE120575")
    assert isinstance(spec_some, Some)
    assert spec_some.unwrap().id == "GSE120575"
    assert prep.get("NON_EXISTENT") is Nothing

    # List accessors
    assert len(prep.ids()) == len(prep)
    assert len(prep.titles()) == len(prep)
    assert len(prep.modalities()) == len(prep)
    assert len(prep.cancer_types()) == len(prep)
    assert len(prep.platforms()) == len(prep)
    assert len(prep.paths()) == len(prep)
    assert len(prep.file_sizes_mb()) == len(prep)
    assert len(prep.organs()) == len(prep)
    assert len(prep.organs(unwrapped=True)) == len(prep)
    assert len(prep.n_samples_or_cells()) == len(prep)
    assert len(prep.n_samples_or_cells(unwrapped=True)) == len(prep)
    assert len(prep.has_response_labels()) == len(prep)
    assert len(prep.is_dysfunctional()) == len(prep)
    assert len(prep.has_qc_specs()) == len(prep)


def test_preprocessed_datasets_to_polars() -> None:
    """Verify conversion of PreprocessedDatasets to Polars DataFrame."""
    prep = list_preprocessed_datasets()
    df = prep.to_polars()

    assert isinstance(df, pl.DataFrame)
    assert len(df) == len(prep)

    expected_cols = {
        "id",
        "title",
        "modality",
        "cancer_type",
        "platform",
        "organ",
        "has_response_labels",
        "is_dysfunctional",
        "n_samples_or_cells",
        "h5ad_path",
        "file_size_mb",
        "has_qc_spec",
    }
    assert expected_cols.issubset(set(df.columns))

    df_alias = prep.to_dataframe()
    assert df.shape == df_alias.shape

    # Repr
    assert "<PreprocessedDatasets len=" in repr(prep)


def test_empty_preprocessed_datasets_to_polars() -> None:
    """Verify that an empty PreprocessedDatasets produces a DataFrame with typed columns."""
    empty = PreprocessedDatasets([])
    df = empty.to_polars()
    assert isinstance(df, pl.DataFrame)
    assert len(df) == 0
    assert df.schema["id"] == pl.String
    assert df.schema["has_response_labels"] == pl.Boolean
    assert df.schema["is_dysfunctional"] == pl.Boolean
    assert df.schema["file_size_mb"] == pl.Float64
    assert df.schema["has_qc_spec"] == pl.Boolean


def test_registered_datasets_preprocessed_method() -> None:
    """Verify RegisteredDatasets.preprocessed() filtering."""
    specs = list_registered_datasets()
    prep = specs.preprocessed()
    assert isinstance(prep, PreprocessedDatasets)
    assert len(prep) == len(list_preprocessed_datasets())

    # Preprocessed on a filtered subset
    melanoma_specs = specs.filter(cancer_type="Melanoma")
    melanoma_prep = melanoma_specs.preprocessed()
    assert all(s.cancer_type.lower() == "melanoma" for s in melanoma_prep)


def test_query_preprocessed_datasets() -> None:
    """Verify high-level query_preprocessed_datasets API."""
    # Default returns PreprocessedDatasets (excluding dysfunctional)
    all_prep = query_preprocessed_datasets()
    assert isinstance(all_prep, PreprocessedDatasets)
    assert len(all_prep) >= 5
    assert "Hugo-iAtlas" not in all_prep.ids()

    # With exclude_dysfunctional=False
    all_prep_with_dys = query_preprocessed_datasets(exclude_dysfunctional=False)
    assert len(all_prep_with_dys) >= 10
    assert "Hugo-iAtlas" in all_prep_with_dys.ids()

    # as_polars=True returns pl.DataFrame
    df = query_preprocessed_datasets(as_polars=True)
    assert isinstance(df, pl.DataFrame)
    assert len(df) == len(all_prep)

    # Filter by specific dataset IDs
    subset = query_preprocessed_datasets(["GSE120575", "Maynard_NSCLC"])
    assert len(subset) == 2
    assert set(subset.ids()) == {"GSE120575", "Maynard_NSCLC"}

    # Filter with aliases
    alias_subset = query_preprocessed_datasets(["SadeFeldman"])
    assert len(alias_subset) == 1
    assert alias_subset.ids()[0] == "GSE120575"

    # Combined filters
    mel_sc_df = query_preprocessed_datasets(
        modality="single_cell",
        cancer_type="Melanoma",
        as_polars=True,
    )
    assert isinstance(mel_sc_df, pl.DataFrame)
    assert len(mel_sc_df) >= 2
    assert "GSE120575" in mel_sc_df["id"].to_list()


def test_dysfunctional_flagging() -> None:
    """Verify that dysfunctional datasets are properly flagged and filtered across APIs."""
    from tme_datasets.query import list_datasets

    # Test get_dataset_spec for functional vs dysfunctional cohorts
    func_spec = get_dataset_spec("GSE120575").unwrap()
    assert func_spec.is_dysfunctional is False

    dys_sc_spec = get_dataset_spec("GSE114727").unwrap()
    assert dys_sc_spec.is_dysfunctional is True

    dys_bulk_spec = get_dataset_spec("Hugo-iAtlas").unwrap()
    assert dys_bulk_spec.is_dysfunctional is True

    # Test list_datasets exclude_dysfunctional behavior
    clean_datasets = list_datasets(exclude_dysfunctional=True)
    all_datasets = list_datasets(exclude_dysfunctional=False)
    assert len(all_datasets) - len(clean_datasets) == 44

    # Single-cell breakdown
    sc_clean = list_datasets(modality=Modality.SINGLE_CELL, exclude_dysfunctional=True)
    sc_all = list_datasets(modality=Modality.SINGLE_CELL, exclude_dysfunctional=False)
    assert len(sc_all) - len(sc_clean) == 36
    assert "GSE114727" not in sc_clean.ids()
    assert "GSE114727" in sc_all.ids()
    assert "GSE120575" in sc_clean.ids()

    # RegisteredDatasets filter method
    reg_clean = list_registered_datasets().filter(exclude_dysfunctional=True)
    reg_all = list_registered_datasets().filter(exclude_dysfunctional=False)
    assert len(reg_all) - len(reg_clean) == 44
