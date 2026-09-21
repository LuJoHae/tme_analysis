"""Unit tests for dataset registry and specs."""

from returns.maybe import Some
from tme_datasets.registry import get_dataset_spec, list_registered_datasets
from tme_datasets.types import Modality


def test_list_registered_datasets() -> None:
    specs = list_registered_datasets()
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
