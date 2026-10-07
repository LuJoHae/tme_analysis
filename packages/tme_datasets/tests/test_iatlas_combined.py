"""Unit tests for load_combined_iatlas_cohorts and iAtlas cohort grouping."""

from returns.result import Failure, Success
from tme_datasets import IATLAS_COMBINED_GROUPS, load_combined_iatlas_cohorts


def test_iatlas_combined_groups_structure() -> None:
    assert "melanoma" in IATLAS_COMBINED_GROUPS
    assert "rcc" in IATLAS_COMBINED_GROUPS
    assert "pancancer" in IATLAS_COMBINED_GROUPS
    assert "all" in IATLAS_COMBINED_GROUPS

    assert len(IATLAS_COMBINED_GROUPS["melanoma"]) == 4
    assert "Hugo-iAtlas" in IATLAS_COMBINED_GROUPS["melanoma"]
    assert "Liu-iAtlas" in IATLAS_COMBINED_GROUPS["melanoma"]

    assert len(IATLAS_COMBINED_GROUPS["rcc"]) == 2
    assert "McDermott-iAtlas" in IATLAS_COMBINED_GROUPS["rcc"]
    assert "Choueiri-iAtlas" in IATLAS_COMBINED_GROUPS["rcc"]


def test_load_combined_iatlas_rcc_cohorts() -> None:
    res = load_combined_iatlas_cohorts("rcc")
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.n_obs > 200
    assert adata.n_vars > 1000
    assert "dataset_id" in adata.obs.columns
    assert set(adata.obs["dataset_id"].unique()) == {"McDermott-iAtlas", "Choueiri-iAtlas"}
    assert "response_binary" in adata.obs.columns


def test_load_combined_iatlas_custom_list() -> None:
    res = load_combined_iatlas_cohorts(["Hugo-iAtlas", "Gide-iAtlas"])
    assert isinstance(res, Success)
    adata = res.unwrap()
    assert adata.n_obs == 27 + 91
    assert set(adata.obs["dataset_id"].unique()) == {"Hugo-iAtlas", "Gide-iAtlas"}


def test_load_combined_iatlas_invalid_inputs() -> None:
    err_res1 = load_combined_iatlas_cohorts("invalid_group_name")
    assert isinstance(err_res1, Failure)
    assert "Unknown iAtlas cohort or grouping" in err_res1.failure()

    err_res2 = load_combined_iatlas_cohorts([])
    assert isinstance(err_res2, Failure)
    assert "No cohort IDs provided" in err_res2.failure()
