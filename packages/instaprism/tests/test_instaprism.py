from __future__ import annotations

from typing import TYPE_CHECKING, SupportsIndex

if TYPE_CHECKING:
    from collections.abc import Sequence

import numpy as np

from instaprism._instaprism import (
    _calculate_fractions,
    _initialize_deconvolution_arrays,
    _normalize_rows_to_stochastic,
    _update_cell_fractions_estimate_by_fixpoint_inplace,
    _update_cell_fractions_estimate_by_sampling_inplace,
    _update_cell_state_gene_expression_by_fixpoint_inplace,
    _update_cell_state_gene_expression_by_sampling_inplace,
    _update_probability_matrix_inplace,
    bayes_prism,
    deconvolution,
    insta_prism,
)


class FakeRng:
    def __init__(self) -> None:
        self.dirichlet_calls: list[np.ndarray] = []
        self.multinomial_calls: list[tuple[int, np.ndarray]] = []

    def dirichlet(self, alpha: np.ndarray) -> np.ndarray:
        self.dirichlet_calls.append(np.array(alpha, dtype=np.float64))
        size = alpha.shape[0]
        return np.asarray(
            np.arange(1, size + 1, dtype=np.float64)
            / np.sum(np.arange(1, size + 1, dtype=np.float64)),
            dtype=np.float64,
        )

    def multinomial(
        self,
        n: int | np.integer | np.ndarray,
        pvals: (np.ndarray | Sequence[int] | Sequence[float]),
        size: SupportsIndex | Sequence[SupportsIndex] | None = None,
    ) -> np.ndarray:
        n_int = int(n)
        pvals_arr = np.asarray(pvals, dtype=np.float64)
        self.multinomial_calls.append((n_int, pvals_arr))
        out = np.zeros_like(pvals_arr, dtype=np.int64)
        out[int(np.argmax(pvals_arr))] = n_int
        return out


def test_deconvolution_returns_expected_shape_and_matches_manual_one_step() -> None:
    bulk = np.array([2.0, 4.0], dtype=np.float64)
    reference = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    result = deconvolution(bulk, reference, n=1, eps=1e-12)
    assert result.shape == reference.shape


def test_deconvolution_zero_iterations_returns_initial_matrix() -> None:
    bulk = np.array([1.0, 2.0], dtype=np.float64)
    reference = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)

    result = deconvolution(bulk, reference, n=0, eps=1e-12)

    expected = np.full(reference.shape, 1.0 / reference.size, dtype=np.float64)
    np.testing.assert_allclose(result, expected)


def test_insta_prism_zero_iterations_initializes_arrays() -> None:
    bulk = np.array([5.0, 6.0], dtype=np.float64)
    reference = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)

    probability_matrix, cell_state_gene_expression, cell_fractions = insta_prism(
        bulk, reference, n_iter=0
    )

    assert probability_matrix.shape == (2, 2)
    assert cell_state_gene_expression.shape == (2, 2)
    np.testing.assert_allclose(cell_fractions, np.array([0.5, 0.5], dtype=np.float64))


def test_initialize_deconvolution_arrays_creates_expected_shapes_and_values() -> None:
    reference = np.ones((3, 2), dtype=np.float64)

    cell_fractions, cell_state_gene_expression, probability_matrix = (
        _initialize_deconvolution_arrays(reference)
    )

    np.testing.assert_allclose(cell_fractions, np.array([0.5, 0.5], dtype=np.float64))
    assert cell_state_gene_expression.shape == (3, 2)
    assert probability_matrix.shape == (3, 2)
    assert cell_fractions.dtype == np.float64
    assert cell_state_gene_expression.dtype == np.float64
    assert probability_matrix.dtype == np.float64


def test_update_probability_matrix_inplace_normalizes_rows() -> None:
    reference = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    theta = np.array([[2.0, 1.0], [1.0, 3.0]], dtype=np.float64)
    probability_matrix = np.empty_like(reference)

    _update_probability_matrix_inplace(reference, theta, probability_matrix)


def test_update_cell_state_gene_expression_by_fixpoint_inplace_multiplies_bulk() -> (
    None
):
    bulk = np.array([2.0, 3.0], dtype=np.float64)
    probability_matrix = np.array([[0.25, 0.75], [0.6, 0.4]], dtype=np.float64)
    cell_state_gene_expression = np.empty((2, 2), dtype=np.float64)

    _update_cell_state_gene_expression_by_fixpoint_inplace(
        cell_state_gene_expression,
        bulk,
        probability_matrix,
    )

    expected = np.array([[0.5, 1.5], [1.8, 1.2]], dtype=np.float64)
    np.testing.assert_allclose(cell_state_gene_expression, expected)


def test_update_cell_fractions_estimate_by_fixpoint_inplace_normalizes_sums() -> None:
    cell_state_gene_expression = np.array(
        [[1.0, 3.0], [2.0, 4.0]],
        dtype=np.float64,
    )
    cell_fractions = np.empty(2, dtype=np.float64)

    _update_cell_fractions_estimate_by_fixpoint_inplace(
        cell_fractions,
        cell_state_gene_expression,
    )


def test_update_cell_fractions_estimate_by_sampling_inplace_uses_rng() -> None:
    cell_state_gene_expression = np.array(
        [[1.0, 2.0], [3.0, 4.0]],
        dtype=np.float64,
    )
    cell_fractions = np.zeros(2, dtype=np.float64)
    rng = FakeRng()

    _update_cell_fractions_estimate_by_sampling_inplace(
        cell_fractions,
        cell_state_gene_expression,
        alpha=0.1,
        rng=rng,
    )

    expected = np.array([1.0 / 3.0, 2.0 / 3.0], dtype=np.float64)
    np.testing.assert_allclose(cell_fractions, expected)
    assert len(rng.dirichlet_calls) == 1
    np.testing.assert_allclose(rng.dirichlet_calls[0], np.array([4.1, 6.1]))


def test_update_cell_state_gene_expression_by_sampling_inplace_uses_rng_per_row() -> (
    None
):
    bulk = np.array([5, 7], dtype=np.float64)
    probability_matrix = np.array([[0.1, 0.9], [0.8, 0.2]], dtype=np.float64)
    cell_state_gene_expression = np.empty((2, 2), dtype=np.float64)
    rng = FakeRng()

    _update_cell_state_gene_expression_by_sampling_inplace(
        cell_state_gene_expression,
        bulk,
        probability_matrix,
        rng,
    )

    expected = np.array([[0.0, 5.0], [7.0, 0.0]], dtype=np.float64)
    np.testing.assert_allclose(cell_state_gene_expression, expected)
    assert len(rng.multinomial_calls) == 2
    assert rng.multinomial_calls[0][0] == 5
    assert rng.multinomial_calls[1][0] == 7


def test_normalize_rows_to_stochastic_normalizes_each_row() -> None:
    x = np.array([[1.0, 1.0], [2.0, 6.0]], dtype=np.float64)

    result = _normalize_rows_to_stochastic(x)

    expected = np.array([[0.5, 0.5], [0.25, 0.75]], dtype=np.float64)
    np.testing.assert_allclose(result, expected)
    np.testing.assert_allclose(result.sum(axis=1), np.ones(2))


def test_calculate_fractions_normalizes_row_sums() -> None:
    x = np.array([[1.0, 1.0], [2.0, 6.0]], dtype=np.float64)

    result = _calculate_fractions(x)

    expected = np.array([0.2, 0.8], dtype=np.float64)
    np.testing.assert_allclose(result, expected)
    np.testing.assert_allclose(result.sum(), 1.0)


def test_bayes_prism_zero_iterations_returns_initialized_arrays() -> None:
    bulk = np.array([5, 7], dtype=np.int64)
    reference = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)

    cell_state_gene_expression, cell_fractions = bayes_prism(
        bulk,
        reference,
        n_iter=0,
        rng=FakeRng(),
    )

    assert cell_state_gene_expression.shape == (2, 2)
    np.testing.assert_allclose(cell_fractions, np.array([0.5, 0.5], dtype=np.float64))


def test_insta_prism_is_deterministic_for_one_iteration() -> None:
    bulk = np.array([2.0, 3.0], dtype=np.float64)
    reference = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)

    probability_matrix, cell_state_gene_expression, cell_fractions = insta_prism(
        bulk,
        reference,
        n_iter=1,
    )

    assert probability_matrix.shape == (2, 2)
    assert cell_state_gene_expression.shape == (2, 2)
    assert cell_fractions.shape == (2,)
    np.testing.assert_allclose(probability_matrix.sum(axis=1), np.ones(2))
    np.testing.assert_allclose(cell_fractions.sum(), 1.0)


def test_bayes_prism_is_deterministic_with_fake_rng() -> None:
    bulk = np.array([2, 3], dtype=np.int64)
    reference = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    rng = FakeRng()

    cell_state_gene_expression, cell_fractions = bayes_prism(
        bulk,
        reference,
        n_iter=1,
        alpha=0.1,
        rng=rng,
    )

    assert cell_state_gene_expression.shape == (2, 2)
    np.testing.assert_allclose(cell_fractions.sum(), 1.0)
    assert len(rng.multinomial_calls) == 2
    assert len(rng.dirichlet_calls) == 1