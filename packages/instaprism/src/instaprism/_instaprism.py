"""Module for performing deconvolution on gene expression data.

This module includes functions to estimate cell type proportions in bulk RNA-seq data
using the methods of BayesPrism and InstaPrism.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, SupportsIndex

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from collections.abc import Sequence


class _RandomGenerator(Protocol):
    """Protocol for random number generators used in deconvolution algorithms.

    This protocol defines the interface for random number generators that provide
    sampling methods required by the Bayesian deconvolution algorithms. It is
    compatible with numpy.random.Generator and similar implementations.

    The protocol ensures type safety when using random sampling operations in
    functions like bayes_prism, allowing for custom random number generator
    implementations while maintaining compatibility with numpy's random module.
    """

    def dirichlet(self, alpha: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]: ...

    def multinomial(
        self,
        n: int | np.integer | npt.NDArray[np.integer],
        pvals: (
            npt.NDArray[np.bool_]
            | npt.NDArray[np.integer]
            | npt.NDArray[np.floating]
            | Sequence[int]
            | Sequence[float]
        ),
        size: SupportsIndex | Sequence[SupportsIndex] | None = None,
    ) -> npt.NDArray[np.int64]: ...


def deconvolution(
    bulk: npt.NDArray[np.float64],
    reference: npt.NDArray[np.float64],
    n: int,
    eps: float,
) -> npt.NDArray[np.float64]:
    """Performs a deconvolution process on the provided bulk dataset.

    This function implements an iterative computational deconvolution to infer
    proportions or contributions of the reference dataset to the bulk dataset.

    Args:
        bulk (np.array): A 1-dimensional array representing the bulk dataset, where
            each element captures an aggregate data point.
        reference (np.array): A 2-dimensional array where each row represents a
            reference signature, and columns indicate features. The shape of this
            array should be (c, g), where 'c' is the number of reference elements
            and 'g' is the number of features/signatures.
        n (int): The number of iterations to perform for the deconvolution process.
        eps (float): A small epsilon value added to avoid division by zero during the
            normalization process.

    Returns:
        np.array: A 2-dimensional array representing the result of the deconvolution
        process. The output matrix has the same shape as the reference matrix, where
        each entry indicates the inferred proportion of contribution.
    """
    c, g = reference.shape

    # Initialize expression matrix
    b = np.array(
        np.ones(shape=(c, g)) / np.repeat(c * g, c * g).reshape((c, g)),
        dtype=np.float64,
    )

    # Iteration scheme
    for _i in range(n):
        b = reference * np.repeat(b.sum(axis=1), g).reshape((c, g))
        b = b / (np.repeat(b.sum(axis=0), c).reshape((g, c)).T + eps)
        b = b * np.repeat(bulk, c).reshape((g, c)).T
    return b


def insta_prism(
    bulk: npt.NDArray[np.float64],  # shape (1, G) or (G,)
    reference: npt.NDArray[np.float64],  # shape (S, G): S cell types x G genes
    n_iter: int = 1000,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Performs deconvolution of gene expression data.

    This function estimates the contribution of each cell type to the observed
    gene expression data (bulk), using an iterative procedure. It updates the
    probability matrix, cell state gene expression, and cell fractions until
    a specified number of iterations is reached.

    Args:
        bulk: np.ndarray
            Gene expression data for a bulk sample. Can be a 1D array with
            shape (G,) or a 2D array with shape (1, G), where G is the number
            of genes.
        reference: np.ndarray
            Reference gene expression profiles for S cell types. Expected as a
            2D array with shape (S, G), where S is the number of cell types
            and G is the number of genes.
        n_iter: int
            The number of iterations for the iterative estimation process.

    Returns:
        np.ndarray
            Tuple containing three arrays:
            - probability_matrix: Estimates for the probability matrix of dimensions
              reflecting the relationships among genes and cell types.
            - cell_state_gene_expression: Gene expression matrix for each cell state.
            - cell_fractions: Estimated proportions of cell types in the bulk sample.
    """
    reference = reference.T  # (G, S)

    cell_fractions, cell_state_gene_expression, probability_matrix = (
        _initialize_deconvolution_arrays(reference)
    )
    for _ in range(n_iter):
        _update_probability_matrix_inplace(
            reference, cell_fractions, probability_matrix
        )
        _update_cell_state_gene_expression_by_fixpoint_inplace(
            cell_state_gene_expression, bulk, probability_matrix
        )
        _update_cell_fractions_estimate_by_fixpoint_inplace(
            cell_fractions, cell_state_gene_expression
        )
    return probability_matrix, cell_state_gene_expression, cell_fractions


def _update_cell_state_gene_expression_by_fixpoint_inplace(
    cell_state_gene_expression: npt.NDArray[np.float64],
    bulk: npt.NDArray[np.float64],
    probability_matrix: npt.NDArray[np.float64],
) -> None:
    """Updates cell state gene expression matrix inplace.

    This function performs element-wise multiplication of the provided bulk dataset
    with the probability matrix, broadcasting as required, and stores the result
    directly into the provided cell_state_gene_expression matrix. The input is
    modified in-place to improve performance and reduce memory overhead.

    Args:
        cell_state_gene_expression: Matrix that stores the final updated cell state
            gene expression values. The matrix is updated in-place.
        bulk: Array containing bulk data which will be used for updating the
            cell state gene expression matrix.
        probability_matrix: Matrix containing probability values for corresponding
            elements that will be used during the update process.

    """
    np.multiply(bulk[:, np.newaxis], probability_matrix, out=cell_state_gene_expression)


def _update_cell_fractions_estimate_by_fixpoint_inplace(
    cell_fractions: npt.NDArray[np.float64],
    cell_state_gene_expression: npt.NDArray[np.float64],
) -> None:
    """Updates the cell fractions estimate in place.

    Args:
        cell_fractions: A NumPy array that will be updated in place to represent
            the normalized proportions of cell fractions. The values in this array
            are overwritten based on calculations from `cell_state_gene_expression`.
        cell_state_gene_expression: A NumPy 2D array containing gene expression values
            where rows represent genes and columns correspond to different cell states.
            The function calculates the sum across rows for each cell state and uses
            these values to update `cell_fractions`.
    """
    cell_fractions[:] = cell_state_gene_expression.sum(axis=0)
    cell_fractions /= cell_fractions.sum()


def _update_probability_matrix_inplace(
    reference: npt.NDArray[np.float64],
    theta: npt.NDArray[np.float64],
    probability_matrix: npt.NDArray[np.float64],
) -> None:
    """Updates the probability matrix in-place.

    This function multiplies the `reference` and `theta` arrays element-wise and
    updates the `probability_matrix` in place. Afterward, the updated
    `probability_matrix` rows are normalized so that the sum of their elements equals 1.
    It is assumed that no rows in the resulting probability matrix will have
    a sum of zero, thereby avoiding division by zero errors.

    Args:
        reference (np.array): A 2D array representing the reference data.
        theta (np.array): A 2D array with the same shape as `reference`.
        probability_matrix (np.array): A 2D array that will be updated in-place.
    """
    probability_matrix[:] = reference * theta  # (G, S)
    probability_matrix /= probability_matrix.sum(
        axis=1, keepdims=True
    )  # no divide-by-zero assumed


def bayes_prism(
    bulk: np.ndarray,  # shape (1, G) or (G,)
    reference: np.ndarray,  # shape (S, G): S cell types x G genes
    n_iter: int = 1000,
    alpha: float = 1e-4,
    rng: _RandomGenerator | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Performs Bayesian deconvolution of bulk gene expression data.

    This function applies an iterative approach to estimate two key outputs:
    cell state gene expression and cell type fractions. Using Bayesian updates,
    it refines these estimates based on given input data.

    Args:
        bulk (np.ndarray): A 1-dimensional array of bulk gene expression data with
            shape `(1, G)` or `(G,)`, where `G` is the number of genes.
        reference (np.ndarray): A 2-dimensional reference matrix with shape `(S, G)`,
            where `S` is the number of cell types, and `G` is the number of genes.
        n_iter (int): Number of iterations for the Bayesian updates. Defaults to 1000.
        alpha (float): Dirichlet prior parameter for cell fractions. Defaults to 1e-4.
        rng (Generator): A NumPy random number generator instance for control over
            the random sampling process. Defaults to a generator initialized with
            `np.random.default_rng(0)`.

    Returns:
        Tuple[np.ndarray, np.ndarray]: A tuple containing:
            - cell_state_gene_expression (np.ndarray): A matrix representing the
              updated cell state gene expression.
            - cell_fractions (np.ndarray): An array representing the estimated
              proportions of each cell type.
    """
    if rng is None:
        rng_impl: _RandomGenerator = np.random.default_rng(0)
    else:
        rng_impl = rng

    reference = reference.T  # (G, S)

    cell_fractions, cell_state_gene_expression, probability_matrix = (
        _initialize_deconvolution_arrays(reference)
    )
    for _ in range(n_iter):
        _update_probability_matrix_inplace(
            reference, cell_fractions, probability_matrix
        )
        _update_cell_state_gene_expression_by_sampling_inplace(
            cell_state_gene_expression, bulk, probability_matrix, rng_impl
        )
        _update_cell_fractions_estimate_by_sampling_inplace(
            cell_fractions, cell_state_gene_expression, alpha, rng_impl
        )
    return cell_state_gene_expression, cell_fractions


def _initialize_deconvolution_arrays(
    reference: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Initializes and returns arrays used for the deconvolution process.

    This function creates and initializes three arrays necessary for the
    deconvolution process. The cell fractions array is initialized with
    equal values summing to 1, while the probability matrix and
    cell state gene expression matrix are pre-allocated with appropriate
    dimensions based on the reference matrix.

    Args:
        reference (numpy.ndarray): A 2D array representing the reference data,
            where rows correspond to genes, and columns correspond to cells.

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: A 1D array representing cell fractions, initialized
              such that each cell has an equal fraction.
            - numpy.ndarray: A 2D array for cell state gene expression, pre-allocated
              based on the dimensions of the reference matrix.
            - numpy.ndarray: A 2D array for the probability matrix, pre-allocated
              based on the reference dimensions.
    """
    gene_number, cell_number = reference.shape  # (G, S)
    cell_fractions = np.full(cell_number, 1.0 / cell_number)
    probability_matrix = np.empty((gene_number, cell_number), dtype=np.float64)
    cell_state_gene_expression = np.empty((gene_number, cell_number), dtype=np.float64)
    return cell_fractions, cell_state_gene_expression, probability_matrix


def _update_cell_fractions_estimate_by_sampling_inplace(
    cell_fractions: npt.NDArray[np.float64],
    cell_state_gene_expression: npt.NDArray[np.float64],
    alpha: float,
    rng: _RandomGenerator,
) -> None:
    """Updates the cell fractions estimate in place by sampling.

    This function modifies the input cell fractions directly, updating its values based
    on a Dirichlet distribution. The parameters of the distribution are derived from
    the sum of the cell state gene expression and the provided alpha value.

    Args:
        cell_fractions: numpy array, modified in place to reflect updated cell fraction
            estimates.
        cell_state_gene_expression: numpy array representing gene expression levels for
            different cell states.
        alpha: float, the concentration parameter added to the summed gene expression to
            form the parameters of the Dirichlet distribution.
        rng: numpy.random.Generator, used to sample from the Dirichlet distribution.
    """
    cell_fractions[:] = rng.dirichlet(cell_state_gene_expression.sum(axis=0) + alpha)


def _update_cell_state_gene_expression_by_sampling_inplace(
    cell_state_gene_expression: npt.NDArray[np.float64],
    bulk: npt.NDArray[np.float64],
    probability_matrix: npt.NDArray[np.float64],
    rng: _RandomGenerator,
) -> None:
    """Updates the `cell_state_gene_expression` matrix by sampling.

    This function modifies the `cell_state_gene_expression` in place based
    on the values of the `bulk` vector and the `probability_matrix`.
    The sampling is performed row-wise in the `cell_state_gene_expression` matrix
    using the multinomial distribution with parameters provided in `bulk` and
    `probability_matrix` for each row. The randomness of the sampling process is
    managed by the provided random number generator `rng`.

    Args:
        cell_state_gene_expression (np.ndarray): A 2D array representing
            the cell state gene expression to be updated in place.
        bulk (np.ndarray): A 1D array containing the total values for each gene to use
            in the multinomial sampling process.
        probability_matrix (np.ndarray): A 2D array representing the probability matrix
            used for the multinomial sampling, where each row corresponds to a gene,
            and the columns correspond to probabilities for each cell state.
        rng (np.random.Generator): A random number generator instance used for
            multinomial sampling to ensure reproducibility across runs.
    """
    for g in range(cell_state_gene_expression.shape[0]):
        cell_state_gene_expression[g, :] = rng.multinomial(
            bulk[g], probability_matrix[g, :]
        )


def _normalize_rows_to_stochastic(
    x: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Normalizes the rows of a 2D array to make them stochastic.

    This function ensures that the rows of the input 2D array sum to 1, effectively
    turning the array into a row-stochastic matrix. Each row is divided by the sum
    of its elements.

    Args:
        x: A 2D NumPy array of integers or floats. The array must have two dimensions
           (ndim=2).

    Returns:
        A 2D NumPy array where each row has been normalized such that the sum of
        elements in each row equals 1. The output array will have the same shape as
        the input array.

    Raises:
        AssertionError: If the input array is not a 2D array (ndim != 2).
    """
    return np.asarray(x / x.sum(axis=1).reshape((x.shape[0], 1)), dtype=np.float64)


def _calculate_fractions(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Calculates row-wise normalized fractions from the input array.

    This function takes a NumPy array, computes the sum along its rows, and normalizes
    these sums to create a fraction representation. The output is a 1D array of
    fractions representing the relative contributions of row-wise sums to
    the overall sum.

    Args:
        x: A NumPy array containing integer or float values. The input array should
           have at least two dimensions.

    Returns:
        A 1D NumPy array of float values representing the normalized fractions derived
        from the input array.
    """
    theta = x.sum(axis=1)
    return np.asarray(theta / theta.sum(), dtype=np.float64)