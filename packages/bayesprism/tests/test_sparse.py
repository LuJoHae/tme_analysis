import pytest
import numpy as np
import scipy.sparse as sp
import torch
from bayesprism.preprocessing import norm_to_one, collapse, filter_bulk_outlier


def test_sparse_dense_norm_to_one_parity() -> None:
    mat = np.array([[10, 0, 20], [0, 5, 15]], dtype=np.float32)
    mat_tensor = torch.from_numpy(mat)

    res_dense = norm_to_one(mat_tensor, pseudo_min=1e-8)

    # Convert sparse to dense before norm_to_one
    sparse_mat = sp.csr_matrix(mat)
    mat_from_sparse = torch.from_numpy(sparse_mat.toarray())
    res_sparse = norm_to_one(mat_from_sparse, pseudo_min=1e-8)

    torch.testing.assert_close(res_dense, res_sparse)


def test_collapse_parity() -> None:
    ref = np.array([
        [10, 20, 30],
        [5, 15, 25],
        [1, 2, 3]
    ], dtype=np.float32)
    labels = ("typeA", "typeA", "typeB")

    res_dense, labels_dense = collapse(ref, labels)
    res_sparse, labels_sparse = collapse(sp.csr_matrix(ref), labels)

    assert labels_dense == labels_sparse == ("typeA", "typeB")
    torch.testing.assert_close(res_dense, res_sparse)
