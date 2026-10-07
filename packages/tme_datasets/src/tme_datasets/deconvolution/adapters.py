"""Adapters to convert deconvolution references for BayesPrism and InstaPrism."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Sequence
import numpy as np
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

if TYPE_CHECKING:
    from .models import DeconvolutionReferenceResult


def export_to_instaprism(
    reference: DeconvolutionReferenceResult,
) -> tuple[np.ndarray, list[str]]:
    """Export reference for InstaPrism deconvolution.

    Returns:
        (reference_matrix, state_labels) where matrix has shape (K_states, G_genes).
    """
    mat, labels, _ = reference.to_numpy(which="state")
    return mat, labels


def export_to_bayesprism(
    reference: DeconvolutionReferenceResult,
    mixture: np.ndarray | None = None,
    bulk_names: Sequence[str] | None = None,
) -> Result[Any, str]:
    """Convert DeconvolutionReferenceResult into BayesPrism S4 equivalent structures.

    If mixture is provided, constructs and returns a `bayesprism.models.Prism` object.
    Otherwise, returns a dictionary containing `RefPhi` instances and `state_to_type_map`.

    Returns:
        Success(Prism) or Success(dict[str, Any]) or Failure(error_message).
    """
    try:
        from bayesprism.models import RefPhi, Prism
    except ImportError:
        return Failure(
            "Package 'bayesprism' is not installed or importable. "
            "Cannot construct BayesPrism models directly."
        )

    try:
        mat_state, states, genes = reference.to_numpy(which="state")
        mat_type, types, _ = reference.to_numpy(which="type")

        ref_phi_state = RefPhi(
            phi=mat_state,
            cell_names=tuple(states),
            gene_names=tuple(genes),
            pseudo_min=1e-8,
        )

        ref_phi_type = RefPhi(
            phi=mat_type,
            cell_names=tuple(types),
            gene_names=tuple(genes),
            pseudo_min=1e-8,
        )

        # Build mapping from cell_type -> tuple[cell_states]
        mapping: dict[str, list[str]] = {t: [] for t in types}
        for row in reference.hierarchy_table.iter_rows(named=True):
            st = str(row["cell_state"])
            tp = str(row["cell_type"])
            if tp in mapping and st not in mapping[tp]:
                mapping[tp].append(st)

        state_to_type_map = {tp: tuple(st_list) for tp, st_list in mapping.items()}

        if mixture is not None:
            n_samples = mixture.shape[0]
            b_names = tuple(bulk_names) if bulk_names else tuple(f"sample_{i}" for i in range(n_samples))
            prism_obj = Prism(
                phi_cell_state=ref_phi_state,
                phi_cell_type=ref_phi_type,
                state_to_type_map=state_to_type_map,
                key=Nothing,
                mixture=mixture,
                bulk_names=b_names,
                gene_names=tuple(genes),
            )
            return Success(prism_obj)

        return Success({
            "phi_cell_state": ref_phi_state,
            "phi_cell_type": ref_phi_type,
            "state_to_type_map": state_to_type_map,
            "malignant_states": reference.malignant_states,
        })
    except Exception as exc:
        return Failure(f"Failed to adapt reference to BayesPrism: {exc}")
