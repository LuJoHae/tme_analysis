from typing import Union, Sequence
import numpy as np
import scipy.sparse as sp
import torch
import polars as pl
from returns.result import Result, Success, Failure
from returns.maybe import Maybe, Some, Nothing

from bayesprism.models import (
    MatrixType,
    RefPhi,
    Prism,
    GibbsControl,
    OptControl,
    GibbsSampler,
    JointPost,
    ThetaPost,
    BayesPrism,
    BayesPrismST,
)
from bayesprism.validation import validate_input
from bayesprism.preprocessing import norm_to_one, collapse, filter_bulk_outlier
from bayesprism.gibbs import run_gibbs
from bayesprism.optimization import update_reference


def new_prism(
    reference: MatrixType,
    cell_type_labels: Sequence[str],
    cell_state_labels: Maybe[Sequence[str]],
    mixture: MatrixType,
    key: Maybe[str] = Nothing,
    gene_names: Maybe[Sequence[str]] = Nothing,
    bulk_names: Maybe[Sequence[str]] = Nothing,
    outlier_cut: float = 0.01,
    outlier_fraction: float = 0.1,
    pseudo_min: float = 1e-8,
) -> Result[Prism, str]:
    """
    Construct a Prism input object from user provided scRNA-seq/GEP reference and bulk mixture data.
    """
    cell_types = tuple(as_str for as_str in cell_type_labels)
    cell_states = cell_state_labels.value_or(cell_types)
    cell_states_tuple = tuple(as_str for as_str in cell_states)

    if len(cell_types) != len(cell_states_tuple):
        return Failure("Error: length of cell_type_labels and cell_state_labels do not match!")

    # Check input validation
    g_names = gene_names.value_or(tuple(f"gene_{i}" for i in range(reference.shape[1])))
    match validate_input(mixture, g_names):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    # Collapse reference by cell states and cell types
    ref_cs_mat, unique_states = collapse(reference, cell_states_tuple)
    ref_ct_mat, unique_types = collapse(reference, cell_types)

    # Filter bulk outliers
    if not isinstance(mixture, torch.Tensor):
        mixture_tensor = torch.from_numpy(np.asarray(mixture)).float()
    else:
        mixture_tensor = mixture.float()

    mixture_filtered, keep_mask = filter_bulk_outlier(
        mixture_tensor, outlier_cut=outlier_cut, outlier_fraction=outlier_fraction
    )

    keep_indices = torch.where(keep_mask)[0]
    ref_cs_filtered = ref_cs_mat[:, keep_indices]
    ref_ct_filtered = ref_ct_mat[:, keep_indices]
    filtered_gene_names = tuple(g_names[i] for i in keep_indices.tolist())

    # Normalize references
    phi_cs = norm_to_one(ref_cs_filtered, pseudo_min=pseudo_min)
    phi_ct = norm_to_one(ref_ct_filtered, pseudo_min=pseudo_min)

    ref_phi_cs = RefPhi(
        phi=phi_cs,
        cell_names=unique_states,
        gene_names=filtered_gene_names,
        pseudo_min=pseudo_min,
    )
    ref_phi_ct = RefPhi(
        phi=phi_ct,
        cell_names=unique_types,
        gene_names=filtered_gene_names,
        pseudo_min=pseudo_min,
    )

    # State to type mapping dictionary
    state_to_type_map: dict[str, tuple[str, ...]] = {}
    for ct in unique_types:
        states_for_ct = tuple(
            dict.fromkeys(
                cell_states_tuple[i]
                for i, t in enumerate(cell_types)
                if t == ct
            )
        )
        state_to_type_map[ct] = states_for_ct

    b_names = bulk_names.value_or(tuple(f"mixture_{n}" for n in range(mixture_filtered.shape[0])))

    prism_obj = Prism(
        phi_cell_state=ref_phi_cs,
        phi_cell_type=ref_phi_ct,
        state_to_type_map=state_to_type_map,
        key=key,
        mixture=mixture_filtered,
        bulk_names=b_names,
        gene_names=filtered_gene_names,
    )

    return Success(prism_obj)


def merge_k(joint_post: JointPost, state_to_type_map: dict[str, tuple[str, ...]]) -> JointPost:
    """Marginalize cell state posterior estimates to cell types."""
    N, G, K_state = joint_post.Z.shape
    cell_types = tuple(state_to_type_map.keys())
    K_type = len(cell_types)

    Z_merged = torch.zeros((N, G, K_type), dtype=torch.float32, device=joint_post.Z.device)
    theta_merged = torch.zeros((N, K_type), dtype=torch.float32, device=joint_post.theta.device)
    theta_cv_merged = torch.zeros((N, K_type), dtype=torch.float32, device=joint_post.theta_cv.device)

    state_names = joint_post.cell_names

    for k, ct in enumerate(cell_types):
        states = state_to_type_map[ct]
        state_indices = [state_names.index(s) for s in states if s in state_names]
        if len(state_indices) > 0:
            Z_merged[:, :, k] = joint_post.Z[:, :, state_indices].sum(dim=2)
            theta_merged[:, k] = joint_post.theta[:, state_indices].sum(dim=1)
            theta_cv_merged[:, k] = joint_post.theta_cv[:, state_indices].mean(dim=1)

    return JointPost(
        Z=Z_merged,
        theta=theta_merged,
        theta_cv=theta_cv_merged,
        constant=joint_post.constant,
        bulk_names=joint_post.bulk_names,
        gene_names=joint_post.gene_names,
        cell_names=cell_types,
    )


def run_prism(
    prism: Prism,
    update_gibbs: bool = True,
    gibbs_control: Maybe[GibbsControl] = Nothing,
    opt_control: Maybe[OptControl] = Nothing,
) -> Result[BayesPrism, str]:
    """
    Main deconvolution workflow for bulk RNA-seq data.
    """
    g_ctrl = gibbs_control.value_or(GibbsControl())
    o_ctrl = opt_control.value_or(OptControl())

    # Step 1: Initial Gibbs sampling over cell states
    sampler_cs = GibbsSampler(
        reference=prism.phi_cell_state,
        X=prism.mixture,
        gibbs_control=g_ctrl,
    )

    match run_gibbs(sampler_cs, final=False):
        case Failure(err):
            return Failure(f"Initial Gibbs sampling failed: {err}")
        case Success(joint_ini_cs):
            if not isinstance(joint_ini_cs, JointPost):
                return Failure("Error: Expected JointPost object from initial Gibbs sampling.")

    # Step 2: Merge states to cell types
    joint_ini_ct = merge_k(joint_ini_cs, prism.state_to_type_map)

    if not update_gibbs:
        bp_obj = BayesPrism(
            prism=prism,
            posterior_initial_cell_state=joint_ini_cs,
            posterior_initial_cell_type=joint_ini_ct,
            gibbs_control=g_ctrl,
            opt_control=o_ctrl,
        )
        return Success(bp_obj)

    # Step 3: Reference refinement & update
    match update_reference(
        Z=joint_ini_ct.Z,
        phi_prime=prism.phi_cell_type,
        state_to_type_map=prism.state_to_type_map,
        key=prism.key,
        opt_control=o_ctrl,
    ):
        case Failure(err):
            return Failure(f"Update reference failed: {err}")
        case Success(psi_ref):
            pass

    # Step 4: Final Gibbs sampling with updated reference
    sampler_update = GibbsSampler(
        reference=psi_ref,
        X=prism.mixture,
        gibbs_control=g_ctrl,
    )

    match run_gibbs(sampler_update, final=True):
        case Failure(err):
            return Failure(f"Final Gibbs sampling failed: {err}")
        case Success(theta_f):
            if not isinstance(theta_f, ThetaPost):
                return Failure("Error: Expected ThetaPost object from final Gibbs sampling.")

    bp_obj = BayesPrism(
        prism=prism,
        posterior_initial_cell_state=joint_ini_cs,
        posterior_initial_cell_type=joint_ini_ct,
        reference_update=Some(psi_ref),
        posterior_theta_f=Some(theta_f),
        gibbs_control=g_ctrl,
        opt_control=o_ctrl,
    )
    return Success(bp_obj)


def run_prism_st(
    prism: Prism,
    gibbs_control: Maybe[GibbsControl] = Nothing,
    opt_control: Maybe[OptControl] = Nothing,
) -> Result[BayesPrismST, str]:
    """
    Deconvolution workflow for Spatial Transcriptomics data using MLE optimizer.
    """
    g_ctrl = gibbs_control.value_or(GibbsControl())
    o_ctrl = opt_control.value_or(OptControl(optimizer="MLE"))

    sampler_ini = GibbsSampler(
        reference=prism.phi_cell_state,
        X=prism.mixture,
        gibbs_control=g_ctrl,
    )

    match run_gibbs(sampler_ini, final=False):
        case Failure(err):
            return Failure(err)
        case Success(joint_ini):
            if not isinstance(joint_ini, JointPost):
                return Failure("Error: Expected JointPost object.")

    match update_reference(
        Z=joint_ini.Z,
        phi_prime=prism.phi_cell_state,
        state_to_type_map=prism.state_to_type_map,
        key=prism.key,
        opt_control=o_ctrl,
    ):
        case Failure(err):
            return Failure(err)
        case Success(psi_ref):
            pass

    sampler_update = GibbsSampler(
        reference=psi_ref,
        X=prism.mixture,
        gibbs_control=g_ctrl,
    )

    match run_gibbs(sampler_update, final=False):
        case Failure(err):
            return Failure(err)
        case Success(joint_update):
            if not isinstance(joint_update, JointPost):
                return Failure("Error: Expected JointPost object.")

    joint_update_ct = merge_k(joint_update, prism.state_to_type_map)

    bp_st = BayesPrismST(
        prism=prism,
        posterior_cell_state=joint_update,
        posterior_cell_type=joint_update_ct,
        reference_update=psi_ref,
        gibbs_control=g_ctrl,
        opt_control=o_ctrl,
    )
    return Success(bp_st)


def get_fraction(
    bp: BayesPrism,
    which_theta: str = "final",
    state_or_type: str = "type",
) -> Result[pl.DataFrame, str]:
    """Extract posterior cell type/state proportion matrix as a Polars DataFrame."""
    match (which_theta, state_or_type):
        case ("first", "state"):
            mat = bp.posterior_initial_cell_state.theta.detach().cpu().numpy()
            cell_names = bp.posterior_initial_cell_state.cell_names
        case ("first", "type"):
            mat = bp.posterior_initial_cell_type.theta.detach().cpu().numpy()
            cell_names = bp.posterior_initial_cell_type.cell_names
        case ("final", _):
            if bp.posterior_theta_f == Nothing:
                return Failure("Error: final theta not present. Run with update_gibbs=True.")
            post_f = bp.posterior_theta_f.unwrap()
            mat = post_f.theta.detach().cpu().numpy()
            cell_names = post_f.cell_names
        case _:
            return Failure(f"Invalid option combination: which_theta={which_theta}, state_or_type={state_or_type}")

    df_dict = {"bulk_id": list(bp.prism.bulk_names)}
    for i, c_name in enumerate(cell_names):
        df_dict[c_name] = list(mat[:, i])

    return Success(pl.DataFrame(df_dict))


def get_exp(
    bp: BayesPrism,
    state_or_type: str,
    cell_name: str,
) -> Result[pl.DataFrame, str]:
    """Extract posterior sample-specific gene expression profile for a target cell state/type."""
    match state_or_type:
        case "state":
            joint_post = bp.posterior_initial_cell_state
        case "type":
            joint_post = bp.posterior_initial_cell_type
        case _:
            return Failure("state_or_type must be either 'state' or 'type'.")

    cell_names = joint_post.cell_names
    if cell_name not in cell_names:
        return Failure(f"Cell name {cell_name} not found in {state_or_type} names.")

    idx = cell_names.index(cell_name)
    mat = joint_post.Z[:, :, idx].detach().cpu().numpy()

    df_dict = {"bulk_id": list(bp.prism.bulk_names)}
    for j, g_name in enumerate(joint_post.gene_names):
        df_dict[g_name] = list(mat[:, j])

    return Success(pl.DataFrame(df_dict))
