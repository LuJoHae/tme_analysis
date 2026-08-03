from typing import Generic, TypeVar, Protocol, Any, Union
import numpy as np
import scipy.sparse as sp
import torch
from pydantic import BaseModel, ConfigDict, Field, field_validator
from returns.maybe import Maybe, Some, Nothing

# Type Aliases for Tensors / Matrices
MatrixType = Union[torch.Tensor, sp.csr_matrix, np.ndarray]

class Reference(Protocol):
    """Protocol for reference objects."""
    @property
    def pseudo_min(self) -> float: ...


class RefPhi(BaseModel):
    """S4 class equivalent refPhi for non-malignant cell types/states."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    phi: MatrixType
    cell_names: tuple[str, ...]
    gene_names: tuple[str, ...]
    pseudo_min: float


class RefTumor(BaseModel):
    """S4 class equivalent refTumor for tumor datasets."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    psi_mal: MatrixType
    psi_env: MatrixType
    key: str
    bulk_names: tuple[str, ...]
    env_cell_names: tuple[str, ...]
    gene_names: tuple[str, ...]
    pseudo_min: float


class Prism(BaseModel):
    """S4 class equivalent prism containing inputs and mappings."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    phi_cell_state: RefPhi
    phi_cell_type: RefPhi
    state_to_type_map: dict[str, tuple[str, ...]]
    key: Maybe[str]
    mixture: MatrixType
    bulk_names: tuple[str, ...]
    gene_names: tuple[str, ...]


class GibbsControl(BaseModel):
    """Controls for Dirichlet-Multinomial Gibbs sampler."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    chain_length: int = 1000
    burn_in: int = 500
    thinning: int = 2
    seed: Maybe[int] = Some(123)
    alpha: float = 1.0
    device: str = "cpu"
    num_threads: int = 1

    @field_validator("seed", mode="before")
    @classmethod
    def validate_seed(cls, v: Any) -> Maybe[int]:
        if isinstance(v, int):
            return Some(v)
        return v


class OptControl(BaseModel):
    """Controls for conjugate gradient / L-BFGS optimization of gamma."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    maxit: int = 100000
    optimizer: str = "MAP"  # "MAP" or "MLE"
    sigma: float = 2.0
    eps: float = 1e-7
    device: str = "cpu"


class NmfControl(BaseModel):
    """Controls for NMF initialization in tumor embedding learning."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    nrun: int = 200
    seed: int = 123


class GibbsSampler(BaseModel):
    """Sampler state holding input parameters."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    reference: Union[RefPhi, RefTumor]
    X: MatrixType
    gibbs_control: GibbsControl


class JointPost(BaseModel):
    """Posterior expectations for Z and theta from Gibbs sampling."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    Z: torch.Tensor  # Shape (N, G, K)
    theta: torch.Tensor  # Shape (N, K)
    theta_cv: torch.Tensor  # Shape (N, K)
    constant: float
    bulk_names: tuple[str, ...]
    gene_names: tuple[str, ...]
    cell_names: tuple[str, ...]


class ThetaPost(BaseModel):
    """Posterior expectations for theta from updated Gibbs sampling."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    theta: torch.Tensor  # Shape (N, K)
    theta_cv: torch.Tensor  # Shape (N, K)
    bulk_names: tuple[str, ...]
    cell_names: tuple[str, ...]


class BayesPrism(BaseModel):
    """S4 class equivalent BayesPrism output for bulk RNA-seq."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    prism: Prism
    posterior_initial_cell_state: JointPost
    posterior_initial_cell_type: JointPost
    reference_update: Maybe[Union[RefPhi, RefTumor]] = Field(default=Nothing)
    posterior_theta_f: Maybe[ThetaPost] = Field(default=Nothing)
    gibbs_control: GibbsControl
    opt_control: OptControl


class BayesPrismST(BaseModel):
    """S4 class equivalent BayesPrismST output for spatial transcriptomics."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    prism: Prism
    posterior_cell_state: JointPost
    posterior_cell_type: JointPost
    reference_update: Union[RefPhi, RefTumor]
    gibbs_control: GibbsControl
    opt_control: OptControl
