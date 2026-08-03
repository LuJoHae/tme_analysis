import pytest
import numpy as np
import torch
from returns.result import Success, Failure
from bayesprism.models import RefPhi, GibbsControl, GibbsSampler, JointPost
from bayesprism.gibbs import run_gibbs, rdirichlet


def test_rdirichlet_properties() -> None:
    alpha = torch.tensor([1.0, 2.0, 3.0], dtype=torch.float32)
    device = torch.device("cpu")
    sample = rdirichlet(alpha, device)

    assert sample.shape == (3,)
    assert torch.all(sample > 0)
    torch.testing.assert_close(sample.sum(), torch.tensor(1.0, dtype=torch.float32))


def test_run_gibbs_execution() -> None:
    phi = torch.tensor([[0.7, 0.3], [0.2, 0.8]], dtype=torch.float32)
    ref = RefPhi(
        phi=phi,
        cell_names=("cell_type_1", "cell_type_2"),
        gene_names=("gene1", "gene2"),
        pseudo_min=1e-8,
    )
    X = torch.tensor([[100, 200], [300, 100]], dtype=torch.float32)
    control = GibbsControl(chain_length=50, burn_in=10, thinning=2, seed=123)

    sampler = GibbsSampler(reference=ref, X=X, gibbs_control=control)
    res = run_gibbs(sampler, final=False)

    match res:
        case Success(joint_post):
            assert isinstance(joint_post, JointPost)
            assert joint_post.theta.shape == (2, 2)
            torch.testing.assert_close(joint_post.theta.sum(dim=1), torch.ones(2, dtype=torch.float32))
        case Failure(err):
            pytest.fail(f"Expected Success, got Failure: {err}")
