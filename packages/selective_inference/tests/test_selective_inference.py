import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from selective_inference import generate_mixture_data, run_mcmc_metropolis_hastings, run_single_run_diagnostics

class DummyArgs:
    def __init__(self, mean_diff=1.5, seed=42, mcmc_steps=20, burn_in=5, model_type="kmeans"):
        self.mean_diff = mean_diff
        self.seed = seed
        self.mcmc_steps = mcmc_steps
        self.burn_in = burn_in
        self.model_type = model_type

def test_generate_mixture_data():
    X, true_labels = generate_mixture_data(
        n_samples_comp1=10, n_samples_comp2=10, n_features=2, 
        mean_diff=1.0, sigma=1.0, random_state=42
    )
    assert X.shape == (20, 2)
    assert len(true_labels) == 20
    assert np.sum(true_labels == 0) == 10
    assert np.sum(true_labels == 1) == 10

def test_run_mcmc_metropolis_hastings_kmeans():
    X, true_labels = generate_mixture_data(
        n_samples_comp1=20, n_samples_comp2=20, n_features=2, 
        mean_diff=0.0, sigma=1.0, random_state=42
    )
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=5)
    kmeans.fit(X)
    pred_labels = kmeans.labels_
    
    mcmc_p, accepted_z, acceptance_rate = run_mcmc_metropolis_hastings(
        X, pred_labels, c1=0, c2=1, sigma=1.0, max_mcmc_steps=50, 
        burn_in=10, random_state=42, model_type="kmeans"
    )
    assert 0.0 <= mcmc_p <= 1.0
    assert len(accepted_z) == 40
    assert 0.0 <= acceptance_rate <= 1.0

def test_run_mcmc_metropolis_hastings_gmm_spherical():
    X, true_labels = generate_mixture_data(
        n_samples_comp1=20, n_samples_comp2=20, n_features=2, 
        mean_diff=0.0, sigma=1.0, random_state=42
    )
    gmm = GaussianMixture(n_components=2, covariance_type="spherical", random_state=42)
    gmm.fit(X)
    pred_labels = gmm.predict(X)
    
    mcmc_p, accepted_z, acceptance_rate = run_mcmc_metropolis_hastings(
        X, pred_labels, c1=0, c2=1, sigma=1.0, max_mcmc_steps=50, 
        burn_in=10, random_state=42, model_type="gmm_spherical"
    )
    assert 0.0 <= mcmc_p <= 1.0
    assert len(accepted_z) == 40
    assert 0.0 <= acceptance_rate <= 1.0

def test_run_single_run_diagnostics_kmeans():
    args = DummyArgs(model_type="kmeans")
    chart, trace, density = run_single_run_diagnostics(args)
    assert chart is not None
    assert trace is not None
    assert density is not None
    assert "K-Means" in chart.title

def test_run_single_run_diagnostics_gmm_spherical():
    args = DummyArgs(model_type="gmm_spherical")
    chart, trace, density = run_single_run_diagnostics(args)
    assert chart is not None
    assert trace is not None
    assert density is not None
    assert "GMM" in chart.title
