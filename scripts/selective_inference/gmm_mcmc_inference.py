import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import logging

class GMMSelectiveInference:
    def __init__(self, n_components=2, random_state=42, max_mcmc_steps=2000, time_limit=None):
        self.n_components = n_components
        self.random_state = random_state
        self.max_mcmc_steps = max_mcmc_steps
        self.gmm = GaussianMixture(n_components=n_components, random_state=random_state)
        self.labels_ = None
        self.time_limit = time_limit # To be handled by caller ideally, but good to store
        
    def fit(self, X):
        self.gmm.fit(X)
        self.labels_ = self.gmm.predict(X)
        return self
        
    def test_difference_in_means(self, X, c1, c2, sigma=1.0, burn_in=500):
        n = X.shape[0]
        d = X.shape[1]
        
        c1_mask = self.labels_ == c1
        c2_mask = self.labels_ == c2
        
        n1 = np.sum(c1_mask)
        n2 = np.sum(c2_mask)
        
        if n1 == 0 or n2 == 0:
            return 1.0, 1.0, {} # Fallback
            
        mean1 = np.mean(X[c1_mask], axis=0)
        mean2 = np.mean(X[c2_mask], axis=0)
        
        diff = mean1 - mean2
        diff_norm = np.linalg.norm(diff)
        
        if diff_norm == 0:
            return 1.0, 1.0, {}
            
        dir_vec = diff / diff_norm
        
        # Construct V matrix
        V = np.zeros_like(X)
        V[c1_mask] = dir_vec / n1
        V[c2_mask] = -dir_vec / n2
        
        phi_obs = np.sum(V * X) # Equivalent to trace(V^T X)
        v_norm_sq = np.sum(V**2)
        std_phi = sigma * np.sqrt(v_norm_sq)
        
        # Naive p-value
        z_val_naive = phi_obs / std_phi
        naive_p = 2 * (1 - norm.cdf(abs(z_val_naive)))
        
        # MCMC Hit-and-Run on the 1D line Z
        # Target distribution: z ~ N(0, std_phi^2) restricted to C(X(z)) == c_obs
        
        accepted_z = []
        z_current = phi_obs
        
        # Proposal standard deviation (tuned for typical acceptance rate)
        tau = std_phi * 0.5 
        
        accept_count = 0
        total_proposals = 0
        
        for step in range(self.max_mcmc_steps):
            z_prop = np.random.normal(loc=z_current, scale=tau)
            
            # Construct X(z_prop)
            X_prop = X + (z_prop - phi_obs) * V
            
            # Run GMM on X_prop
            gmm_prop = GaussianMixture(n_components=self.n_components, random_state=self.random_state)
            gmm_prop.fit(X_prop)
            labels_prop = gmm_prop.predict(X_prop)
            
            # Check if assignments match perfectly
            if np.array_equal(labels_prop, self.labels_):
                # Metropolis-Hastings acceptance ratio
                # Since proposal is symmetric, ratio is just target(z_prop) / target(z_current)
                log_ratio = -0.5 * ((z_prop/std_phi)**2 - (z_current/std_phi)**2)
                if np.log(np.random.uniform()) < log_ratio:
                    z_current = z_prop
                    accept_count += 1
            
            total_proposals += 1
            if step >= burn_in:
                accepted_z.append(z_current)
                
        accepted_z = np.array(accepted_z)
        
        if len(accepted_z) == 0:
            selective_p = 1.0
        else:
            # Empirical p-value: fraction of samples more extreme than observed
            selective_p = 2 * min(np.mean(accepted_z >= phi_obs), np.mean(accepted_z <= phi_obs))
            
        plot_data = {
            'accepted_z': accepted_z,
            'phi_obs': phi_obs,
            'std_phi': std_phi,
            'acceptance_rate': accept_count / total_proposals
        }
        
        return naive_p, selective_p, plot_data
