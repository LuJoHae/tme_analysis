import numpy as np
from scipy.stats import norm

class KMeansSelectiveInference:
    def __init__(self, n_clusters, max_iter=100, random_state=None):
        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.random_state = random_state
        self.cluster_centers_ = None
        self.labels_ = None
        self.assignment_path_ = []
        
    def fit(self, X):
        if self.random_state is not None:
            np.random.seed(self.random_state)
        n_samples, n_features = X.shape
        
        # Initialize centers randomly from data points
        initial_indices = np.random.choice(n_samples, self.n_clusters, replace=False)
        centers = X[initial_indices].copy()
        
        self.assignment_path_ = []
        labels = np.zeros(n_samples, dtype=int)
        
        for iteration in range(self.max_iter):
            # Compute distances
            # distances shape: (n_samples, n_clusters)
            distances = np.linalg.norm(X[:, np.newaxis] - centers, axis=2)**2
            new_labels = np.argmin(distances, axis=1)
            
            self.assignment_path_.append(new_labels.copy())
            
            if np.array_equal(labels, new_labels):
                break
                
            labels = new_labels
            
            # Update centers
            for k in range(self.n_clusters):
                if np.sum(labels == k) > 0:
                    centers[k] = np.mean(X[labels == k], axis=0)
                    
        self.cluster_centers_ = centers
        self.labels_ = labels
        return self

    def _compute_v_min_max(self, X, c1, c2):
        """
        Compute the truncation limits [V_min, V_max] for the contrast testing
        difference in means between cluster c1 and c2.
        We condition only on the final cluster assignments to approximate the polyhedral region,
        which bounds the contrast phi.
        """
        n_samples, n_features = X.shape
        labels = self.labels_
        
        # Contrast matrix nu (n_samples, n_features)
        # We test the difference in means vector v_diff
        v_diff = self.cluster_centers_[c1] - self.cluster_centers_[c2]
        
        n1 = np.sum(labels == c1)
        n2 = np.sum(labels == c2)
        
        nu = np.zeros_like(X)
        if n1 > 0:
            nu[labels == c1] = v_diff / n1
        if n2 > 0:
            nu[labels == c2] = -v_diff / n2
        
        # Norm of contrast
        nu_norm_sq = np.sum(nu**2)
        if nu_norm_sq == 0:
            return -np.inf, np.inf, 0.0, 0.0
            
        # phi is the observed test statistic
        phi_obs = np.sum(nu * X)
        
        # X_perp is the projection orthogonal to nu
        X_perp = X - nu * (phi_obs / nu_norm_sq)
        
        V_min = -np.inf
        V_max = np.inf
        
        # For each point i in cluster k, it must be closer to center k than to any center j
        # || x_i - mu_k ||^2 <= || x_i - mu_j ||^2
        # where x_i(phi) = X_perp_i + nu_i * phi
        # mu_k(phi) = mean_{l in C_k} (X_perp_l + nu_l * phi)
        
        for iteration_labels in [self.labels_]: # Just use final labels for approximation to save compute
            # Compute centers as functions of phi
            # center_k(phi) = center_k_perp + center_k_nu * phi
            centers_perp = np.zeros((self.n_clusters, n_features))
            centers_nu = np.zeros((self.n_clusters, n_features))
            
            for k in range(self.n_clusters):
                mask = iteration_labels == k
                if np.sum(mask) > 0:
                    centers_perp[k] = np.mean(X_perp[mask], axis=0)
                    centers_nu[k] = np.mean(nu[mask], axis=0)
                    
            for i in range(n_samples):
                k = iteration_labels[i]
                
                x_i_perp = X_perp[i]
                x_i_nu = nu[i]
                
                for j in range(self.n_clusters):
                    if k == j:
                        continue
                        
                    # Condition: || x_i(phi) - mu_k(phi) ||^2 <= || x_i(phi) - mu_j(phi) ||^2
                    # Let A(phi) = x_i(phi) - mu_k(phi) = (x_i_perp - mu_k_perp) + phi * (x_i_nu - mu_k_nu)
                    # Let B(phi) = x_i(phi) - mu_j(phi) = (x_i_perp - mu_j_perp) + phi * (x_i_nu - mu_j_nu)
                    # Condition: ||A(phi)||^2 - ||B(phi)||^2 <= 0
                    
                    a_perp = x_i_perp - centers_perp[k]
                    a_nu = x_i_nu - centers_nu[k]
                    
                    b_perp = x_i_perp - centers_perp[j]
                    b_nu = x_i_nu - centers_nu[j]
                    
                    # ||A(phi)||^2 = ||a_perp||^2 + 2 phi <a_perp, a_nu> + phi^2 ||a_nu||^2
                    # We want: c2 * phi^2 + c1 * phi + c0 <= 0
                    
                    c2 = np.sum(a_nu**2) - np.sum(b_nu**2)
                    c1 = 2 * (np.sum(a_perp * a_nu) - np.sum(b_perp * b_nu))
                    c0 = np.sum(a_perp**2) - np.sum(b_perp**2)
                    
                    # Solve quadratic inequality c2 * phi^2 + c1 * phi + c0 <= 0
                    # For the observed phi, this should be <= 0
                    
                    if np.abs(c2) < 1e-10:
                        # Linear inequality c1 * phi + c0 <= 0
                        if c1 > 1e-10:
                            # phi <= -c0 / c1
                            V_max = min(V_max, -c0 / c1)
                        elif c1 < -1e-10:
                            # phi >= -c0 / c1
                            V_min = max(V_min, -c0 / c1)
                    else:
                        # Roots of c2 * phi^2 + c1 * phi + c0 = 0
                        discriminant = c1**2 - 4 * c2 * c0
                        if discriminant >= 0:
                            root1 = (-c1 - np.sqrt(discriminant)) / (2 * c2)
                            root2 = (-c1 + np.sqrt(discriminant)) / (2 * c2)
                            
                            r_min = min(root1, root2)
                            r_max = max(root1, root2)
                            
                            if c2 > 0:
                                # c2 > 0 means parabola opens upwards, <= 0 between roots
                                V_min = max(V_min, r_min)
                                V_max = min(V_max, r_max)
                            else:
                                # c2 < 0 means parabola opens downwards, <= 0 outside roots
                                # This creates a union of intervals. We keep the interval containing phi_obs
                                if phi_obs <= r_min:
                                    V_max = min(V_max, r_min)
                                elif phi_obs >= r_max:
                                    V_min = max(V_min, r_max)
                                    
        return V_min, V_max, phi_obs, np.sqrt(nu_norm_sq)

    def test_difference_in_means(self, X, c1, c2, sigma=1.0):
        """
        Test the difference in means between cluster c1 and c2.
        sigma is the noise variance of the data. 
        """
        V_min, V_max, phi_obs, nu_norm = self._compute_v_min_max(X, c1, c2)
        
        if nu_norm == 0:
            return 1.0, 1.0, {'V_min': -np.inf, 'V_max': np.inf, 'phi_obs': 0, 'std_phi': 1.0}
            
        # The variance of phi is sigma^2 * ||nu||^2 / ||nu||^4 = sigma^2 / ||nu||^2
        # wait, phi = sum(nu * X) / ||nu||^2. Variance is sum(nu^2) * sigma^2 / ||nu||^4 = sigma^2 / ||nu||^2
        variance_phi = (sigma**2) / (nu_norm**2)
        std_phi = np.sqrt(variance_phi)
        
        # Calculate naive p-value (assuming normally distributed without truncation)
        z_obs = np.abs(phi_obs) / std_phi
        naive_p = 2 * (1 - norm.cdf(z_obs))
        
        # Calculate selective p-value using truncated normal
        # Survival function of truncated normal
        # P(|Z| > z_obs | Z in [V_min/std, V_max/std])
        
        alpha_min = V_min / std_phi
        alpha_max = V_max / std_phi
        z_val = phi_obs / std_phi
        
        # CDF of truncated normal at z
        denom = norm.cdf(alpha_max) - norm.cdf(alpha_min)
        if denom <= 0:
            return naive_p, 1.0, {'V_min': V_min, 'V_max': V_max, 'phi_obs': phi_obs, 'std_phi': std_phi} # Fallback
            
        cdf_z = (norm.cdf(z_val) - norm.cdf(alpha_min)) / denom
        
        # 2-sided selective p-value
        selective_p = 2 * min(cdf_z, 1 - cdf_z)
        
        plot_data = {
            'V_min': V_min,
            'V_max': V_max,
            'phi_obs': phi_obs,
            'std_phi': std_phi
        }
        
        return naive_p, selective_p, plot_data
