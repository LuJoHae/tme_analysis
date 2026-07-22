#!/usr/bin/env python3
import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from firthmodels import FirthLogisticRegression

# Set aesthetic style
sns.set_theme(style="whitegrid")
plt.rcParams.update({"font.size": 11, "figure.autolayout": True})

OUTPUT_DIR = os.path.join("output", "firth_logistic_regression")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def generate_data(n_samples=60, beta_0=0.0, beta_1=5.0, alpha_param=1.0, beta_param=5.0, seed=None):
    """
    Generate synthetic data with continuous X in [0, 1] skewed towards low values
    and binary outcome Y following a logistic function P(Y=1|X) = expit(beta_0 + beta_1 * X).
    """
    if seed is not None:
        np.random.seed(seed)

    # X generated from Beta distribution (low values predominant when alpha=1, beta=5)
    X = np.random.beta(a=alpha_param, b=beta_param, size=n_samples)
    
    # Calculate true probabilities
    p = expit(beta_0 + beta_1 * X)
    
    # Sample binary Y
    Y = np.random.binomial(n=1, p=p)
    
    return X.reshape(-1, 1), Y, p


def fit_both_models(X, Y):
    """
    Fit Firth Logistic Regression and standard MLE Logistic Regression.
    Returns (firth_model, mle_model, (b0_firth, b1_firth), (b0_mle, b1_mle)).
    """
    # Firth model
    firth = FirthLogisticRegression()
    firth.fit(X, Y)
    b0_firth = float(np.ravel(firth.intercept_)[0])
    b1_firth = float(np.ravel(firth.coef_)[0])

    # Standard MLE model (unregularized)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mle = LogisticRegression(penalty=None, solver="lbfgs", max_iter=2000)
        mle.fit(X, Y)
    b0_mle = float(np.ravel(mle.intercept_)[0])
    b1_mle = float(np.ravel(mle.coef_)[0])

    return firth, mle, (b0_firth, b1_firth), (b0_mle, b1_mle)


def plot_single_dataset(beta_0=0.0, beta_1=5.0, seed=42):
    """
    Plot single dataset fit comparing Firth vs MLE vs True Logistic Curve.
    """
    X, Y, _ = generate_data(n_samples=60, beta_0=beta_0, beta_1=beta_1, seed=seed)
    firth, mle, (b0_firth, b1_firth), (b0_mle, b1_mle) = fit_both_models(X, Y)

    x_grid = np.linspace(0, 1, 300).reshape(-1, 1)
    p_true = expit(beta_0 + beta_1 * x_grid)
    p_firth = firth.predict_proba(x_grid)[:, 1]
    p_mle = mle.predict_proba(x_grid)[:, 1]

    fig, ax = plt.subplots(figsize=(9, 6))

    # Add jittered points for visualization
    jitter = np.random.normal(0, 0.02, size=len(Y))
    ax.scatter(X[Y == 0], Y[Y == 0] + jitter[Y == 0], color="#2b5c8f", alpha=0.7, edgecolors="none", label="Observed Y=0")
    ax.scatter(X[Y == 1], Y[Y == 1] + jitter[Y == 1], color="#d95f02", alpha=0.7, edgecolors="none", label="Observed Y=1")

    # Plot probability curves
    ax.plot(x_grid, p_true, "k--", linewidth=2.5, label=f"True Curve (β₀={beta_0:.1f}, β₁={beta_1:.1f})")
    ax.plot(x_grid, p_firth, "#1b9e77", linewidth=2.2, label=f"Firth Fit (β̂₀={b0_firth:.2f}, β̂₁={b1_firth:.2f})")
    ax.plot(x_grid, p_mle, "#e7298a", linewidth=2.2, linestyle="-.", label=f"Standard MLE (β̂₀={b0_mle:.2f}, β̂₁={b1_mle:.2f})")

    ax.set_xlabel("Independent Variable X (Continuous ∈ [0, 1], Skewed Low)")
    ax.set_ylabel("P(Y = 1 | X)")
    ax.set_ylim(-0.08, 1.08)
    ax.set_title(f"Logistic Regression Fit on Skewed Data (N=60)\nTrue: β₀={beta_0}, β₁={beta_1}", fontsize=13, fontweight="bold")
    ax.legend(loc="upper left", frameon=True)

    out_path = os.path.join(OUTPUT_DIR, "single_dataset_fitted_curves.png")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_varying_skewness(beta_0=0.0, beta_1=5.0, seed=123):
    """
    Show datasets under different degrees of X skewness (Uniform vs Beta(1,5) vs Beta(0.5,5)).
    """
    skewness_configs = [
        ("Uniform X ~ U(0,1)", 1.0, 1.0),
        ("Moderate Skew X ~ Beta(1, 5)", 1.0, 5.0),
        ("High Skew X ~ Beta(0.5, 5)", 0.5, 5.0),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

    for idx, (label, a, b) in enumerate(skewness_configs):
        ax = axes[idx]
        X, Y, _ = generate_data(n_samples=70, beta_0=beta_0, beta_1=beta_1, alpha_param=a, beta_param=b, seed=seed + idx)
        firth, mle, (b0_firth, b1_firth), (b0_mle, b1_mle) = fit_both_models(X, Y)

        x_grid = np.linspace(0, 1, 200).reshape(-1, 1)
        p_true = expit(beta_0 + beta_1 * x_grid)
        p_firth = firth.predict_proba(x_grid)[:, 1]
        p_mle = mle.predict_proba(x_grid)[:, 1]

        ax.scatter(X, Y, c=Y, cmap="coolwarm", alpha=0.6, edgecolors="none")
        ax.plot(x_grid, p_true, "k--", label="True")
        ax.plot(x_grid, p_firth, "#1b9e77", label=f"Firth (β₁={b1_firth:.2f})")
        ax.plot(x_grid, p_mle, "#e7298a", linestyle="-.", label=f"MLE (β₁={b1_mle:.2f})")

        ax.set_title(label, fontsize=11, fontweight="bold")
        ax.set_xlabel("X")
        if idx == 0:
            ax.set_ylabel("Y / P(Y=1)")
        ax.legend(loc="upper left", fontsize=9)

    fig.suptitle(f"Effect of Covariate Skewness on Model Fits (N=70, True β₀={beta_0}, β₁={beta_1})", fontsize=14, fontweight="bold")
    out_path = os.path.join(OUTPUT_DIR, "datasets_varying_skewness.png")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


def run_monte_carlo(n_runs=500, n_samples=60, beta_0=0.0, beta_1=5.0, seed=42):
    """
    Monte Carlo simulation to estimate parameter bias for Firth vs MLE.
    """
    np.random.seed(seed)
    
    results = []

    for i in range(n_runs):
        X, Y, _ = generate_data(n_samples=n_samples, beta_0=beta_0, beta_1=beta_1)
        
        # Ensure non-trivial binary outcome
        if len(np.unique(Y)) < 2:
            continue

        try:
            _, _, (b0_firth, b1_firth), (b0_mle, b1_mle) = fit_both_models(X, Y)
            # Clip extreme numerical divergence from complete separation in standard MLE for plotting sanity
            b0_mle_c = np.clip(b0_mle, -30, 30)
            b1_mle_c = np.clip(b1_mle, -30, 30)
            
            results.append({
                "run": i,
                "b0_firth": b0_firth,
                "b1_firth": b1_firth,
                "b0_mle": b0_mle_c,
                "b1_mle": b1_mle_c,
            })
        except Exception:
            continue

    df = pd.DataFrame(results)

    # Compute Bias
    bias_b0_firth = df["b0_firth"].mean() - beta_0
    bias_b1_firth = df["b1_firth"].mean() - beta_1
    bias_b0_mle = df["b0_mle"].mean() - beta_0
    bias_b1_mle = df["b1_mle"].mean() - beta_1

    print("=== Monte Carlo Simulation Results (N_samples={}, N_runs={}) ===".format(n_samples, len(df)))
    print(f"Intercept β₀ (True = {beta_0}):")
    print(f"  Firth Mean: {df['b0_firth'].mean():.3f} | Bias: {bias_b0_firth:+.3f}")
    print(f"  MLE   Mean: {df['b0_mle'].mean():.3f} | Bias: {bias_b0_mle:+.3f}")
    print(f"Slope β₁ (True = {beta_1}):")
    print(f"  Firth Mean: {df['b1_firth'].mean():.3f} | Bias: {bias_b1_firth:+.3f}")
    print(f"  MLE   Mean: {df['b1_mle'].mean():.3f} | Bias: {bias_b1_mle:+.3f}")

    # Plot Parameter Distributions
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Intercept beta_0
    sns.kdeplot(df["b0_firth"], ax=axes[0], color="#1b9e77", fill=True, alpha=0.3, label=f"Firth (Bias={bias_b0_firth:+.2f})")
    sns.kdeplot(df["b0_mle"], ax=axes[0], color="#e7298a", fill=True, alpha=0.3, label=f"MLE (Bias={bias_b0_mle:+.2f})")
    axes[0].axvline(beta_0, color="black", linestyle="--", linewidth=2, label=f"True β₀ ({beta_0})")
    axes[0].set_title("Intercept Parameter Estimate (β̂₀)", fontsize=12, fontweight="bold")
    axes[0].set_xlabel("Estimated Intercept β̂₀")
    axes[0].legend(loc="upper left")

    # Slope beta_1
    sns.kdeplot(df["b1_firth"], ax=axes[1], color="#1b9e77", fill=True, alpha=0.3, label=f"Firth (Bias={bias_b1_firth:+.2f})")
    sns.kdeplot(df["b1_mle"], ax=axes[1], color="#e7298a", fill=True, alpha=0.3, label=f"MLE (Bias={bias_b1_mle:+.2f})")
    axes[1].axvline(beta_1, color="black", linestyle="--", linewidth=2, label=f"True β₁ ({beta_1})")
    axes[1].set_title("Slope Parameter Estimate (β̂₁)", fontsize=12, fontweight="bold")
    axes[1].set_xlabel("Estimated Slope β̂₁")
    axes[1].legend(loc="upper right")

    fig.suptitle(f"Monte Carlo Parameter Estimate Distributions ({len(df)} runs, N={n_samples})", fontsize=14, fontweight="bold")
    out_path = os.path.join(OUTPUT_DIR, "parameter_distributions_monte_carlo.png")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")

    return df


def run_sample_size_bias_trend(sample_sizes=[15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100], n_runs=300, beta_0=0.0, beta_1=5.0, seed=101):
    """
    Analyze how parameter estimation bias and standard deviation change with sample size for Firth vs MLE (range 10 to 100).
    """
    np.random.seed(seed)
    
    bias_records = []

    for n in sample_sizes:
        b0_firth_list, b1_firth_list = [], []
        b0_mle_list, b1_mle_list = [], []

        for _ in range(n_runs):
            X, Y, _ = generate_data(n_samples=n, beta_0=beta_0, beta_1=beta_1)
            if len(np.unique(Y)) < 2:
                continue

            try:
                _, _, (b0_f, b1_f), (b0_m, b1_m) = fit_both_models(X, Y)
                b0_firth_list.append(b0_f)
                b1_firth_list.append(b1_f)
                b0_mle_list.append(np.clip(b0_m, -40, 40))
                b1_mle_list.append(np.clip(b1_m, -40, 40))
            except Exception:
                continue

        bias_records.append({
            "Sample Size N": n,
            "Firth Intercept Bias": np.mean(b0_firth_list) - beta_0,
            "Firth Intercept SD": np.std(b0_firth_list),
            "MLE Intercept Bias": np.mean(b0_mle_list) - beta_0,
            "MLE Intercept SD": np.std(b0_mle_list),
            "Firth Slope Bias": np.mean(b1_firth_list) - beta_1,
            "Firth Slope SD": np.std(b1_firth_list),
            "MLE Slope Bias": np.mean(b1_mle_list) - beta_1,
            "MLE Slope SD": np.std(b1_mle_list),
        })

    bias_df = pd.DataFrame(bias_records)

    # Plot bias vs sample size with shaded standard deviation bands (±1 SD)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Intercept Bias & SD
    axes[0].plot(bias_df["Sample Size N"], bias_df["Firth Intercept Bias"], "o-", color="#1b9e77", linewidth=2.5, label="Firth Mean Bias")
    axes[0].fill_between(
        bias_df["Sample Size N"],
        bias_df["Firth Intercept Bias"] - bias_df["Firth Intercept SD"],
        bias_df["Firth Intercept Bias"] + bias_df["Firth Intercept SD"],
        color="#1b9e77", alpha=0.15, label="Firth ±1 SD"
    )
    axes[0].plot(bias_df["Sample Size N"], bias_df["MLE Intercept Bias"], "s--", color="#e7298a", linewidth=2.5, label="Standard MLE Mean Bias")
    axes[0].fill_between(
        bias_df["Sample Size N"],
        bias_df["MLE Intercept Bias"] - bias_df["MLE Intercept SD"],
        bias_df["MLE Intercept Bias"] + bias_df["MLE Intercept SD"],
        color="#e7298a", alpha=0.15, label="Standard MLE ±1 SD"
    )
    axes[0].axhline(0, color="black", linestyle="--", linewidth=1.5)
    axes[0].set_title("Intercept Bias (β̂₀ - β₀) with ±1 SD Band", fontsize=12, fontweight="bold")
    axes[0].set_xlabel("Sample Size N")
    axes[0].set_ylabel("Parameter Estimation Bias")
    axes[0].legend(loc="best")

    # Slope Bias & SD
    axes[1].plot(bias_df["Sample Size N"], bias_df["Firth Slope Bias"], "o-", color="#1b9e77", linewidth=2.5, label="Firth Mean Bias")
    axes[1].fill_between(
        bias_df["Sample Size N"],
        bias_df["Firth Slope Bias"] - bias_df["Firth Slope SD"],
        bias_df["Firth Slope Bias"] + bias_df["Firth Slope SD"],
        color="#1b9e77", alpha=0.15, label="Firth ±1 SD"
    )
    axes[1].plot(bias_df["Sample Size N"], bias_df["MLE Slope Bias"], "s--", color="#e7298a", linewidth=2.5, label="Standard MLE Mean Bias")
    axes[1].fill_between(
        bias_df["Sample Size N"],
        bias_df["MLE Slope Bias"] - bias_df["MLE Slope SD"],
        bias_df["MLE Slope Bias"] + bias_df["MLE Slope SD"],
        color="#e7298a", alpha=0.15, label="Standard MLE ±1 SD"
    )
    axes[1].axhline(0, color="black", linestyle="--", linewidth=1.5)
    axes[1].set_title("Slope Bias (β̂₁ - β₁) with ±1 SD Band", fontsize=12, fontweight="bold")
    axes[1].set_xlabel("Sample Size N")
    axes[1].set_ylabel("Parameter Estimation Bias")
    axes[1].legend(loc="best")

    fig.suptitle("Parameter Estimation Bias & Standard Deviation vs. Sample Size (Skewed Covariate X)", fontsize=14, fontweight="bold")
    out_path = os.path.join(OUTPUT_DIR, "bias_vs_sample_size.png")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


def main():
    print("Starting Firth Logistic Regression vs Standard MLE comparison...")
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Step 1: Single dataset fit & probability curves
    plot_single_dataset(beta_0=0.0, beta_1=5.0, seed=42)

    # Step 2: Varying degrees of covariate skewness
    plot_varying_skewness(beta_0=0.0, beta_1=5.0, seed=123)

    # Step 3: Monte Carlo simulation for bias comparison
    run_monte_carlo(n_runs=500, n_samples=60, beta_0=0.0, beta_1=5.0, seed=42)

    # Step 4: Bias vs Sample Size
    run_sample_size_bias_trend(sample_sizes=[15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100], n_runs=300, beta_0=0.0, beta_1=5.0, seed=101)

    print("All tasks completed successfully!")


if __name__ == "__main__":
    main()
