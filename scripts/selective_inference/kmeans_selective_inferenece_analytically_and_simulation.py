import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import selective_inference as si

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def main():
    parser = argparse.ArgumentParser(description="K-Means Selective Inference simulation (Analytical & MCMC)")
    parser.add_argument("--n-simulations", type=int, default=50, help="Number of Monte Carlo simulations to run")
    parser.add_argument("--mcmc-steps", type=int, default=1500, help="Number of MCMC steps")
    parser.add_argument("--burn-in", type=int, default=300, help="MCMC burn-in steps")
    parser.add_argument("--mean-diff", type=float, default=1.5, help="True mean difference (0 for Null)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()
    
    # Styling plots
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'figure.titlesize': 18
    })
    
    out_dir = "output/kmeans_selective_inference_analytical_and_sim_plots"
    os.makedirs(out_dir, exist_ok=True)
    logging.info(f"Output directory established: {out_dir}")
    
    si.run_single_run_diagnostics(args, out_dir)

    si.run_ecdf_simulation(args, out_dir)    
    
    print("\n=======================================================")
    print("K-Means Selective Inference Simulation Successfully Completed!")
    print(f"All plots and data saved to: {out_dir}/")
    print("=======================================================")

if __name__ == "__main__":
    main()
