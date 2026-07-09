#!/usr/bin/env python3
"""
Calculate and plot the most frequent mutations for each iAtlas cohort using datalair.
"""

import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Ensure local packages are in python path if running directly
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))

import datalair
from ml_pipelines import IAtlasMostFrequentMutations


def plot_mutation_frequencies(df_top, ds_name, total_samples, output_path):
    """
    Generate a premium horizontal bar chart showing the most frequent mutations in the cohort.
    """
    fig, ax = plt.subplots(figsize=(10, 6.5), dpi=300)
    
    # Use a premium sequential color palette
    colors = sns.color_palette("viridis", len(df_top))
    
    # Horizontal bar plot
    bars = ax.barh(df_top['Gene'], df_top['Frequency_Pct'], color=colors, edgecolor='none', height=0.6)
    
    # Invert y-axis to have the most frequent at the top
    ax.invert_yaxis()
    
    # Customize grid and spines
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
    ax.yaxis.grid(False)
    
    for spine in ['top', 'right', 'left', 'bottom']:
        ax.spines[spine].set_visible(False)
        
    # Add labels and title
    ax.set_xlabel('Mutation Frequency (%)', fontsize=12, fontweight='bold', labelpad=10, color='#2C3E50')
    ax.set_ylabel('Gene (Hugo Symbol)', fontsize=12, fontweight='bold', labelpad=10, color='#2C3E50')
    ax.set_title(f'Most Frequent Mutations in {ds_name}\n(Cohort Size n = {total_samples})', 
                 fontsize=14, fontweight='bold', pad=15, color='#2C3E50')
                 
    # Adjust tick labels
    ax.tick_params(axis='both', which='major', labelsize=10, length=0)
    for label in ax.get_yticklabels():
        label.set_fontweight('semibold')
        label.set_color('#34495E')
        
    # Add text labels on the bars
    for bar, (_, row) in zip(bars, df_top.iterrows()):
        width = bar.get_width()
        mut_count = int(row['Mutated_Samples'])
        pct = row['Frequency_Pct']
        label_text = f" {pct:.1f}% ({mut_count}/{total_samples})"
        # Draw text at the end of the bar
        ax.text(width, bar.get_y() + bar.get_height()/2, label_text,
                va='center', ha='left', fontsize=9.5, fontweight='bold', color='#2C3E50')
                
    # Set X-axis limit with some padding for the labels
    max_freq = df_top['Frequency_Pct'].max()
    ax.set_xlim(0, min(100, max_freq + 12))
    
    plt.tight_layout()
    plt.savefig(output_path, format='svg', bbox_inches='tight')
    plt.close(fig)


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    # Resolve project root and output directories
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output"
    mut_output_dir = output_dir / "iAtlas-mutations"
    mut_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Derive the mutation frequencies using the new datalair Dataset
    print("Deriving mutation frequencies dataset using datalair...")
    ds_mut = IAtlasMostFrequentMutations()
    lair.safe_derive(ds_mut)
    
    filepaths = lair.get_dataset_filepaths(ds_mut)
    df_summary = pd.read_csv(filepaths["most_frequent_mutations.csv"])
    
    # Copy/Save the CSV file to the dedicated output directory
    summary_path = mut_output_dir / "most_frequent_mutations.csv"
    df_summary.to_csv(summary_path, index=False)
    print(f"Saved complete mutation frequencies table to {summary_path}")
    
    # Plotting
    cohorts = df_summary['Cohort'].unique()
    print(f"Found {len(cohorts)} cohorts with mutation data. Plotting top mutations...")
    
    for cohort in cohorts:
        cohort_df = df_summary[df_summary['Cohort'] == cohort]
        total_samples = cohort_df.iloc[0]['Total_Samples']
        
        # Take top 15 most frequent mutated genes
        df_top = cohort_df.head(15).copy()
        df_top = df_top.rename(columns={
            'Hugo_Symbol': 'Gene',
            'Mutation_Frequency_Pct': 'Frequency_Pct'
        })
        
        plot_path = mut_output_dir / f"{cohort}_mutations.svg"
        print(f"  -> Plotting top mutations to {plot_path.name}...")
        plot_mutation_frequencies(df_top, cohort, total_samples, plot_path)
        
    # Display top 5 genes for each cohort in the console output
    print("\nSummary of Top 5 Mutated Genes per Cohort:")
    for cohort in cohorts:
        print(f"\n--- {cohort} ---")
        cohort_df = df_summary[df_summary['Cohort'] == cohort].head(5)
        for _, row in cohort_df.iterrows():
            print(f"  {row['Hugo_Symbol']}: {row['Mutation_Frequency_Pct']:.1f}% ({int(row['Mutated_Samples'])}/{int(row['Total_Samples'])})")


if __name__ == "__main__":
    main()
