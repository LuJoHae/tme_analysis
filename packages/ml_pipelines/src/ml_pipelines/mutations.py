import datalair
import numpy as np
import pandas as pd
import ici_datasets
from pathlib import Path

class IAtlasMostFrequentMutations(datalair.Dataset):
    """
    Datalair Dataset containing the calculated mutation frequencies across all iAtlas cohorts.
    """
    def __init__(self) -> None:
        super().__init__(namespace="IAtlasMostFrequentMutations", dataset_name="most_frequent_mutations")

    def derive(self, lair: datalair.Lair) -> None:
        df_summary = calculate_iatlas_mutation_frequencies(lair)
        output_dir = lair.get_path(self)
        output_dir.mkdir(parents=True, exist_ok=True)
        df_summary.to_csv(output_dir / "most_frequent_mutations.csv", index=False)


def calculate_iatlas_mutation_frequencies(lair: datalair.Lair) -> pd.DataFrame:
    """
    Calculate somatic mutation frequencies for each iAtlas cohort.
    
    Filters mutations using the standard TMB pipeline criteria:
    - VAF >= 0.05
    - Variant_Classification is protein-altering/coding.
    """
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    iatlas_dataset_names = [s for s in dataset_class.get_dataset_name() if "iAtlas" in s]
    
    all_summary_data = []
    
    for ds_name in iatlas_dataset_names:
        try:
            data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, dataset_class, ds_name)
            mut_file = data_dir / "data_mutations.txt"
            
            if not mut_file.exists():
                continue
                
            # Load clinical samples using the existing function to establish cohort size
            df_clinical = ici_datasets.cbioportal_datasets.load_data_clinical(data_dir)
            total_samples = len(df_clinical)
            
            if total_samples == 0:
                continue
                
            # Read mutation data
            df_mut = pd.read_csv(mut_file, sep="\t", low_memory=False)
            
            # Compute Variant Allele Frequency (VAF)
            df_mut['VAF'] = np.where(df_mut['t_depth'] > 0, df_mut['t_alt_count'] / df_mut['t_depth'], 0)
            
            # Apply TMB filters: coding classifications and VAF >= 0.05
            qualifying_classifications = {
                'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation',
                'Frame_Shift_Ins', 'Frame_Shift_Del', 'In_Frame_Ins',
                'In_Frame_Del', 'Translation_Start_Site', 'Splice_Site'
            }
            
            qualifying_subset = df_mut[
                (df_mut['Variant_Classification'].isin(qualifying_classifications)) &
                (df_mut['VAF'] >= 0.05)
            ].copy()
            
            # Restrict to mutations that belong to clinical samples
            qualifying_subset = qualifying_subset[qualifying_subset['Tumor_Sample_Barcode'].isin(df_clinical.index)]
            
            if len(qualifying_subset) == 0:
                continue
                
            # Calculate mutation frequency per gene (Hugo_Symbol)
            gene_counts = (
                qualifying_subset.groupby('Hugo_Symbol')['Tumor_Sample_Barcode']
                .nunique()
                .reset_index()
            )
            gene_counts.columns = ['Gene', 'Mutated_Samples']
            gene_counts['Total_Samples'] = total_samples
            gene_counts['Frequency_Pct'] = (gene_counts['Mutated_Samples'] / total_samples) * 100.0
            
            # Sort by frequency descending
            gene_counts = gene_counts.sort_values(by='Frequency_Pct', ascending=False).reset_index(drop=True)
            
            # Save all mutated genes for this cohort to output list
            for _, row in gene_counts.iterrows():
                all_summary_data.append({
                    'Cohort': ds_name,
                    'Hugo_Symbol': row['Gene'],
                    'Mutated_Samples': row['Mutated_Samples'],
                    'Total_Samples': row['Total_Samples'],
                    'Mutation_Frequency_Pct': row['Frequency_Pct']
                })
                
        except Exception as e:
            print(f"Error processing cohort {ds_name} during mutation frequency calculation: {e}")
            
    if all_summary_data:
        return pd.DataFrame(all_summary_data)
    else:
        return pd.DataFrame(columns=['Cohort', 'Hugo_Symbol', 'Mutated_Samples', 'Total_Samples', 'Mutation_Frequency_Pct'])
