import sys
import importlib
import pandas as pd
import datalair
import ici_datasets

sys.path.append("scripts")
iatlas_tmb = importlib.import_module("iAtlas-TMB")

def main():
    lair = datalair.Lair("/storage/halu/lair")
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    iatlas_dataset_names = [s for s in dataset_class.get_dataset_name() if "iAtlas" in s]
    
    for ds_name in iatlas_dataset_names:
        print(f"\nDataset: {ds_name}")
        data_dir = iatlas_tmb.get_dataset_dir(lair, dataset_class, ds_name)
        
        for fname in ["data_clinical_patient.txt", "data_clinical_sample.txt"]:
            cf = data_dir / fname
            if cf.exists():
                try:
                    df = pd.read_csv(cf, sep='\t', comment='#')
                    found_cols = [c for c in df.columns if 'cancer' in c.lower() or 'disease' in c.lower() or 'type' in c.lower()]
                    if found_cols:
                        print(f"  {fname}: Found cols {found_cols}")
                        for col in found_cols:
                            print(f"    {col}: {df[col].unique()[:5]}")
                except Exception as e:
                    pass

if __name__ == "__main__":
    main()
