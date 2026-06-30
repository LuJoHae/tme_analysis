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
    
    for ds_name in iatlas_dataset_names[:3]:
        print(f"\nDataset: {ds_name}")
        data_dir = iatlas_tmb.get_dataset_dir(lair, dataset_class, ds_name)
        
        for fname in ["data_clinical_patient.txt", "data_clinical_sample.txt"]:
            cf = data_dir / fname
            if cf.exists():
                try:
                    df = pd.read_csv(cf, sep='\t', comment='#')
                    print(f"  {fname}:")
                    for col in ['RESPONSE', 'RESPONDER', 'CLINICAL_BENEFIT']:
                        if col in df.columns:
                            print(f"    {col}: {df[col].unique().tolist()}")
                except Exception as e:
                    print(f"  {fname}: Could not read - {e}")

if __name__ == "__main__":
    main()
