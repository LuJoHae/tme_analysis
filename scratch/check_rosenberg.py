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
        if "Rosenberg" in ds_name or "Cloughesy" in ds_name or "Gide" in ds_name or "Anders" in ds_name:
            print(f"\nChecking dataset: {ds_name}")
            try:
                data_dir = iatlas_tmb.get_dataset_dir(lair, dataset_class, ds_name)
                print(f"Data directory: {data_dir}")
                print("Files in directory:")
                if data_dir.exists():
                    for f in data_dir.iterdir():
                        print(f"  - {f.name}")
                else:
                    print("Directory does not exist!")
                    
                # Check clinical sample
                clin_sample = data_dir / "data_clinical_sample.txt"
                if clin_sample.exists():
                    df_s = pd.read_csv(clin_sample, sep='\t', comment='#')
                    for col in df_s.columns:
                        if 'TMB' in col.upper() or 'MUT' in col.upper() or 'BURDEN' in col.upper() or 'LOAD' in col.upper():
                            print(f"  Found potential TMB column in clinical sample: {col}")
                            
                # Check clinical patient
                clin_patient = data_dir / "data_clinical_patient.txt"
                if clin_patient.exists():
                    df_p = pd.read_csv(clin_patient, sep='\t', comment='#')
                    for col in df_p.columns:
                        if 'TMB' in col.upper() or 'MUT' in col.upper() or 'BURDEN' in col.upper() or 'LOAD' in col.upper():
                            print(f"  Found potential TMB column in clinical patient: {col}")
                            
            except Exception as e:
                print(f"Error: {e}")

if __name__ == "__main__":
    main()
