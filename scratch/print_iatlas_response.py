import datalair
import icir
import anndata as ad

def main():
    lair = datalair.Lair("/storage/halu/lair")
    print("Loading original iAtlas datasets...")
    for name, path in lair.get_dataset_filepaths(icir.datasets.ImmuneCheckpointTherapyResponseProcessedGeneNormalizedClinicalDataNormalized()).items():
        adata = ad.read_h5ad(path, backed='r')
        resp_cols = [c for c in adata.obs.columns if 'resp' in c.lower() or 'clin' in c.lower() or 'benefit' in c.lower()]
        print(f"Dataset {name}: {resp_cols}")
        if len(resp_cols) > 0:
            print(f"  Sample values from {resp_cols[0]}:", list(adata.obs[resp_cols[0]].unique()[:5]))

if __name__ == "__main__":
    main()
