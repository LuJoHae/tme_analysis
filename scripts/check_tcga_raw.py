import datalair
from ici_datasets.other_datasets import load_tcga
import numpy as np

def main():
    lair = datalair.Lair("/storage/halu/lair")
    adata = load_tcga(lair)
    print("TCGA shape:", adata.shape)
    
    # Check X type and some values
    X_sample = adata.X[:5, :5]
    if hasattr(X_sample, "toarray"):
        X_sample = X_sample.toarray()
    print("Sample raw values (first 5x5):\n", X_sample)
    
    # Calculate row sums of raw X
    row_sums = []
    for i in range(min(5, adata.shape[0])):
        row = adata.X[i]
        if hasattr(row, "toarray"):
            row = row.toarray()
        row_sums.append(row.sum())
    print("Row sums of first 5 samples:", row_sums)
    
    # Check min/max/mean
    non_zero = adata.X[adata.X > 0]
    if hasattr(non_zero, "toarray"):
        non_zero = non_zero.toarray()
    print("Min positive value:", np.min(non_zero))
    
if __name__ == "__main__":
    main()
