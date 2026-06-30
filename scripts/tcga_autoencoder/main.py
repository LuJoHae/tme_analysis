import os
import torch
from torch.utils.data import DataLoader, TensorDataset
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import torch.optim as optim
import torch.nn.functional as F

import tcga
import datalair
from ici_datasets.other_datasets import load_tcga
from .model import TCGAVAE
from .train import train_vae
from .plot import plot_loss, plot_latent_space, plot_reconstruction_error, plot_latent_density

def run_pipeline():
    print("Loading TCGA Data...")
    lair = datalair.Lair("/storage/halu/lair")
    
    # Safe-derive TCGA dataset before loading
    print("Safe-deriving TCGA dataset...")
    ds = tcga.AllProjectsAdata()
    lair.safe_derive(ds)
    
    adata = load_tcga(lair)
    adata.var_names_make_unique()
    print(f"Original shape: {adata.shape}")
    
    print("Converting gene expression matrix to float32...")
    adata.X = adata.X.astype(np.float32)
    
    # Preprocessing
    print("Normalizing and logging data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
        
    print("Selecting highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    adata = adata[:, adata.var.highly_variable].copy()
    
    print("Scaling data...")
    sc.pp.scale(adata, max_value=10)
    
    print("Running PCA with 500 components...")
    sc.tl.pca(adata, n_comps=500)
    print(f"Shape after PCA: {adata.obsm['X_pca'].shape}")
    
    X = adata.obsm['X_pca']
    
    label_col = 'dataset' if 'dataset' in adata.obs else ('project' if 'project' in adata.obs else None)
    if label_col:
        labels = adata.obs[label_col].values
    else:
        labels = np.array(['Unknown'] * adata.shape[0])
        
    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    y_encoded = le.fit_transform(labels)
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    tensor_x = torch.Tensor(X)
    tensor_y = torch.tensor(y_encoded, dtype=torch.long)
    dataset = TensorDataset(tensor_x, tensor_y)
    dataloader = DataLoader(dataset, batch_size=128, shuffle=True)
    
    input_dim = X.shape[1]
    model = TCGAVAE(input_dim=input_dim).to(device)
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    
    print("Starting training...")
    import time
    start_time = time.time()
    epochs = 100
    
    output_dir = Path("output/tcga_autoencoder")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    loss_history = train_vae(model, dataloader, optimizer, epochs, device, 
                             eval_data=tensor_x, labels=labels, output_dir=output_dir)
    end_time = time.time()
    
    runtime = end_time - start_time
    print(f"Training completed in {runtime:.2f} seconds.")
    
    print("Running inference...")
    model.eval()
    with torch.no_grad():
        tensor_x_device = tensor_x.to(device)
        recon_x, mu, logvar = model(tensor_x_device)
        
        z = mu.cpu().numpy()
        recon_x_cpu = recon_x.cpu()
        
        mse_errors = F.mse_loss(recon_x_cpu, tensor_x, reduction='none').mean(dim=1).numpy()
        
    print("Generating remaining plots...")
    
    plot_loss(loss_history, output_dir)
    plot_latent_space(z, labels, output_dir)
    plot_reconstruction_error(mse_errors, output_dir)
    plot_latent_density(z, output_dir)
    
    import csv
    with open(output_dir / "runtime.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["runtime_seconds"])
        writer.writerow([runtime])
        
    torch.save(model.state_dict(), output_dir / "weights.pt")
    print(f"Pipeline completed. Results saved to {output_dir}")

if __name__ == "__main__":
    run_pipeline()
