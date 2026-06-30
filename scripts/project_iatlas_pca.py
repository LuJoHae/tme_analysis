import os
import pathlib
import datalair
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score, accuracy_score
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader

from ici_datasets.other_datasets import HarmonizedTcgaIAtlas
from ici_datasets.bagaev_datasets import Signature
import icir
from gene_utils import read_gene_sets

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

class ResponseNet(nn.Module):
    def __init__(self, input_dim):
        super(ResponseNet, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(128, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(64, 1)
        )
    def forward(self, x):
        return self.net(x)

def plot_cohort_curves(y_trues, y_preds, cohorts, output_dir, prefix, title_suffix):
    unique_cohorts = np.unique(cohorts)
    colors = sns.color_palette("tab20", n_colors=len(unique_cohorts))
    
    # ROC Plot
    fig_roc, ax_roc = plt.subplots(figsize=(10, 8))
    # PR Plot
    fig_pr, ax_pr = plt.subplots(figsize=(10, 8))
    
    for i, cohort in enumerate(unique_cohorts):
        mask = (cohorts == cohort)
        y_t = y_trues[mask]
        y_p = y_preds[mask]
        
        if len(np.unique(y_t)) > 1:
            # ROC
            fpr, tpr, _ = roc_curve(y_t, y_p)
            roc_auc = auc(fpr, tpr)
            ax_roc.plot(fpr, tpr, color=colors[i], lw=2, label=f"{cohort} (AUC = {roc_auc:.2f})")
            
            # PR
            prec, rec, _ = precision_recall_curve(y_t, y_p)
            pr_auc = average_precision_score(y_t, y_p)
            ax_pr.plot(rec, prec, color=colors[i], lw=2, label=f"{cohort} (AP = {pr_auc:.2f})")
        else:
            print(f"Skipping curves for {cohort} ({title_suffix}) as it lacks binary response diversity.")

    ax_roc.plot([0, 1], [0, 1], 'k--', lw=1)
    ax_roc.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], xlabel='False Positive Rate', ylabel='True Positive Rate', title=f"ROC by Cohort - {title_suffix}")
    ax_roc.legend(loc="lower right", fontsize=8, ncol=2)
    fig_roc.tight_layout()
    fig_roc.savefig(output_dir / f"cohort_roc_{prefix}.svg", dpi=300)
    plt.close(fig_roc)
    
    ax_pr.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], xlabel='Recall', ylabel='Precision', title=f"PR Curve by Cohort - {title_suffix}")
    ax_pr.legend(loc="lower left", fontsize=8, ncol=2)
    fig_pr.tight_layout()
    fig_pr.savefig(output_dir / f"cohort_pr_{prefix}.svg", dpi=300)
    plt.close(fig_pr)


def train_and_evaluate_nn(X, y, cohorts, output_dir, prefix, title):
    print(f"\n--- 5-Fold CV PyTorch NN for {title} ---")
    
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    y_preds_all = np.zeros(len(y))
    epochs = 50
    unique_cohorts = np.unique(cohorts)
    
    # Dictionary to track epoch validation accuracies: epoch -> cohort -> list of fold accuracies
    epoch_cohort_accs = {ep: {c: [] for c in unique_cohorts} for ep in range(epochs)}

    for fold_i, (train_index, test_index) in enumerate(skf.split(X, y)):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        cohort_test = cohorts[test_index]

        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        train_dataset = TensorDataset(torch.tensor(X_train_scaled, dtype=torch.float32), torch.tensor(y_train, dtype=torch.float32))
        test_dataset = TensorDataset(torch.tensor(X_test_scaled, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
        
        train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True, drop_last=True)
        test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
        
        model = ResponseNet(input_dim=X.shape[1]).to(device)
        criterion = nn.BCEWithLogitsLoss()
        optimizer = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-2)

        for epoch in range(epochs):
            model.train()
            for batch_X, batch_y in train_loader:
                batch_X, batch_y = batch_X.to(device), batch_y.to(device)
                optimizer.zero_grad()
                outputs = model(batch_X).squeeze()
                if outputs.dim() == 0: outputs = outputs.unsqueeze(0)
                if batch_y.dim() == 0: batch_y = batch_y.unsqueeze(0)
                loss = criterion(outputs, batch_y)
                loss.backward()
                optimizer.step()
                
            # Track epoch validation accuracy
            model.eval()
            y_val_preds = []
            with torch.no_grad():
                for batch_X, _ in test_loader:
                    batch_X = batch_X.to(device)
                    logits = model(batch_X).squeeze()
                    if logits.dim() == 0: logits = logits.unsqueeze(0)
                    probs = torch.sigmoid(logits).cpu().numpy()
                    y_val_preds.extend(probs)
            
            y_val_preds = np.array(y_val_preds)
            y_val_bin = (y_val_preds > 0.5).astype(int)
            
            for cohort in unique_cohorts:
                mask = (cohort_test == cohort)
                if np.sum(mask) > 0:
                    acc = accuracy_score(y_test[mask], y_val_bin[mask])
                    epoch_cohort_accs[epoch][cohort].append(acc)

        # Final predictions for this fold
        y_preds_all[test_index] = y_val_preds

    # 1. Plot Epoch vs Accuracy per Cohort
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = sns.color_palette("tab20", n_colors=len(unique_cohorts))
    
    for i, cohort in enumerate(unique_cohorts):
        mean_accs = []
        for ep in range(epochs):
            fold_accs = epoch_cohort_accs[ep][cohort]
            if len(fold_accs) > 0:
                mean_accs.append(np.mean(fold_accs))
            else:
                mean_accs.append(np.nan)
        ax.plot(range(1, epochs+1), mean_accs, color=colors[i], label=cohort, lw=2)
        
    ax.set(xlabel="Epoch", ylabel="Mean Validation Accuracy", title=f"Epoch Accuracy per Cohort (5-Fold CV) - {title}")
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8, ncol=1)
    fig.tight_layout()
    fig.savefig(output_dir / f"epoch_acc_5fold_{prefix}.svg", dpi=300)
    plt.close(fig)

    # 2. Plot OOF ROC and PR Curves per Cohort
    plot_cohort_curves(y, y_preds_all, cohorts, output_dir, f"5fold_{prefix}", f"5-Fold CV ({prefix})")


def train_loco_cv(X, y, cohorts, output_dir, prefix, title):
    print(f"\n--- LOCO CV PyTorch NN for {title} ---")
    
    unique_cohorts = np.unique(cohorts)
    y_preds_all = np.zeros(len(y))
    
    for cohort in unique_cohorts:
        test_mask = (cohorts == cohort)
        train_mask = ~test_mask
        
        if np.sum(train_mask) == 0 or np.sum(test_mask) == 0:
            continue
            
        X_train, X_test = X[train_mask], X[test_mask]
        y_train, y_test = y[train_mask], y[test_mask]
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        train_dataset = TensorDataset(torch.tensor(X_train_scaled, dtype=torch.float32), torch.tensor(y_train, dtype=torch.float32))
        test_dataset = TensorDataset(torch.tensor(X_test_scaled, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
        
        train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True, drop_last=True)
        test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
        
        model = ResponseNet(input_dim=X.shape[1]).to(device)
        criterion = nn.BCEWithLogitsLoss()
        optimizer = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-2)

        epochs = 50
        for epoch in range(epochs):
            model.train()
            for batch_X, batch_y in train_loader:
                batch_X, batch_y = batch_X.to(device), batch_y.to(device)
                optimizer.zero_grad()
                outputs = model(batch_X).squeeze()
                if outputs.dim() == 0: outputs = outputs.unsqueeze(0)
                if batch_y.dim() == 0: batch_y = batch_y.unsqueeze(0)
                loss = criterion(outputs, batch_y)
                loss.backward()
                optimizer.step()
                
        model.eval()
        y_val_preds = []
        with torch.no_grad():
            for batch_X, _ in test_loader:
                batch_X = batch_X.to(device)
                logits = model(batch_X).squeeze()
                if logits.dim() == 0: logits = logits.unsqueeze(0)
                probs = torch.sigmoid(logits).cpu().numpy()
                y_val_preds.extend(probs)
                
        y_preds_all[test_mask] = np.array(y_val_preds)
        
    # Plot LOCO ROC and PR Curves per Cohort
    plot_cohort_curves(y, y_preds_all, cohorts, output_dir, f"loco_{prefix}", f"LOCO CV ({prefix})")


def get_iatlas_clinical_data(lair, adata_iatlas):
    print("Loading original iAtlas clinical data, TMB, and Cancer Type from cBioPortal datasets...")
    responses_map = {}
    tmb_map = {}
    cancer_type_map = {}
    
    import sys
    import importlib
    import ici_datasets
    if "scripts" not in sys.path:
        sys.path.append("scripts")
    iatlas_tmb = importlib.import_module("iAtlas-TMB")
    from gene_utils import calculate_maf_tmb
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    iatlas_dataset_names = [s for s in dataset_class.get_dataset_name() if "iAtlas" in s]
    
    for ds_name in iatlas_dataset_names:
        try:
            data_dir = iatlas_tmb.get_dataset_dir(lair, dataset_class, ds_name)
            clin_sample = data_dir / "data_clinical_sample.txt"
            clin_patient = data_dir / "data_clinical_patient.txt"
            mut_file = data_dir / "data_mutations.txt"
            
            if mut_file.exists():
                try:
                    data_mutations = pd.read_csv(mut_file, sep="\t", low_memory=False)
                    tmb_results = calculate_maf_tmb(data_mutations)
                    for _, row in tmb_results.iterrows():
                        tmb_map[row['Tumor_Sample_Barcode']] = row['TMB_Score']
                except Exception as e:
                    pass
            
            patient_resp = {}
            if clin_patient.exists():
                df_p = pd.read_csv(clin_patient, sep='\t', comment='#')
                resp_col = 'RESPONDER' if 'RESPONDER' in df_p.columns else ('CLINICAL_BENEFIT' if 'CLINICAL_BENEFIT' in df_p.columns else None)
                if resp_col and 'PATIENT_ID' in df_p.columns:
                    for _, row in df_p.iterrows():
                        val = row[resp_col]
                        if pd.notna(val):
                            if str(val).lower() == 'true':
                                patient_resp[str(row['PATIENT_ID'])] = 1
                            elif str(val).lower() == 'false':
                                patient_resp[str(row['PATIENT_ID'])] = 0
            
            if clin_sample.exists():
                df_s = pd.read_csv(clin_sample, sep='\t', comment='#')
                if 'SAMPLE_ID' in df_s.columns:
                    for _, row in df_s.iterrows():
                        sample_id = str(row['SAMPLE_ID'])
                        resp = np.nan
                        if 'PATIENT_ID' in row and str(row['PATIENT_ID']) in patient_resp:
                            resp = patient_resp[str(row['PATIENT_ID'])]
                        responses_map[sample_id] = resp
                        
                        if 'CANCER_TYPE' in df_s.columns:
                            cancer_type_map[sample_id] = row['CANCER_TYPE']
                            
                        # Fallback for TMB if MAF parsing failed or TMB is pre-calculated
                        if 'TMB_NONSYNONYMOUS' in df_s.columns and pd.notna(row['TMB_NONSYNONYMOUS']):
                            if sample_id not in tmb_map:
                                val = row['TMB_NONSYNONYMOUS']
                                if str(val).strip().upper() not in ["NA", "NAN", ""]:
                                    try:
                                        tmb_map[sample_id] = float(val)
                                    except ValueError:
                                        pass
                
        except Exception as e:
            pass

    mapped_responses = []
    mapped_tmb = []
    mapped_cancer = []
    for sample in adata_iatlas.obs_names:
        clean_sample = sample
        if clean_sample not in responses_map and '-' in clean_sample:
            clean_sample = clean_sample.rsplit('-', 1)[0]
        mapped_responses.append(responses_map.get(clean_sample, np.nan))
        mapped_tmb.append(tmb_map.get(clean_sample, np.nan))
        mapped_cancer.append(cancer_type_map.get(clean_sample, 'Unknown'))
    
    return np.array(mapped_responses), np.array(mapped_tmb), np.array(mapped_cancer)


def project_and_plot(adata_tcga, adata_iatlas, output_dir, prefix, title):
    print(f"\n--- Processing {title} ---")
    
    X_tcga = adata_tcga.X.toarray() if hasattr(adata_tcga.X, 'toarray') else adata_tcga.X
    X_iatlas = adata_iatlas.X.toarray() if hasattr(adata_iatlas.X, 'toarray') else adata_iatlas.X

    print("Standardizing data based on TCGA distribution for PCA...")
    scaler = StandardScaler(with_mean=True, with_std=True)
    X_tcga_scaled = scaler.fit_transform(X_tcga)
    X_iatlas_scaled = scaler.transform(X_iatlas)

    pca = PCA(n_components=50, random_state=42)
    tcga_pca = pca.fit_transform(X_tcga_scaled)
    iatlas_pca = pca.transform(X_iatlas_scaled)

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.scatter(tcga_pca[:, 0], tcga_pca[:, 1], color='grey', alpha=0.2, s=15, label='TCGA (Background)', edgecolors='none')

    iatlas_cohorts = adata_iatlas.obs['dataset'].values
    unique_cohorts = np.unique(iatlas_cohorts)
    colors = sns.color_palette("tab20", n_colors=len(unique_cohorts))
    
    for i, cohort in enumerate(unique_cohorts):
        mask = (iatlas_cohorts == cohort)
        ax.scatter(iatlas_pca[mask, 0], iatlas_pca[mask, 1], color=colors[i], alpha=0.7, s=25, label=cohort, edgecolors='none')

    ax.set_title(title, fontsize=16)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)", fontsize=12)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)", fontsize=12)
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0., fontsize=10, markerscale=2)
    
    pca_output = output_dir / f"projected_pca_{prefix}.svg"
    fig.tight_layout()
    fig.savefig(pca_output, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # NN Training Preparation
    resp = adata_iatlas.obs['response'].values
    tmb = adata_iatlas.obs['TMB_Score'].values
    valid_mask = (~pd.isna(resp)) & (~pd.isna(tmb))
    
    if not valid_mask.any():
        print("No valid response and TMB labels found! Skipping NN training.")
        return
        
    # Append TMB and Cancer Type as additional features
    X_nn = X_iatlas[valid_mask]
    tmb_nn = adata_iatlas.obs['TMB_Score'].values[valid_mask].reshape(-1, 1)
    
    cancer_types_nn = adata_iatlas.obs['cancer_type'].values[valid_mask]
    cancer_dummies = pd.get_dummies(cancer_types_nn, dummy_na=False).values
    
    X_nn = np.hstack([X_nn, tmb_nn, cancer_dummies])
    
    y_nn = resp[valid_mask].astype(int)
    cohorts_nn = iatlas_cohorts[valid_mask]
    
    print(f"Training NN with {X_nn.shape[0]} valid samples and {X_nn.shape[1]} features (including TMB)...")
    train_and_evaluate_nn(X_nn, y_nn, cohorts_nn, output_dir, prefix, title)
    train_loco_cv(X_nn, y_nn, cohorts_nn, output_dir, prefix, title)


def main():
    output_dir = pathlib.Path("output/project_iatlas_pca")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    lair = datalair.Lair("/storage/halu/lair")
    lair.assert_ok_satus()
    
    print("Loading Bagaev Signature Genes...")
    ds_sig = Signature()
    lair.safe_derive(ds_sig)
    filepaths_sig = lair.get_dataset_filepaths(ds_sig)
    bagaev_signature = read_gene_sets(filepaths_sig["gene_signatures.gmt"])
    bagaev_genes = sorted(list(set.union(*[x.genes for x in bagaev_signature.values()])))
    print(f"Found {len(bagaev_genes)} unique Bagaev genes.")

    ds_harmonized = HarmonizedTcgaIAtlas()
    ds_path = lair.get_path(ds_harmonized)
    
    if not ds_path.exists():
        raise FileNotFoundError(f"Harmonized datasets not found at {ds_path}. Please run `scripts/TCGA-background.py` first.")

    def process_pair(tcga_file, iatlas_file, prefix, title):
        adata_tcga = ad.read_h5ad(ds_path / tcga_file)
        adata_iatlas = ad.read_h5ad(ds_path / iatlas_file)
        
        responses, tmb_scores, cancer_types = get_iatlas_clinical_data(lair, adata_iatlas)
        adata_iatlas.obs['response'] = responses
        adata_iatlas.obs['TMB_Score'] = tmb_scores
        adata_iatlas.obs['cancer_type'] = cancer_types
        
        valid_cohorts = []
        for cohort in adata_iatlas.obs['dataset'].unique():
            mask = adata_iatlas.obs['dataset'] == cohort
            cohort_tmb = adata_iatlas.obs.loc[mask, 'TMB_Score']
            
            # If all samples in this cohort are NaN, median will be NaN
            if np.isnan(cohort_tmb.median()):
                print(f"Dropping cohort {cohort} as it completely lacks mutation data.")
                continue
                
            valid_cohorts.append(cohort)
            
        adata_iatlas = adata_iatlas[adata_iatlas.obs['dataset'].isin(valid_cohorts)].copy()
            
        target_genes = [g for g in bagaev_genes if g in adata_tcga.var_names]
        adata_tcga = adata_tcga[:, target_genes].copy()
        adata_iatlas = adata_iatlas[:, target_genes].copy()
        
        project_and_plot(adata_tcga, adata_iatlas, output_dir, prefix, title)

    # --- Uncorrected Pipeline ---
    process_pair("tcga_harmonized_uncorrected.h5ad", "iatlas_harmonized_uncorrected.h5ad", "uncorrected", "Uncorrected")

    # --- ComBat Corrected Pipeline ---
    process_pair("tcga_harmonized_combat.h5ad", "iatlas_harmonized_combat.h5ad", "combat", "ComBat Corrected")

if __name__ == "__main__":
    main()
