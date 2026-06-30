import sys
import importlib.util
from pathlib import Path
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

import datalair
from gene_utils import read_gene_sets
import ici_datasets
from ici_datasets.bagaev_datasets import Signature
from ici_datasets.cbioportal_datasets import CBioPortalDataset

# Dynamically load iAtlas-TMB to reuse data loading
spec = importlib.util.spec_from_file_location("iAtlas_TMB", Path(__file__).parent / "iAtlas-TMB.py")
iAtlas_TMB = importlib.util.module_from_spec(spec)
sys.modules["iAtlas_TMB"] = iAtlas_TMB
spec.loader.exec_module(iAtlas_TMB)

# === 1. Data Loading and Preprocessing ===
get_dataset_dir = iAtlas_TMB.get_dataset_dir

def load_and_transform_data_mrna(data_dir):
    return iAtlas_TMB.load_and_transform_data_mrna(data_dir)

def load_data_clinical(data_dir):
    filepath_clinical_sample = data_dir / "data_clinical_sample.txt"
    data_clinical_sample = pd.read_csv(filepath_clinical_sample, sep="\t", index_col=0, skiprows=4)
    normalized = {c.strip().lower(): c for c in data_clinical_sample.columns}
    response_col = None
    for cand in ["response", "best_response", "recist"]:
        if cand in normalized:
            response_col = normalized[cand]
            break
    if response_col is None:
        raise KeyError("No response column found")
    
    data_clinical_sample = data_clinical_sample.set_index("SAMPLE_ID")
    responses = data_clinical_sample[response_col].map({
        "Complete Response": "R", "Partial Response": "R",
        "Stable Disease": "NR", "Progressive Disease": "NR",
        "CR": "R", "PR": "R", "SD": "NR", "PD": "NR",
    })
    return pd.DataFrame({"response": responses})

def load_data():
    lair = datalair.Lair("/storage/halu/lair")
    ds = Signature()
    lair.safe_derive(ds)
    filepaths = lair.get_dataset_filepaths(ds)
    bagaev_signature = read_gene_sets(filepaths["gene_signatures.gmt"])
    bagaev_genes = set.union(*[x.genes for x in bagaev_signature.values()])

    iatlas_dataset_names = [s for s in CBioPortalDataset.get_dataset_name() if "iAtlas" in s]

    envs_X = []
    envs_Y = []
    env_names = []

    for ds_name in iatlas_dataset_names:
        try:
            data_dir = get_dataset_dir(lair, CBioPortalDataset, ds_name)
            df_clinical = load_data_clinical(data_dir)
            df_mrna = load_and_transform_data_mrna(data_dir)
            
            # Normalize to log1p TPM
            df_mrna = np.log1p(df_mrna)
            
            # Filter valid clinical responses
            df_clinical = df_clinical[df_clinical['response'].isin(['R', 'NR'])]
            if df_clinical.empty:
                continue
                
            common_samples = df_mrna.columns.intersection(df_clinical.index)
            if len(common_samples) == 0:
                continue
                
            df_mrna = df_mrna[common_samples]
            df_clinical = df_clinical.loc[common_samples]
            
            envs_X.append(df_mrna.T)
            envs_Y.append(df_clinical['response'].map({'R': 1, 'NR': 0}).values)
            env_names.append(ds_name)
            print(f"Loaded {ds_name}: {len(common_samples)} samples")
        except Exception as e:
            print(f"Skipped {ds_name}: {e}")

    # Isolate intersection of genes present in ALL environments AND Bagaev signature
    if not envs_X:
        raise ValueError("No valid iAtlas datasets found.")
        
    common_genes = set(envs_X[0].columns)
    for X in envs_X[1:]:
        common_genes = common_genes.intersection(X.columns)
    common_genes = list(common_genes.intersection(bagaev_genes))
    print(f"\nUsing {len(common_genes)} common Bagaev signature genes.")

    for i in range(len(envs_X)):
        X_val = envs_X[i][common_genes].values.astype(np.float32)
        # Local cohort standardization to mitigate batch effects
        scaler = StandardScaler()
        envs_X[i] = scaler.fit_transform(X_val)

    return envs_X, envs_Y, env_names, common_genes

# === 2. IRM PyTorch Implementation ===
class LinearModel(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Linear(input_dim, 1)
    def forward(self, x):
        return self.net(x)

def compute_irm_penalty(logits, y, pos_weight=None):
    scale = torch.tensor(1.).requires_grad_()
    loss = F.binary_cross_entropy_with_logits(logits * scale, y.float(), pos_weight=pos_weight)
    grad = torch.autograd.grad(loss, [scale], create_graph=True)[0]
    return torch.sum(grad**2)

def train_irm(envs_train, input_dim=0, lr=1e-3, steps=100, penalty_weight=1e4, penalty_anneal_iters=20):
    model = LinearModel(input_dim)
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    
    # Calculate global pos_weight for class imbalance
    all_Y = np.concatenate([Y for X, Y in envs_train if len(Y) > 0])
    num_pos = np.sum(all_Y == 1)
    num_neg = np.sum(all_Y == 0)
    pos_weight = torch.tensor([num_neg / num_pos], dtype=torch.float32) if num_pos > 0 else torch.tensor([1.0])

    model.train()
    for step in range(steps):
        train_nll = 0.
        train_penalty = 0.
        for X, Y in envs_train:
            if len(X) == 0: continue
            X_t = torch.tensor(X)
            Y_t = torch.tensor(Y).view(-1, 1)
            
            logits = model(X_t)
            nll = F.binary_cross_entropy_with_logits(logits, Y_t.float(), pos_weight=pos_weight)
            
            if penalty_weight > 0:
                penalty = compute_irm_penalty(logits, Y_t, pos_weight=pos_weight)
            else:
                penalty = 0.0
            
            train_nll += nll
            train_penalty += penalty
            
        train_nll /= len(envs_train)
        train_penalty /= len(envs_train)
        
        weight = penalty_weight if step >= penalty_anneal_iters else (1.0 if penalty_weight > 0 else 0.0)
        loss = train_nll + weight * train_penalty
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    return model

def evaluate(model, X, Y):
    model.eval()
    with torch.no_grad():
        if len(X) == 0:
            return np.nan, np.nan
        X_t = torch.tensor(X)
        logits = model(X_t)
        probs = torch.sigmoid(logits).numpy().flatten()
    if len(np.unique(Y)) > 1:
        roc = roc_auc_score(Y, probs)
        pr = average_precision_score(Y, probs)
    else:
        roc, pr = np.nan, np.nan
    return roc, pr

# === 3. Evaluation Schemes ===
def run_loco_cv(envs_X, envs_Y, env_names, input_dim, is_erm=False):
    all_results = []
    if is_erm:
        hparam_grid = [
            {'lr': 1e-3, 'penalty_weight': 0.0},
            {'lr': 1e-2, 'penalty_weight': 0.0}
        ]
    else:
        hparam_grid = [
            {'lr': 1e-3, 'penalty_weight': 1e2},
            {'lr': 1e-3, 'penalty_weight': 1e4},
            {'lr': 1e-2, 'penalty_weight': 1e2},
        ]

    for test_idx in range(len(envs_X)):
        test_env_X = envs_X[test_idx]
        test_env_Y = envs_Y[test_idx]
        test_name = env_names[test_idx]
        
        # Training environments
        train_envs = [(envs_X[i], envs_Y[i]) for i in range(len(envs_X)) if i != test_idx]
        if len(train_envs) < 2:
            print(f"Skipping {test_name}: Not enough environments for IRM.")
            continue
            
        print(f"\n--- LOCO CV: Testing on {test_name} ---")
        best_hparams = None
        best_val_auc = -1
        
        # Inner loop for HPO (Hold-out one cohort from training environments)
        for hparams in hparam_grid:
            val_auc_sum = 0
            valid_folds = 0
            for val_idx in range(len(train_envs)):
                inner_train_envs = [train_envs[i] for i in range(len(train_envs)) if i != val_idx]
                inner_val_X, inner_val_Y = train_envs[val_idx]
                
                if len(inner_train_envs) < 1: continue
                model = train_irm(inner_train_envs, input_dim=input_dim, **hparams)
                auc, _ = evaluate(model, inner_val_X, inner_val_Y)
                if not np.isnan(auc):
                    val_auc_sum += auc
                    valid_folds += 1
            avg_val_auc = val_auc_sum / valid_folds if valid_folds > 0 else 0
            if avg_val_auc > best_val_auc:
                best_val_auc = avg_val_auc
                best_hparams = hparams
                
        # Final train with best hparams
        model = train_irm(train_envs, input_dim=input_dim, **best_hparams)
        test_auc, test_pr = evaluate(model, test_env_X, test_env_Y)
        
        print(f"Test AUC: {test_auc:.4f}, Test PR: {test_pr:.4f} (Hparams: {best_hparams})")
        eval_name = 'ERM_LOCO' if is_erm else 'IRM_LOCO'
        all_results.append({'Evaluation': eval_name, 'Test_Cohort': test_name, 'ROC_AUC': test_auc, 'PR_AUC': test_pr})
    
    return all_results

def run_kfold_cv(envs_X, envs_Y, env_names, input_dim):
    all_results = []
    
    X_all = np.vstack(envs_X)
    Y_all = np.concatenate(envs_Y)
    cohort_labels = np.concatenate([[env_names[i]] * len(envs_X[i]) for i in range(len(envs_X))])
    
    print(f"\n--- K-Fold CV (5 Folds) ---")
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    for fold, (train_idx, test_idx) in enumerate(kf.split(X_all)):
        X_te, Y_te = X_all[test_idx], Y_all[test_idx]
        
        # Reconstruct environments based on cohorts
        train_cohorts = cohort_labels[train_idx]
        unique_cohorts = np.unique(train_cohorts)
        
        train_envs = []
        for c in unique_cohorts:
            idx = (train_cohorts == c)
            # Only use environment if it has both R and NR or at least some samples
            if len(np.unique(Y_all[train_idx][idx])) > 1:
                train_envs.append((X_all[train_idx][idx], Y_all[train_idx][idx]))
        
        if len(train_envs) < 2:
            train_envs = [(X_all[train_idx], Y_all[train_idx])]
            
        model = train_irm(train_envs, input_dim=input_dim, lr=1e-3, penalty_weight=1e2)
        test_auc, test_pr = evaluate(model, X_te, Y_te)
        
        print(f"Fold {fold} - Test AUC: {test_auc:.4f}, Test PR: {test_pr:.4f}")
        all_results.append({'Evaluation': 'IRM_5-Fold', 'Test_Cohort': f"Fold_{fold}", 'ROC_AUC': test_auc, 'PR_AUC': test_pr})

    return all_results

if __name__ == "__main__":
    envs_X, envs_Y, env_names, common_genes = load_data()
    
    results = []
    print("\n=== Running ERM Baseline (LOCO CV) ===")
    results.extend(run_loco_cv(envs_X, envs_Y, env_names, input_dim=len(common_genes), is_erm=True))
    
    print("\n=== Running IRM (LOCO CV) ===")
    results.extend(run_loco_cv(envs_X, envs_Y, env_names, input_dim=len(common_genes), is_erm=False))
    
    results.extend(run_kfold_cv(envs_X, envs_Y, env_names, input_dim=len(common_genes)))
    
    out_dir = Path("output") / Path(__file__).stem
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "irm_results.csv"
    pd.DataFrame(results).to_csv(out_file, index=False)
    print(f"\nResults saved to {out_file}")
