import sys
import importlib.util
from pathlib import Path
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss, log_loss
from sklearn.model_selection import KFold
import anndata as ad
import scanpy as sc
import random

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

get_dataset_dir = iAtlas_TMB.get_dataset_dir

TISSUE_MAP = {
    'Riaz-iAtlas': 'Melanoma',
    'Liu-iAtlas': 'Melanoma',
    'Gide-iAtlas': 'Melanoma',
    'Hugo-iAtlas': 'Melanoma',
    'Choueiri-iAtlas': 'RCC',
    'McDermott-iAtlas': 'RCC',
    'Cloughesy-iAtlas': 'GBM',
    'Padron-iAtlas': 'PAAD',
}

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
            
            df_mrna = np.log1p(df_mrna)
            
            df_clinical = df_clinical[df_clinical['response'].isin(['R', 'NR'])]
            if df_clinical.empty: continue
                
            common_samples = df_mrna.columns.intersection(df_clinical.index)
            if len(common_samples) == 0: continue
                
            df_mrna = df_mrna[common_samples]
            df_clinical = df_clinical.loc[common_samples]
            
            envs_X.append(df_mrna.T)
            envs_Y.append(df_clinical['response'].map({'R': 1, 'NR': 0}).values)
            env_names.append(ds_name)
            print(f"Loaded {ds_name}: {len(common_samples)} samples")
        except Exception as e:
            print(f"Skipped {ds_name}: {e}")

    if not envs_X: raise ValueError("No valid iAtlas datasets found.")
        
    common_genes = set(envs_X[0].columns)
    for X in envs_X[1:]:
        common_genes = common_genes.intersection(X.columns)
    common_genes = list(common_genes.intersection(bagaev_genes))
    print(f"\nUsing {len(common_genes)} common Bagaev signature genes.")

    for i in range(len(envs_X)):
        envs_X[i] = envs_X[i][common_genes].values.astype(np.float32)

    # ComBat Harmonization
    all_X = np.vstack(envs_X)
    cohort_labels = np.concatenate([[env_names[i]] * len(envs_X[i]) for i in range(len(envs_X))])
    
    try:
        adata = ad.AnnData(X=all_X)
        adata.obs['cohort'] = pd.Categorical(cohort_labels)
        sc.pp.combat(adata, key='cohort')
        all_X_harmonized = adata.X.astype(np.float32)
        print("ComBat Harmonization successful.")
    except Exception as e:
        print(f"ComBat failed: {e}. Falling back to unharmonized data.")
        all_X_harmonized = all_X.astype(np.float32)

    # Compute Meta-Features (Signature Scores)
    sig_features = []
    for sig_name, sig_obj in bagaev_signature.items():
        sig_indices = [idx for idx, g in enumerate(common_genes) if g in sig_obj.genes]
        if len(sig_indices) > 0:
            sig_scores = np.mean(all_X_harmonized[:, sig_indices], axis=1)
            sig_features.append(sig_scores)
    
    all_X_meta = np.column_stack(sig_features)
    print(f"Extracted {all_X_meta.shape[1]} Bagaev signature meta-features.")

    envs_X_raw = []
    envs_X_meta = []
    idx_start = 0
    for i in range(len(envs_X)):
        n = len(envs_X[i])
        envs_X_raw.append(all_X_harmonized[idx_start:idx_start+n])
        envs_X_meta.append(all_X_meta[idx_start:idx_start+n])
        idx_start += n

    return envs_X_raw, envs_X_meta, envs_Y, env_names, len(common_genes), all_X_meta.shape[1]

class FocalLoss(nn.Module):
    def __init__(self, gamma=2.0):
        super().__init__()
        self.gamma = gamma

    def forward(self, logits, targets, pos_weight=None):
        bce_loss = F.binary_cross_entropy_with_logits(logits, targets, reduction='none', pos_weight=pos_weight)
        probs = torch.sigmoid(logits)
        pt = torch.where(targets == 1, probs, 1 - probs)
        focal_loss = (1 - pt) ** self.gamma * bce_loss
        return focal_loss.mean()

class BCELossWrap(nn.Module):
    def forward(self, logits, targets, pos_weight=None):
        return F.binary_cross_entropy_with_logits(logits, targets, pos_weight=pos_weight)

class LinearModel(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Linear(input_dim, 1)
    def forward(self, x):
        return self.net(x)

class MLPModel(nn.Module):
    def __init__(self, input_dim, hidden_size=32):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_size),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(hidden_size, 1)
        )
    def forward(self, x):
        return self.net(x)

def get_model(model_type, input_dim, hidden_size=32):
    if model_type == 'linear': return LinearModel(input_dim)
    return MLPModel(input_dim, hidden_size)

def get_loss_fn(loss_type):
    if loss_type == 'focal': return FocalLoss()
    return BCELossWrap()

def compute_irm_penalty(logits, y, pos_weight=None):
    scale = torch.tensor(1.).requires_grad_()
    loss = F.binary_cross_entropy_with_logits(logits * scale, y.float(), pos_weight=pos_weight)
    grad = torch.autograd.grad(loss, [scale], create_graph=True)[0]
    return torch.sum(grad**2)

def train_irm(envs_train, config, input_dim=0, lr=1e-3, steps=100, penalty_weight=10.0, penalty_anneal_iters=20, weight_decay=1e-2):
    model = get_model(config['model'], input_dim, hidden_size=config.get('hidden_size', 32))
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    criterion = get_loss_fn(config['loss'])
    
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
            nll = criterion(logits, Y_t.float(), pos_weight=pos_weight)
            
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
        
    l1_norm = sum(p.abs().sum().item() for p in model.parameters())
    l2_norm = sum((p**2).sum().item() for p in model.parameters())**0.5
    
    return model, {'l1_norm': l1_norm, 'l2_norm': l2_norm, 'final_penalty': train_penalty.item() if isinstance(train_penalty, torch.Tensor) else train_penalty, 'dro_q': None}

def train_group_dro(envs_train, config, input_dim=0, lr=1e-3, steps=100, step_size=0.01, weight_decay=1e-2, **kwargs):
    model = get_model(config['model'], input_dim, hidden_size=config.get('hidden_size', 32))
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    criterion = get_loss_fn(config['loss'])
    
    num_envs = len([X for X, Y in envs_train if len(X) > 0])
    if num_envs == 0: return model, {}
    q = torch.ones(num_envs) / num_envs
    
    all_Y = np.concatenate([Y for X, Y in envs_train if len(Y) > 0])
    num_pos = np.sum(all_Y == 1)
    num_neg = np.sum(all_Y == 0)
    pos_weight = torch.tensor([num_neg / num_pos], dtype=torch.float32) if num_pos > 0 else torch.tensor([1.0])

    model.train()
    for step in range(steps):
        losses = []
        for X, Y in envs_train:
            if len(X) == 0: continue
            X_t = torch.tensor(X)
            Y_t = torch.tensor(Y).view(-1, 1)
            
            logits = model(X_t)
            nll = criterion(logits, Y_t.float(), pos_weight=pos_weight)
            losses.append(nll)
            
        losses = torch.stack(losses)
        
        with torch.no_grad():
            q = q * torch.exp(step_size * losses)
            q = q / q.sum()
            
        loss = torch.sum(q * losses)
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
    l1_norm = sum(p.abs().sum().item() for p in model.parameters())
    l2_norm = sum((p**2).sum().item() for p in model.parameters())**0.5
    q_dict = {f"q_env_{i}": val.item() for i, val in enumerate(q)}
    
    return model, {'l1_norm': l1_norm, 'l2_norm': l2_norm, 'final_penalty': None, 'dro_q': q_dict}

def train_vrex(envs_train, config, input_dim=0, lr=1e-3, steps=100, vrex_penalty=10.0, vrex_anneal_iters=20, weight_decay=1e-2, **kwargs):
    model = get_model(config['model'], input_dim, hidden_size=config.get('hidden_size', 32))
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    criterion = get_loss_fn(config['loss'])
    
    all_Y = np.concatenate([Y for X, Y in envs_train if len(Y) > 0])
    num_pos = np.sum(all_Y == 1)
    num_neg = np.sum(all_Y == 0)
    pos_weight = torch.tensor([num_neg / num_pos], dtype=torch.float32) if num_pos > 0 else torch.tensor([1.0])

    model.train()
    for step in range(steps):
        losses = []
        for X, Y in envs_train:
            if len(X) == 0: continue
            X_t = torch.tensor(X)
            Y_t = torch.tensor(Y).view(-1, 1)
            
            logits = model(X_t)
            nll = criterion(logits, Y_t.float(), pos_weight=pos_weight)
            losses.append(nll)
            
        losses = torch.stack(losses)
        mean_loss = losses.mean()
        var_loss = losses.var() if len(losses) > 1 else torch.tensor(0.0)
        
        weight = vrex_penalty if step >= vrex_anneal_iters else (1.0 if vrex_penalty > 0 else 0.0)
        loss = mean_loss + weight * var_loss
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
    l1_norm = sum(p.abs().sum().item() for p in model.parameters())
    l2_norm = sum((p**2).sum().item() for p in model.parameters())**0.5
    
    return model, {'l1_norm': l1_norm, 'l2_norm': l2_norm, 'final_penalty': var_loss.item() if isinstance(var_loss, torch.Tensor) else var_loss, 'dro_q': None}


def evaluate(model, X, Y, return_probs=False, mc_dropout=False):
    if mc_dropout:
        model.train() # keep dropout enabled
    else:
        model.eval()
        
    with torch.no_grad():
        if len(X) == 0:
            if return_probs: return np.nan, np.nan, np.nan, np.nan, []
            return np.nan, np.nan, np.nan, np.nan
        X_t = torch.tensor(X)
        
        if mc_dropout:
            all_probs = []
            for _ in range(20):
                logits = model(X_t)
                all_probs.append(torch.sigmoid(logits).numpy().flatten())
            probs = np.mean(all_probs, axis=0)
        else:
            logits = model(X_t)
            probs = torch.sigmoid(logits).numpy().flatten()
            
    if len(np.unique(Y)) > 1:
        roc = roc_auc_score(Y, probs)
        pr = average_precision_score(Y, probs)
        brier = brier_score_loss(Y, probs)
        ll = log_loss(Y, probs, labels=[0, 1])
    else:
        roc, pr, brier, ll = np.nan, np.nan, np.nan, np.nan
    if return_probs:
        return roc, pr, brier, ll, probs
    return roc, pr, brier, ll

def run_loco_cv(envs_X, envs_Y, env_names, input_dim, config):
    all_results = []
    mode = config['eval']
    
    if mode == 'erm':
        train_fn = train_irm
        hparams = {'lr': config['lr'], 'penalty_weight': 0.0, 'weight_decay': config['weight_decay']}
    elif mode == 'irm':
        train_fn = train_irm
        hparams = {'lr': config['lr'], 'penalty_weight': config['penalty_weight'], 'weight_decay': config['weight_decay']}
    elif mode == 'dro':
        train_fn = train_group_dro
        hparams = {'lr': config['lr'], 'step_size': config['dro_step_size'], 'weight_decay': config['weight_decay']}
    elif mode == 'vrex':
        train_fn = train_vrex
        hparams = {'lr': config['lr'], 'vrex_penalty': config['vrex_penalty'], 'weight_decay': config['weight_decay']}

    for test_idx in range(len(envs_X)):
        test_env_X = envs_X[test_idx]
        test_env_Y = envs_Y[test_idx]
        test_name = env_names[test_idx]
        
        tissue = TISSUE_MAP.get(test_name, "Unknown")
            
        train_envs = [(envs_X[i], envs_Y[i]) for i in range(len(envs_X)) if i != test_idx]
        if len(train_envs) < 1:
            continue
            
        model, stats = train_fn(train_envs, config, input_dim=input_dim, **hparams)
        
        train_metrics = {}
        for env_idx, (tr_X, tr_Y) in enumerate(train_envs):
            t_roc, t_pr, t_brier, t_ll = evaluate(model, tr_X, tr_Y, mc_dropout=config['mc_dropout'])
            train_metrics[f'train_env_{env_idx}_roc'] = t_roc
            train_metrics[f'train_env_{env_idx}_pr'] = t_pr
            
        test_auc, test_pr, test_brier, test_ll, test_probs = evaluate(model, test_env_X, test_env_Y, return_probs=True, mc_dropout=config['mc_dropout'])
        
        res_dict = {
            'Test_Cohort': test_name, 'Tissue': tissue, 
            'ROC_AUC': test_auc, 'PR_AUC': test_pr, 'Brier': test_brier, 'LogLoss': test_ll,
            'probs': test_probs, 'labels': test_env_Y,
            'L1_Norm': stats.get('l1_norm'), 'L2_Norm': stats.get('l2_norm'),
            'Final_Penalty': stats.get('final_penalty')
        }
        res_dict.update(train_metrics)
        if stats.get('dro_q') is not None:
            res_dict.update(stats['dro_q'])
            
        all_results.append(res_dict)
    
    return all_results

if __name__ == "__main__":
    envs_X_raw, envs_X_meta, envs_Y, env_names, num_raw, num_meta = load_data()
    
    # Generate random configurations
    random.seed(42)
    np.random.seed(42)
    torch.manual_seed(42)
    
    param_probs = {
        'model': {'mlp': 1.0},
        'features': {'raw': 1.0},
        'loss': {'bce': 1.0},
        'eval': {'erm': 0.33, 'irm': 0.33, 'dro': 0.34},
        'mc_dropout': {True: 1.0},
        'lr': {1e-4: 0.2, 5e-4: 0.3, 1e-3: 0.3, 5e-3: 0.2},
        'weight_decay': {1e-3: 0.2, 1e-2: 0.5, 5e-2: 0.2, 1e-1: 0.1},
        'hidden_size': {16: 0.2, 32: 0.5, 64: 0.3},
        'penalty_weight': {1.0: 0.33, 10.0: 0.33, 100.0: 0.34},
        'dro_step_size': {0.01: 0.5, 0.1: 0.5},
        'vrex_penalty': {1.0: 0.5, 10.0: 0.5}
    }
    
    N_CONFIGS = 100
    configs_to_try = []
    for _ in range(N_CONFIGS):
        config = {}
        for param, probs in param_probs.items():
            choices = list(probs.keys())
            weights = list(probs.values())
            
            # Handle mc_dropout logic constraint
            if param == 'mc_dropout' and config.get('model') == 'linear':
                config[param] = False
                continue
            config[param] = random.choices(choices, weights=weights, k=1)[0]
        configs_to_try.append(config)
        
    print(f"\nGenerated {N_CONFIGS} random configurations to evaluate.")
    
    results = []
    for idx, config in enumerate(configs_to_try):
        print(f"\n=== Running Config {idx+1}/{N_CONFIGS}: {config} ===")
        envs_X = envs_X_raw if config['features'] == 'raw' else envs_X_meta
        input_dim = num_raw if config['features'] == 'raw' else num_meta
        
        res = run_loco_cv(envs_X, envs_Y, env_names, input_dim=input_dim, config=config)
        
        for r in res:
            r.update(config)
            results.append(r)
            
    out_dir = Path("output") / Path(__file__).stem
    out_dir.mkdir(parents=True, exist_ok=True)
    
    for r in results:
        r.pop('probs', None)
        r.pop('labels', None)
        
    out_file = out_dir / "random_search_results.csv"
    pd.DataFrame(results).to_csv(out_file, index=False)
    print(f"\nRandom search results saved to {out_file}")
