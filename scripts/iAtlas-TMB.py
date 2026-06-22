from itertools import product
import os
import sys
import datalair
import matplotlib.pyplot as plt
import pprint as pprint
from pandas import DataFrame
import ici_datasets # pyrefly: ignore [missing-import]
import seaborn as sns
from pathlib import Path
import numpy as np
from scipy.stats import mannwhitneyu
from pprint import pprint
from gene_utils import calculate_maf_tmb, read_gene_sets, ssgsea_formula
from pyensembl import EnsemblRelease
import pandas as pd
from pathlib import Path
from typing import Callable
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import CalibratedClassifierCV
from joblib import Parallel, delayed
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_curve, roc_auc_score, accuracy_score
from sklearn.decomposition import PCA
from sklearn.compose import ColumnTransformer
from ici_datasets.bagaev_datasets import Signature


def transform_expression_to_tpm(counts_df: pd.DataFrame, release: int = 110) -> pd.DataFrame:
    """
    Computes TPM from raw counts using pyensembl for local gene metrics lookup.

    Parameters:
    counts_df: DataFrame of raw counts (rows = Hugo symbols, cols = samples).
    release: Ensembl release version number.
    """
    # 1. Initialize the local Ensembl database manager
    data_source = EnsemblRelease(release)

    def extract_kb_span(symbol: str) -> float | None:
        """Maps a Hugo symbol to its genomic span in kilobases."""
        try:
            genes = data_source.genes_by_name(symbol)
            # Map to the first transcript variant if multiple exist, else return None
            return (genes[0].end - genes[0].start + 1) / 1e3 if genes else None
        except ValueError:
            return None

    # 2. Filter-map operation over the index keyspace using assignment expression
    gene_lengths = pd.Series({
        symbol: length
        for symbol in counts_df.index
        if (length := extract_kb_span(symbol)) is not None
    })

    # 3. Restrict matrix operation to the intersection of existing annotations
    common_genes = counts_df.index.intersection(gene_lengths.index)
    counts = counts_df.loc[common_genes]
    lengths = gene_lengths.loc[common_genes]

    # 4. Matrix Transformations: RPK followed by TPM normalization
    rpk = counts.div(lengths, axis=0)
    return rpk.div(rpk.sum(axis=0), axis=1) * 1e6

def transform_rpkm_to_tpm(rpkm_df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts an RPKM/FPKM normalized DataFrame to TPM.
    """
    # Normalize each column such that its sum is strictly 10^6
    tpm = rpkm_df.div(rpkm_df.sum(axis=0), axis=1) * 1e6

    return tpm

def transform_tpm_to_tpm(df: pd.DataFrame) -> pd.DataFrame: return df


def load_and_transform_data_mrna(data_dir: str | Path) -> pd.DataFrame:
    p_dir = Path(data_dir)
    dispatch_map: dict[str, tuple[Callable[[pd.DataFrame], pd.DataFrame], str]] = {
        "data_mrna_seq_expression.txt": (transform_expression_to_tpm, "none"),
        "data_mrna_seq_tpm.txt": (transform_tpm_to_tpm, "tpm"),
        "data_mrna_seq_rpkm.txt": (transform_rpkm_to_tpm, "rpkm"),
    }

    existing_files = tuple(filter(lambda f: (p_dir / f).is_file(), dispatch_map.keys()))
    if len(existing_files) != 1:
        raise ValueError(
            f"Invariant violation: Expected exactly 1 file, but found {len(existing_files)} in {p_dir}."
        )
    target = existing_files[0]
    transform_fn, normalization = dispatch_map[target]

    data_mrna = transform_fn(pd.read_csv(p_dir / target, sep='\t', index_col=0))
    data_mrna.attrs["original_normalization"] = normalization
    if data_mrna.isna().any(axis=None):
        raise ValueError("Data integrity violation: DataFrame contains NA values.")
    return data_mrna


def load_and_process_data(data_dir):
    data_clinical = load_data_clinical(data_dir)
    data_mrna_seq_expression = load_and_transform_data_mrna(data_dir)
    return data_clinical, data_mrna_seq_expression


def load_data_clinical(data_dir) -> DataFrame:
    # Known aliases for the response column across cBioPortal datasets
    RESPONSE_COL_CANDIDATES = [
    "RESPONSE", "Response", "response",
    "BEST_RESPONSE", "Best_Response", "best_response",
    "RECIST", "recist",
]

    RESPONSE_VALUE_MAP = {
        "Complete Response": "R",
        "Partial Response": "R",
        "Stable Disease": "NR",
        "Progressive Disease": "NR",
        # Add other common encodings as needed:
        "CR": "R",
        "PR": "R",
        "SD": "NR",
        "PD": "NR",
    }

    def _resolve_column(df, candidates):
        """Return the actual column name matching any candidate (case/space-insensitive)."""
        normalized = {c.strip().lower(): c for c in df.columns}
        for cand in candidates:
            key = cand.strip().lower()
            if key in normalized:
                return normalized[key]
        raise KeyError(
            f"None of {candidates} found. Available columns: {list(df.columns)}"
        )

    filepath_clinical_sample = data_dir / "data_clinical_sample.txt"
    data_clinical_sample = pd.read_csv(filepath_clinical_sample, sep="\t", index_col=0, skiprows=4)

    response_col = _resolve_column(data_clinical_sample, RESPONSE_COL_CANDIDATES)
    if response_col is None:
        raise KeyError(
            "No response column found. Available columns: "
            f"{list(data_clinical_sample.columns)}"
        )

    data_clinical_sample = data_clinical_sample.set_index("SAMPLE_ID")
    data_clinical_sample = data_clinical_sample[response_col]


    filepath_mutations = data_dir / "data_mutations.txt"
    if filepath_mutations.exists():
        data_mutations = pd.read_csv(data_dir / "data_mutations.txt", sep="\t", low_memory=False)
        tmb_results = calculate_maf_tmb(data_mutations)
        tmb_results = tmb_results.set_index("Tumor_Sample_Barcode")
        data_clinical_sample = pd.concat([tmb_results, data_clinical_sample], join="inner", axis=1)

    data_clinical_sample["response"] = data_clinical_sample[response_col].map(RESPONSE_VALUE_MAP)
    return data_clinical_sample


def get_dataset_dir(lair, dataset_class, ds_name):
    ds = dataset_class(name=ds_name)
    filepaths = lair.get_dataset_filepaths(ds)
    unpacked_file_key = next(f for f in filepaths.keys() if not f.endswith('.tar.gz'))
    filepath = filepaths[unpacked_file_key]
    return filepath / filepath.name

# def load_and_process_data(data_dir):
#     data_mutations = pd.read_csv(data_dir / "data_mutations.txt", sep="\t")
#     tmb_results = calculate_maf_tmb(data_mutations)
#     tmb_results = tmb_results.set_index("Tumor_Sample_Barcode")
#     data_clinical_sample = pd.read_csv(data_dir / "data_clinical_sample.txt", sep="\t",
#                                        index_col=0, skiprows=4)
#     data_clinical_sample = data_clinical_sample.set_index("SAMPLE_ID")
#     data_clinical_sample = data_clinical_sample["RESPONSE"]
#     df = pd.concat([tmb_results, data_clinical_sample], join="inner", axis=1)
#     df["response"] = df["RESPONSE"].map({
#         "Complete Response": "R",
#         "Partial Response": "R",
#         "Stable Disease": "NR",
#         "Progressive Disease": "NR"
#     })
#     return df

def plot_tmb_vs_response(df, ax, ds_name):
    ax.set_title(ds_name)
    sns.boxplot(data=df, x="response", y="TMB_Score", ax=ax)
    response_counts = df["response"].value_counts()
    r_values = df.loc[df["response"] == "R", "TMB_Score"].dropna()
    nr_values = df.loc[df["response"] == "NR", "TMB_Score"].dropna()
    if len(r_values) > 0 and len(nr_values) > 0:
        pvalue_text = plot_significance(ax, df, nr_values, r_values)
    else:
        pvalue_text = "p=N/A"
    counts_text = "\n".join([f"{k}: {v}" for k, v in response_counts.items()])
    counts_text = f"{counts_text}\n{pvalue_text}"
    ax.text(0.95, 0.95, counts_text, transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


def plot_significance(ax, df, nr_values, r_values) -> str:
    stat, pvalue = mannwhitneyu(r_values, nr_values, alternative="two-sided")
    sig_marker = convert_to_sig_marker(pvalue)
    pvalue_text = f"p={pvalue:.3g} {sig_marker}"
    # Draw significance bracket between the two boxes
    y_max = df["TMB_Score"].max()
    y_min = df["TMB_Score"].min()
    bracket_y = y_max + (y_max - y_min) * 0.05
    bracket_h = (y_max - y_min) * 0.02
    ax.plot([0, 0, 1, 1],
            [bracket_y, bracket_y + bracket_h, bracket_y + bracket_h, bracket_y],
            lw=1.5, c="black")
    ax.text(0.5, bracket_y + bracket_h, pvalue_text,
            ha="center", va="bottom", fontsize=9,
            color="red" if pvalue < 0.05 else "black")
    return pvalue_text


def convert_to_sig_marker(pvalue) -> str:
    if pvalue < 0.001:
        sig_marker = "***"
    elif pvalue < 0.01:
        sig_marker = "**"
    elif pvalue < 0.05:
        sig_marker = "*"
    else:
        sig_marker = "ns"
    return sig_marker


def plot_tmb_iatlas(lair):
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    iatlas_dataset_names = [s for s in dataset_class.get_dataset_name() if "iAtlas" in s]
    fig, axes = plt.subplots(2, len(iatlas_dataset_names) // 2, figsize=(20, 10))
    axes = np.ndarray.flatten(axes)
    for ds_name, ax in zip(iatlas_dataset_names, axes, strict=True):
        ax.set_title(ds_name)
        try:
            data_dir = get_dataset_dir(lair, dataset_class, ds_name)
            df_clinical, _ = load_and_process_data(data_dir)
            plot_tmb_vs_response(df_clinical, ax, ds_name)
        except FileNotFoundError:
            print(ds_name, "\tFile not found")
            continue
        except KeyError:
            continue
    return fig


def predict_response(df_clinical, data_mrna_seq_expression, bagaev_signature, ax, ds_name, use_features: str = "both", n_components: int = 2, balance_classes: bool = True, classifier_type: str = "rf"):
    # Cast clinical index to string to prevent numeric/string type mismatches
    df_clinical = df_clinical.copy()
    df_clinical.index = df_clinical.index.astype(str)

    # 1. Align samples
    common_samples = df_clinical.index.intersection(data_mrna_seq_expression.columns)
    
    # 2. Filter for samples that have valid response and features
    cols_to_check = ["response"]
    if use_features in ("tmb", "both"):
        cols_to_check.append("TMB_Score")
    
    df_clinical_valid = df_clinical.loc[common_samples].dropna(subset=cols_to_check)
    valid_samples = df_clinical_valid.index
    
    if len(valid_samples) < 5:
        print(f"Skipping {ds_name} - too few samples: {len(valid_samples)}")
        ax.text(0.5, 0.5, f"Too few samples (n={len(valid_samples)})", ha="center", va="center")
        return None, None, None, None
    
    y = df_clinical_valid["response"]
    if len(np.unique(y)) < 2:
        print(f"Skipping {ds_name} - single class present in response")
        ax.text(0.5, 0.5, "Single class in response", ha="center", va="center")
        return None, None, None, None
        
    # 3. Compute Bagaev signature scores if needed
    if use_features in ("signatures", "both"):
        expr_df = data_mrna_seq_expression.loc[:, valid_samples].T
        sig_scores = ssgsea_formula(expr_df, bagaev_signature)
        
    # 4. Construct features based on the use_features selection
    if use_features == "tmb":
        features = df_clinical_valid[["TMB_Score"]]
    elif use_features == "signatures":
        features = sig_scores
    elif use_features == "both":
        features = pd.concat([df_clinical_valid["TMB_Score"], sig_scores], axis=1)
    else:
        raise ValueError(f"Invalid use_features choice: {use_features}")
        
    features.columns = features.columns.astype(str)
    
    # 5. Define ML pipeline incorporating StandardScaler and PCA on signatures
    class_weight = "balanced" if balance_classes else None
    if classifier_type == "rf":
        model = RandomForestClassifier(n_estimators=100, random_state=42, class_weight=class_weight)
    elif classifier_type == "svm":
        model = CalibratedClassifierCV(SVC(class_weight=class_weight, random_state=42), ensemble=False)
    elif classifier_type == "lr":
        model = LogisticRegression(random_state=42, class_weight=class_weight, max_iter=1000)
    else:
        raise ValueError(f"Invalid classifier_type: {classifier_type}")

    if use_features == "tmb":
        # No PCA applicable for single feature
        clf = make_pipeline(
            StandardScaler(),
            model
        )
    elif use_features == "signatures":
        # PCA on all signatures
        clf = make_pipeline(
            StandardScaler(),
            PCA(n_components=n_components),
            model
        )
    elif use_features == "both":
        # PCA on signatures only, standard scaling on TMB
        signature_cols = [col for col in features.columns if col != "TMB_Score"]
        preprocessor = ColumnTransformer(
            transformers=[
                ("sig_pca", make_pipeline(StandardScaler(), PCA(n_components=n_components)), signature_cols),
                ("tmb_scale", StandardScaler(), ["TMB_Score"])
            ]
        )
        clf = make_pipeline(
            preprocessor,
            model
        )
    
    # 6. Stratified Cross-Validation to get predicted probabilities for "R"
    y_encoded = (y == "R").astype(int)
    cv = StratifiedKFold(n_splits=min(5, len(valid_samples)), shuffle=True, random_state=42)
    y_proba = cross_val_predict(clf, features, y_encoded, cv=cv, method="predict_proba")[:, 1]
    
    # 6b. Calculate metrics
    y_pred = (y_proba >= 0.5).astype(int)
    accuracy = accuracy_score(y_encoded, y_pred)
    auc = roc_auc_score(y_encoded, y_proba)
    
    # 7. Plot ROC curve
    fpr, tpr, _ = roc_curve(y_encoded, y_proba)
    
    ax.plot(fpr, tpr, label=f"AUC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.5)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"{ds_name} (n={len(valid_samples)})")
    ax.legend(loc="lower right")
    
    # 8. Fit final model on all data
    clf.fit(features, y_encoded)
    
    return accuracy, auc, y_proba, clf


def predict_response_from_mrna_and_tmb_iatlas(lair, use_features: str = "both", n_components: int = 2, balance_classes: bool = True, classifier_type: str = "rf"):
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    iatlas_dataset_names = [s for s in dataset_class.get_dataset_name() if "iAtlas" in s]
    
    # Load Bagaev signatures once
    ds_sig = Signature()
    lair.safe_derive(ds_sig)
    filepaths = lair.get_dataset_filepaths(ds_sig)
    bagaev_signature = read_gene_sets(filepaths["gene_signatures.gmt"])
    
    fig, axes = plt.subplots(2, len(iatlas_dataset_names) // 2, figsize=(20, 10))
    axes = np.ndarray.flatten(axes)
    dataset_results = []
    for ds_name, ax in zip(iatlas_dataset_names, axes, strict=True):
        ax.set_title(ds_name)
        try:
            data_dir = get_dataset_dir(lair, dataset_class, ds_name)
            df_clinical, data_mrna_seq_expression = load_and_process_data(data_dir)
            accuracy, auc, y_pred, clf = predict_response(
                df_clinical, data_mrna_seq_expression, bagaev_signature, ax, ds_name, use_features=use_features, n_components=n_components, balance_classes=balance_classes, classifier_type=classifier_type
            )
            if accuracy is not None:
                dataset_results.append({
                    "dataset": ds_name,
                    "accuracy": accuracy,
                    "auc": auc
                })
        except FileNotFoundError:
            print(ds_name, "\tFile not found")
            continue
        except KeyError:
            continue

    return fig, dataset_results


def run_and_save_pipeline(lair_path, output_dir, use_features, n_components, balance_classes, classifier_type):
    import datalair
    import matplotlib.pyplot as plt
    lair = datalair.Lair(lair_path)
    
    if use_features == "tmb":
        filename = f"predict_response_iatlas_tmb_{classifier_type}.svg"
    else:
        filename = f"predict_response_iatlas_{use_features}_pca{n_components}_balance{balance_classes}_{classifier_type}.svg"
        
    print(f"Running: features={use_features}, PCA={n_components}, balance={balance_classes}, classifier={classifier_type}")
    
    fig, dataset_results = predict_response_from_mrna_and_tmb_iatlas(
        lair, 
        use_features=use_features, 
        n_components=n_components, 
        balance_classes=balance_classes, 
        classifier_type=classifier_type
    )
    fig.tight_layout()
    plot_output_dir = output_dir / "iAtlas-TMB"
    plot_output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_output_dir / filename)
    plt.close(fig)
    
    # Enrich results with parameters
    for res in dataset_results:
        res.update({
            "use_features": use_features,
            "n_components": n_components if use_features != "tmb" else None,
            "balance_classes": balance_classes,
            "classifier_type": classifier_type
        })
    return dataset_results


def load_and_combine_cohorts(lair, use_features: str = "both"):
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    iatlas_dataset_names = [s for s in dataset_class.get_dataset_name() if "iAtlas" in s]
    
    # Load Bagaev signatures once
    ds_sig = Signature()
    lair.safe_derive(ds_sig)
    filepaths = lair.get_dataset_filepaths(ds_sig)
    bagaev_signature = read_gene_sets(filepaths["gene_signatures.gmt"])
    
    cohort_features_list = []
    cohort_y_list = []
    cohort_names_list = []
    
    for ds_name in iatlas_dataset_names:
        try:
            data_dir = get_dataset_dir(lair, dataset_class, ds_name)
            df_clinical, data_mrna_seq_expression = load_and_process_data(data_dir)
            
            # Align samples
            df_clinical = df_clinical.copy()
            df_clinical.index = df_clinical.index.astype(str)
            common_samples = df_clinical.index.intersection(data_mrna_seq_expression.columns)
            
            cols_to_check = ["response"]
            if use_features in ("tmb", "both"):
                cols_to_check.append("TMB_Score")
                
            df_clinical_valid = df_clinical.loc[common_samples].dropna(subset=cols_to_check)
            valid_samples = df_clinical_valid.index
            
            if len(valid_samples) < 5:
                continue
                
            y = df_clinical_valid["response"]
            if len(np.unique(y)) < 2:
                continue
                
            if use_features in ("signatures", "both"):
                expr_df = data_mrna_seq_expression.loc[:, valid_samples].T
                sig_scores = ssgsea_formula(expr_df, bagaev_signature)
                
            if use_features == "tmb":
                features = df_clinical_valid[["TMB_Score"]].copy()
            elif use_features == "signatures":
                features = sig_scores.copy()
            elif use_features == "both":
                features = pd.concat([df_clinical_valid["TMB_Score"], sig_scores], axis=1)
                
            features.columns = features.columns.astype(str)
            
            cohort_features_list.append(features)
            cohort_y_list.append(y)
            cohort_names_list.extend([ds_name] * len(valid_samples))
            
        except FileNotFoundError:
            print(ds_name, "\tFile not found")
            continue
        except KeyError:
            continue
            
    if not cohort_features_list:
        raise ValueError("No cohorts successfully loaded and processed.")
        
    X_all = pd.concat(cohort_features_list, axis=0).reset_index(drop=True)
    y_all = pd.concat(cohort_y_list, axis=0).reset_index(drop=True)
    cohorts_all = np.array(cohort_names_list)
    
    return X_all, y_all, cohorts_all


def standardize_cohort_wise(X, cohorts):
    X_std = X.copy()
    for cohort in np.unique(cohorts):
        mask = (cohorts == cohort)
        cohort_data = X.loc[mask]
        mean = cohort_data.mean(axis=0)
        std = cohort_data.std(axis=0)
        std[std == 0] = 1.0
        X_std.loc[mask] = (cohort_data - mean) / std
    return X_std


def predict_response_combined(X_all, y_all, cohorts_all, ax, use_features: str = "both", n_components: int = 2, balance_classes: bool = True, classifier_type: str = "rf", eval_mode: str = "combined_cv", cohort_standardization: bool = True):
    # Apply cohort-wise standardization if enabled
    if cohort_standardization:
        X_all = standardize_cohort_wise(X_all, cohorts_all)

    # 1. Define base model
    class_weight = "balanced" if balance_classes else None
    if classifier_type == "rf":
        model = RandomForestClassifier(n_estimators=100, random_state=42, class_weight=class_weight)
    elif classifier_type == "svm":
        model = CalibratedClassifierCV(SVC(class_weight=class_weight, random_state=42), ensemble=False)
    elif classifier_type == "lr":
        model = LogisticRegression(random_state=42, class_weight=class_weight, max_iter=1000)
    else:
        raise ValueError(f"Invalid classifier_type: {classifier_type}")

    # 2. Define preprocessing pipeline
    if use_features == "tmb":
        clf = make_pipeline(StandardScaler(), model)
    elif use_features == "signatures":
        clf = make_pipeline(StandardScaler(), PCA(n_components=n_components), model)
    elif use_features == "both":
        signature_cols = [col for col in X_all.columns if col != "TMB_Score"]
        preprocessor = ColumnTransformer(
            transformers=[
                ("sig_pca", make_pipeline(StandardScaler(), PCA(n_components=n_components)), signature_cols),
                ("tmb_scale", StandardScaler(), ["TMB_Score"])
            ]
        )
        clf = make_pipeline(preprocessor, model)

    y_encoded_all = (y_all == "R").astype(int)
    y_proba_all = np.zeros(len(y_all))

    # 3. Perform Cross-Validation based on eval_mode
    if eval_mode == "combined_cv":
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        y_proba_all = cross_val_predict(clf, X_all, y_encoded_all, cv=cv, method="predict_proba")[:, 1]
    elif eval_mode == "loco_cv":
        unique_cohorts = np.unique(cohorts_all)
        for target_cohort in unique_cohorts:
            test_mask = (cohorts_all == target_cohort)
            train_mask = ~test_mask
            
            clf.fit(X_all.iloc[train_mask], y_encoded_all.iloc[train_mask])
            y_proba_all[test_mask] = clf.predict_proba(X_all.iloc[test_mask])[:, 1]
    else:
        raise ValueError(f"Invalid eval_mode: {eval_mode}")

    # 4. Compute metrics per cohort and plot ROC curves
    unique_cohorts = np.unique(cohorts_all)
    cohort_results = []
    
    cmap = plt.colormaps["tab10"]
    
    for idx, cohort in enumerate(unique_cohorts):
        mask = (cohorts_all == cohort)
        y_true_cohort = y_encoded_all.iloc[mask]
        y_proba_cohort = y_proba_all[mask]
        
        y_pred_cohort = (y_proba_cohort >= 0.5).astype(int)
        accuracy = accuracy_score(y_true_cohort, y_pred_cohort)
        
        # In case a cohort has only one class, handle ROC AUC exception
        if len(np.unique(y_true_cohort)) < 2:
            auc = np.nan
            print(f"Cohort {cohort} has single class in CV fold evaluation. AUC set to NaN.")
        else:
            auc = roc_auc_score(y_true_cohort, y_proba_cohort)
            
        cohort_results.append({
            "dataset": cohort,
            "accuracy": accuracy,
            "auc": auc
        })
        
        # Plot ROC curve for this cohort on the shared axis
        if not np.isnan(auc):
            fpr, tpr, _ = roc_curve(y_true_cohort, y_proba_cohort)
            ax.plot(fpr, tpr, label=f"{cohort} (AUC={auc:.3f})", color=cmap(idx % 10))
            
    ax.plot([0, 1], [0, 1], "k--", alpha=0.5)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"Combined Training ({eval_mode})")
    ax.legend(loc="lower right", fontsize="small")
    
    # Fit final model on all data
    clf.fit(X_all, y_encoded_all)
    
    return cohort_results, clf


def run_and_save_combined_pipeline(lair_path, output_dir, use_features, n_components, balance_classes, classifier_type, eval_mode, cohort_standardization):
    import datalair
    import matplotlib.pyplot as plt
    lair = datalair.Lair(lair_path)
    
    # Load and combine cohort datasets
    X_all, y_all, cohorts_all = load_and_combine_cohorts(lair, use_features=use_features)
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    cohort_results, clf = predict_response_combined(
        X_all, y_all, cohorts_all, ax,
        use_features=use_features,
        n_components=n_components,
        balance_classes=balance_classes,
        classifier_type=classifier_type,
        eval_mode=eval_mode,
        cohort_standardization=cohort_standardization
    )
    
    if use_features == "tmb":
        filename = f"predict_response_iatlas_combined_{eval_mode}_tmb_std{cohort_standardization}_{classifier_type}.svg"
    else:
        filename = f"predict_response_iatlas_combined_{eval_mode}_{use_features}_pca{n_components}_balance{balance_classes}_std{cohort_standardization}_{classifier_type}.svg"
        
    print(f"Running Combined: eval_mode={eval_mode}, std={cohort_standardization}, features={use_features}, PCA={n_components}, balance={balance_classes}, classifier={classifier_type}")
    
    fig.tight_layout()
    plot_output_dir = output_dir / "iAtlas-TMB"
    plot_output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_output_dir / filename)
    plt.close(fig)
    
    # Enrich results with parameters
    for res in cohort_results:
        res.update({
            "use_features": use_features,
            "n_components": n_components if use_features != "tmb" else None,
            "balance_classes": balance_classes,
            "classifier_type": classifier_type,
            "eval_mode": eval_mode,
            "cohort_standardization": cohort_standardization
        })
        
    return cohort_results


def main():
    lair = datalair.Lair("/storage/halu/lair")
    lair.assert_ok_satus()
    output_dir = Path(".") / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    fig = plot_tmb_iatlas(lair)
    fig.tight_layout()
    plot_output_dir = output_dir / "iAtlas-TMB"
    plot_output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_output_dir / "tmb_iatlas.svg")
    plt.close(fig)

    tasks = []
    
    # 1. TMB-only tasks
    for clf_type in ["rf", "svm", "lr"]:
        tasks.append({
            "use_features": "tmb",
            "n_components": 2,
            "balance_classes": True,
            "classifier_type": clf_type
        })

    # 2. PCA-dependent tasks
    for mode in ["both", "signatures"]:
        for n_comp, balance, clf_type in product([2, 3, 7], [True, False], ["rf", "svm", "lr"]):
            tasks.append({
                "use_features": mode,
                "n_components": n_comp,
                "balance_classes": balance,
                "classifier_type": clf_type
            })

    # Run tasks in parallel
    results_list = Parallel(n_jobs=-1)(
        delayed(run_and_save_pipeline)(
            lair_path="/storage/halu/lair",
            output_dir=output_dir,
            use_features=t["use_features"],
            n_components=t["n_components"],
            balance_classes=t["balance_classes"],
            classifier_type=t["classifier_type"]
        )
        for t in tasks
    )

    # Flatten list of lists
    flat_results = [res for sublist in results_list for res in sublist]

    # Convert to DataFrame and save to CSV
    results_df = pd.DataFrame(flat_results)
    results_df.to_csv(output_dir / "metrics_summary.csv", index=False)
    print(f"Saved metrics summary to {output_dir / 'metrics_summary.csv'}")

    # 3. Combined sweeps tasks
    combined_tasks = []
    for eval_mode in ["combined_cv", "loco_cv"]:
        for cohort_std in [True, False]:
            # TMB-only combined tasks
            for clf_type in ["rf", "svm", "lr"]:
                for balance in [True, False]:
                    combined_tasks.append({
                        "use_features": "tmb",
                        "n_components": 2,
                        "balance_classes": balance,
                        "classifier_type": clf_type,
                        "eval_mode": eval_mode,
                        "cohort_standardization": cohort_std
                    })
            # PCA-dependent combined tasks
            for mode in ["both", "signatures"]:
                for n_comp, balance, clf_type in product([2, 3, 7], [True, False], ["rf", "svm", "lr"]):
                    combined_tasks.append({
                        "use_features": mode,
                        "n_components": n_comp,
                        "balance_classes": balance,
                        "classifier_type": clf_type,
                        "eval_mode": eval_mode,
                        "cohort_standardization": cohort_std
                    })

    # Run combined tasks in parallel
    print(f"Running {len(combined_tasks)} combined cohort sweeps in parallel...")
    combined_results_list = Parallel(n_jobs=-1)(
        delayed(run_and_save_combined_pipeline)(
            lair_path="/storage/halu/lair",
            output_dir=output_dir,
            use_features=t["use_features"],
            n_components=t["n_components"],
            balance_classes=t["balance_classes"],
            classifier_type=t["classifier_type"],
            eval_mode=t["eval_mode"],
            cohort_standardization=t["cohort_standardization"]
        )
        for t in combined_tasks
    )

    # Flatten combined list of lists
    flat_combined_results = [res for sublist in combined_results_list for res in sublist]

    # Convert to DataFrame and save to CSV
    combined_df = pd.DataFrame(flat_combined_results)
    combined_df.to_csv(output_dir / "metrics_summary_combined.csv", index=False)
    print(f"Saved combined metrics summary to {output_dir / 'metrics_summary_combined.csv'}")

    import json

    # 4. Report Best Models
    report_lines = []
    report_lines.append("=========================================")
    report_lines.append("         BEST MODELS REPORT              ")
    report_lines.append("=========================================\n")

    # A. Best Individual Cohort Models
    report_lines.append("--- BEST MODELS PER COHORT (INDIVIDUAL PIPELINE) ---")
    df_sorted = results_df.sort_values(by=["auc", "accuracy"], ascending=[False, False])
    best_indiv_df = df_sorted.groupby("dataset").first().reset_index()
    
    indiv_dict = {}
    for _, row in best_indiv_df.iterrows():
        ds = row["dataset"]
        auc = row["auc"]
        acc = row["accuracy"]
        config = {
            "use_features": row["use_features"],
            "n_components": int(row["n_components"]) if pd.notna(row["n_components"]) else None,
            "balance_classes": bool(row["balance_classes"]),
            "classifier_type": row["classifier_type"],
            "auc": float(auc) if pd.notna(auc) else None,
            "accuracy": float(acc) if pd.notna(acc) else None
        }
        indiv_dict[ds] = config
        report_lines.append(
            f"Dataset: {ds}\n"
            f"  Best Classifier : {config['classifier_type'].upper()} (features: {config['use_features']}, PCA: {config['n_components']}, balance: {config['balance_classes']})\n"
            f"  Performance     : AUC = {config['auc']:.4f}, Accuracy = {config['accuracy']:.4f}\n"
        )

    # B. Best Combined Cohort Models
    report_lines.append("--- BEST GLOBAL MODELS (COMBINED PIPELINE) ---")
    config_cols = ["use_features", "n_components", "balance_classes", "classifier_type", "eval_mode", "cohort_standardization"]
    mean_combined = combined_df.groupby(config_cols)[["auc", "accuracy"]].mean().reset_index()
    mean_combined_sorted = mean_combined.sort_values(by=["auc", "accuracy"], ascending=[False, False])
    
    combined_dict = {}
    for eval_mode in ["combined_cv", "loco_cv"]:
        eval_df = mean_combined_sorted[mean_combined_sorted["eval_mode"] == eval_mode]
        if not eval_df.empty:
            row = eval_df.iloc[0]
            auc = row["auc"]
            acc = row["accuracy"]
            config = {
                "use_features": row["use_features"],
                "n_components": int(row["n_components"]) if pd.notna(row["n_components"]) else None,
                "balance_classes": bool(row["balance_classes"]),
                "classifier_type": row["classifier_type"],
                "cohort_standardization": bool(row["cohort_standardization"]),
                "mean_auc": float(auc) if pd.notna(auc) else None,
                "mean_accuracy": float(acc) if pd.notna(acc) else None
            }
            combined_dict[eval_mode] = config
            report_lines.append(
                f"Evaluation Mode: {eval_mode}\n"
                f"  Best Classifier : {config['classifier_type'].upper()} (features: {config['use_features']}, PCA: {config['n_components']}, balance: {config['balance_classes']}, cohort_std: {config['cohort_standardization']})\n"
                f"  Mean Performance: Mean AUC = {config['mean_auc']:.4f}, Mean Accuracy = {config['mean_accuracy']:.4f}\n"
            )
            
    report_text = "\n".join(report_lines)
    print(report_text)
    
    with open(output_dir / "best_models_report.txt", "w") as f:
        f.write(report_text)
    print(f"Saved best models text report to {output_dir / 'best_models_report.txt'}")
    
    with open(output_dir / "best_model_individual.json", "w") as f:
        json.dump(indiv_dict, f, indent=4)
    print(f"Saved best individual models JSON to {output_dir / 'best_model_individual.json'}")
    
    with open(output_dir / "best_model_combined.json", "w") as f:
        json.dump(combined_dict, f, indent=4)
    print(f"Saved best combined models JSON to {output_dir / 'best_model_combined.json'}")


if __name__ == "__main__":
    main()
