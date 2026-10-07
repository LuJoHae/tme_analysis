"""Provider loaders for the 17 Pan-Cancer Single-Cell Reference Atlas datasets.

All loaders follow functional programming principles:
- Pure functions returning Result[AnnData, str].
- Non-destructive metadata standardization with canonical columns:
  cell_id, patient, sample, organ, cancer_type, cancer_code, cell_type_author, dataset.
- Sparse CSR expression matrix construction.
- Verification and metadata tagging of expression type (raw integer counts vs normalized).
- Support for optional subsetting (by cancer_type, patient, etc.) and subsampling.
"""

from __future__ import annotations

import gzip
from pathlib import Path
import re
import shutil
import tempfile
from typing import Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import scipy.io as sio
import scipy.sparse as sp
from returns.result import Failure, Result, Success

from ..download.fetcher import download_single_file, unpack_tar
from ..download.geo import download_arrayexpress_files, download_geo_supplementary
from ..logging import get_logger
from ..preprocessing.matrix_inspection import tag_expression_metadata

logger = get_logger("providers.atlas_single_cell")


def _is_dir_empty(p: Path) -> bool:
    """Return True if path does not exist or contains no files/directories."""
    return not p.exists() or not any(p.iterdir())


def _resolve_10x_file(directory: Path, prefix: str, candidate_suffixes: Sequence[str]) -> Path | None:
    """Find the first existing candidate file matching prefix + suffix."""
    for sfx in candidate_suffixes:
        candidate = directory / f"{prefix}{sfx}"
        if candidate.is_file():
            return candidate
    return None


def _to_csr(mat: object) -> sp.csr_matrix:
    """Converts a sparse or dense matrix (genes x cells) into cells x genes CSR matrix."""
    if sp.issparse(mat):
        return mat.T.tocsr()  # type: ignore[attr-defined]
    return sp.csr_matrix(np.asarray(mat).T)


def _standardize_obs(
    adata: ad.AnnData,
    dataset: str,
    organ: str,
    cancer_type: str,
    cancer_code: str,
    patient_col: str | None = None,
    sample_col: str | None = None,
    cell_type_col: str | None = None,
) -> ad.AnnData:
    """Standardize .obs columns while preserving all existing columns purely."""
    obs = adata.obs.copy()
    obs["cell_id"] = adata.obs_names.astype(str)
    obs["dataset"] = dataset
    obs["organ"] = organ
    obs["cancer_type"] = cancer_type
    obs["cancer_code"] = cancer_code

    obs["patient"] = obs[patient_col].astype(str) if patient_col and patient_col in obs.columns else "Unknown"
    obs["sample"] = obs[sample_col].astype(str) if sample_col and sample_col in obs.columns else "Unknown"
    obs["cell_type_author"] = obs[cell_type_col].astype(str) if cell_type_col and cell_type_col in obs.columns else "Unknown"

    new_adata = adata.copy()
    new_adata.obs = obs
    new_adata.obs_names = obs["cell_id"]
    new_adata.obs_names.name = "cell_id"
    return new_adata


def _apply_subset_and_subsample(
    adata: ad.AnnData,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
    seed: int = 42,
) -> ad.AnnData:
    """Purely filters and/or downsamples an AnnData object."""
    res_adata = adata

    # 1. Apply column matching subset if provided
    if subset:
        mask = np.ones(res_adata.n_obs, dtype=bool)
        for col, val in subset.items():
            if col in res_adata.obs.columns:
                mask = mask & (res_adata.obs[col].astype(str).str.lower() == str(val).lower())
            else:
                logger.warning("Subset column '%s' not present in .obs (available: %s)", col, list(res_adata.obs.columns))
        if np.any(mask):
            res_adata = res_adata[mask].copy()
        else:
            logger.warning("Subset filter %s matched 0 cells, returning original.", subset)

    # 2. Downsample if subsample_n is specified and smaller than n_obs
    if subsample_n is not None and 0 < subsample_n < res_adata.n_obs:
        rng = np.random.default_rng(seed)
        chosen_indices = np.sort(rng.choice(res_adata.n_obs, size=subsample_n, replace=False))
        res_adata = res_adata[chosen_indices].copy()

    return res_adata


# =========================================================================
# 1. Pelka et al. 2021 (GSE178341) - Colorectal Cancer
# =========================================================================
def load_pelka_crc(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Pelka et al. 2021 Colorectal Cancer atlas (GSE178341, ~65k cells)."""
    expected = [
        "GSE178341_crc10x_full_c295v4_submit.h5",
        "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz",
        "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz",
    ]
    if auto_download:
        match download_geo_supplementary("GSE178341", raw_dir, expected_files=expected):
            case Failure(err):
                return Failure(f"Failed to fetch Pelka raw files: {err}")

    h5_path = raw_dir / "GSE178341_crc10x_full_c295v4_submit.h5"
    meta_path = raw_dir / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"
    cluster_path = raw_dir / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"

    if not h5_path.exists():
        return Failure(f"Pelka 10x H5 file not found at {h5_path}")

    try:
        import scanpy as sc
        adata = sc.read_10x_h5(h5_path)
        if meta_path.exists() and cluster_path.exists():
            meta_df = pd.read_csv(meta_path)
            clust_df = pd.read_csv(cluster_path)
            merged_meta = pd.concat([meta_df.set_index("cellID"), clust_df.set_index("sampleID")], axis=1)
            adata.obs = merged_meta.reindex(adata.obs_names)

        # Standardize features
        if "gene_ids" in adata.var.columns:
            adata.var["gene_name"] = adata.var_names
            clean_ids = [gid.split(".")[0] for gid in adata.var["gene_ids"].astype(str)]
            adata.var_names = clean_ids
            adata.var_names_make_unique()

        adata = _standardize_obs(
            adata,
            dataset="Pelka2021",
            organ="Colon",
            cancer_type="Colorectal Cancer",
            cancer_code="COAD",
            patient_col="PatientTypeID",
            sample_col="Sample",
            cell_type_col="Cluster",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Pelka et al. dataset: {exc}")


# =========================================================================
# 2. Azizi et al. 2018 (GSE114727) - Breast Cancer
# =========================================================================
def load_azizi_brca(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Azizi et al. 2018 Breast Cancer atlas (GSE114727, ~45k cells)."""
    tar_path = raw_dir / "GSE114727_RAW.tar"
    if auto_download and not tar_path.exists():
        match download_geo_supplementary("GSE114727", raw_dir, expected_files=["GSE114727_RAW.tar"]):
            case Failure(err):
                return Failure(f"Failed to download Azizi archive: {err}")

    extract_dir = raw_dir / "extracted_azizi"
    if _is_dir_empty(extract_dir) and tar_path.exists():
        match unpack_tar(tar_path, extract_dir):
            case Failure(err):
                return Failure(f"Failed to unpack Azizi archive: {err}")

    try:
        mtx_files = sorted(list(extract_dir.glob("*matrix.mtx*")))
        adatas: list[ad.AnnData] = []

        for mtx_p in mtx_files:
            prefix = mtx_p.name.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
            bc_p = _resolve_10x_file(extract_dir, prefix, ["_barcodes.tsv.gz", "_barcodes.tsv"])
            genes_p = _resolve_10x_file(extract_dir, prefix, ["_genes.tsv.gz", "_genes.tsv", "_features.tsv.gz", "_features.tsv"])

            if bc_p is None or genes_p is None:
                continue

            mat = _to_csr(sio.mmread(mtx_p))
            barcodes = [f"{prefix}_{b}" for b in pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()]
            genes_df = pd.read_csv(genes_p, header=None, sep="\t")
            ensembl_ids = genes_df[0].astype(str).tolist()
            gene_symbols = genes_df[1].astype(str).tolist() if len(genes_df.columns) > 1 else ensembl_ids

            var_df = pd.DataFrame(index=ensembl_ids)
            var_df["gene_name"] = gene_symbols

            sub_adata = ad.AnnData(X=mat.astype(np.float32), obs=pd.DataFrame(index=barcodes), var=var_df)
            parts = prefix.split("_")
            sub_adata.obs["geo_id"] = parts[0] if len(parts) > 0 else "Unknown"
            sub_adata.obs["patient"] = parts[1] if len(parts) > 1 else "Unknown"
            sub_adata.obs["tissue"] = parts[2] if len(parts) > 2 else "Breast"
            adatas.append(sub_adata)

        if not adatas:
            # Fallback to counts.csv.gz if mtx not present
            csv_files = sorted(list(extract_dir.glob("*counts.csv.gz")))
            for cp in csv_files:
                df = pd.read_csv(cp, index_col=0).fillna(0)
                sub_adata = ad.AnnData(X=sp.csr_matrix(df.values.T.astype(np.float32)), obs=pd.DataFrame(index=df.columns), var=pd.DataFrame(index=df.index))
                sub_adata.var["gene_name"] = sub_adata.var_names
                adatas.append(sub_adata)

        if not adatas:
            return Failure(f"No valid single-cell matrices extracted from {extract_dir}")

        combined = ad.concat(adatas, axis=0, join="outer")
        combined = _standardize_obs(
            combined,
            dataset="Azizi2018",
            organ="Breast",
            cancer_type="Breast Cancer",
            cancer_code="BRCA",
            patient_col="patient",
            sample_col="geo_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Azizi et al. dataset: {exc}")


# =========================================================================
# 3. Qian et al. 2020 (E-MTAB-8107) - Pan-Cancer
# =========================================================================
def load_qian_pancancer(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Qian et al. 2020 Pan-Cancer blueprint (E-MTAB-8107, ~200k cells)."""
    if auto_download:
        match download_arrayexpress_files("E-MTAB-8107", raw_dir):
            case Failure(err):
                # If full listing fails, check if files exist locally
                pass

    counts_files = sorted(list(raw_dir.glob("*.counts.csv*")))
    sdrf_files = list(raw_dir.glob("*sdrf*.txt"))

    if not counts_files:
        return Failure(f"Qian counts files not found in {raw_dir}")

    try:
        sdrf_df = pd.read_csv(sdrf_files[0], sep="\t") if sdrf_files else pd.DataFrame()
        adatas: list[ad.AnnData] = []

        for cp in counts_files:
            df = pl.read_csv(cp)
            genes = df.columns[0]
            gene_names = df[genes].to_list()
            cell_ids = [c for c in df.columns if c != genes]
            mat = sp.csr_matrix(df.select(cell_ids).to_numpy().T.astype(np.float32))

            sub_adata = ad.AnnData(
                X=mat,
                obs=pd.DataFrame(index=cell_ids),
                var=pd.DataFrame(index=gene_names),
            )
            sub_adata.var["gene_name"] = gene_names
            source_name = cp.name.split(".")[0]
            sub_adata.obs["sample"] = source_name

            if not sdrf_df.empty and "Source Name" in sdrf_df.columns:
                match_rows = sdrf_df[sdrf_df["Source Name"].astype(str) == source_name]
                if len(match_rows) > 0:
                    for k in match_rows.columns:
                        sub_adata.obs[k] = match_rows[k].iloc[0]

            adatas.append(sub_adata)

        combined = ad.concat(adatas, axis=0, join="outer")
        combined = _standardize_obs(
            combined,
            dataset="Qian2020",
            organ="Multi-organ",
            cancer_type="Pan-Cancer",
            cancer_code="Pan-Cancer",
            patient_col="Characteristics[individual]" if "Characteristics[individual]" in combined.obs.columns else None,
            sample_col="sample",
            cell_type_col="Characteristics[cell type]" if "Characteristics[cell type]" in combined.obs.columns else None,
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Qian et al. dataset: {exc}")


# =========================================================================
# 4. Cheng et al. 2021 (GSE154763) - Pan-Cancer T-Cells
# =========================================================================
def load_cheng_pancancer(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Cheng et al. 2021 Pan-Cancer T-cell atlas (GSE154763, ~390k cells)."""
    if auto_download:
        match download_geo_supplementary("GSE154763", raw_dir):
            case Failure(err):
                pass

    expr_files = sorted(list(raw_dir.glob("*normalized_expression.csv.gz")))
    if not expr_files:
        return Failure(f"Cheng expression files not found in {raw_dir}")

    try:
        valid_pairs: list[tuple[Path, Path, str]] = []
        for ef in expr_files:
            meta_f = Path(str(ef).replace("normalized_expression.csv.gz", "metadata.csv.gz"))
            if meta_f.exists():
                cancer_tag = ef.name.split("_")[1] if len(ef.name.split("_")) > 1 else "Pan"
                valid_pairs.append((ef, meta_f, cancer_tag))

        if not valid_pairs:
            return Failure(f"No valid expression/metadata pairs found in {raw_dir}")

        # Determine unified gene space from CSV headers without loading full files
        gene_set: set[str] = set()
        for ef, _, _ in valid_pairs:
            header_df = pd.read_csv(ef, index_col=0, nrows=0)
            gene_set.update(header_df.columns)

        all_genes = sorted(list(gene_set))
        gene_to_col = {g: i for i, g in enumerate(all_genes)}
        n_all_genes = len(all_genes)

        obs_dfs: list[pd.DataFrame] = []
        csr_matrices: list[sp.csr_matrix] = []
        import gc

        for ef, meta_f, cancer_tag in valid_pairs:
            meta_df = pd.read_csv(meta_f, index_col=0)
            meta_df["cancer_type_study"] = cancer_tag
            obs_dfs.append(meta_df)

            # Read expression CSV
            df_expr = pd.read_csv(ef, index_col=0)
            col_indices = np.array([gene_to_col[g] for g in df_expr.columns], dtype=np.int32)
            raw_csr = sp.csr_matrix(df_expr.to_numpy(dtype=np.float32))
            del df_expr
            gc.collect()

            # Align sparse columns directly to unified gene space
            aligned_indices = col_indices[raw_csr.indices]
            aligned_csr = sp.csr_matrix(
                (raw_csr.data, aligned_indices, raw_csr.indptr),
                shape=(raw_csr.shape[0], n_all_genes),
                dtype=np.float32,
            )
            del raw_csr
            csr_matrices.append(aligned_csr)

        combined_X = sp.vstack(csr_matrices, format="csr")
        del csr_matrices
        combined_obs = pd.concat(obs_dfs, axis=0)
        del obs_dfs
        gc.collect()

        var_df = pd.DataFrame(index=all_genes)
        var_df["gene_name"] = all_genes
        combined = ad.AnnData(X=combined_X, obs=combined_obs, var=var_df)

        combined = _standardize_obs(
            combined,
            dataset="Cheng2021",
            organ="Multi-organ",
            cancer_type="Pan-Cancer T-Cells",
            cancer_code="Pan-Cancer",
            patient_col="patient" if "patient" in combined.obs.columns else None,
            sample_col="sample" if "sample" in combined.obs.columns else None,
            cell_type_col="cell_type" if "cell_type" in combined.obs.columns else "cluster",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Cheng et al. dataset: {exc}")


# =========================================================================
# 5. Leader et al. 2021 (GSE154826) - Non-Small Cell Lung Cancer
# =========================================================================
def load_leader_nsclc(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Leader et al. 2021 NSCLC atlas (GSE154826, ~360k cells)."""
    meta_p = raw_dir / "leader_cell_metadata.csv"
    annots_p = raw_dir / "leader_annots_list.csv"

    if auto_download:
        match download_geo_supplementary("GSE154826", raw_dir):
            case Failure(err):
                pass
        if not meta_p.exists():
            download_single_file(
                "https://raw.githubusercontent.com/effiken/Leader_et_al/master/input_tables/cell_metadata.csv",
                meta_p,
            )
        if not annots_p.exists():
            download_single_file(
                "https://raw.githubusercontent.com/effiken/Leader_et_al/master/input_tables/annots_list.csv",
                annots_p,
            )

    tar_files = sorted([f for f in raw_dir.glob("*.tar*") if f.is_file() and not f.name.endswith(".part")])
    if not tar_files:
        return Failure(f"Leader tar archives not found in {raw_dir}")

    extract_dir = raw_dir / "extracted_leader"
    if _is_dir_empty(extract_dir):
        for tf in tar_files:
            unpack_tar(tf, extract_dir)

    try:
        # Load author cell annotations if available
        has_cell_meta = meta_p.is_file() and annots_p.is_file()
        df_cells = pd.read_csv(meta_p) if has_cell_meta else None
        df_annots = pd.read_csv(annots_p) if has_cell_meta else None

        cluster_map: dict[int, dict[str, str]] = {}
        if df_annots is not None:
            cluster_map = df_annots.set_index("cluster")[["lineage", "sub_lineage"]].to_dict("index")

        def _resolve_cell_type(cluster_id: int) -> str:
            info = cluster_map.get(cluster_id, {})
            sub = info.get("sub_lineage")
            lin = info.get("lineage")
            if pd.notna(sub) and str(sub).strip():
                return str(sub).strip()
            if pd.notna(lin) and str(lin).strip():
                return str(lin).strip()
            return f"Cluster_{cluster_id}"

        if df_cells is not None:
            df_cells["cell_type_author"] = [
                _resolve_cell_type(int(c)) for c in df_cells["cluster_ID"]
            ]
            df_cells["barcode"] = [
                cid.split("_", 1)[1] if "_" in str(cid) else str(cid)
                for cid in df_cells["cell_ID"]
            ]

        # Load sample metadata if available
        samp_p = raw_dir / "GSE154826_sample_annots.csv.gz"
        if not samp_p.exists():
            samp_p = raw_dir / "GSE154826_sample_annots.csv"
        df_samp = pd.read_csv(samp_p) if samp_p.is_file() else None

        batch_to_samples: dict[int, list[dict[str, object]]] = {}
        if df_samp is not None:
            for _, r in df_samp.iterrows():
                bid = int(r["amp_batch_ID"])
                batch_to_samples.setdefault(bid, []).append(r.to_dict())

        mtx_files = sorted(list(extract_dir.glob("*matrix.mtx*")))
        if not mtx_files:
            return Failure(f"No matrix files found in {extract_dir}")

        adatas: list[ad.AnnData] = []
        for mp in mtx_files:
            prefix = mp.name.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
            bc_p = _resolve_10x_file(extract_dir, prefix, ["_barcodes.tsv.gz", "_barcodes.tsv"])
            feat_p = _resolve_10x_file(extract_dir, prefix, ["_features.tsv.gz", "_features.tsv", "_genes.tsv.gz", "_genes.tsv"])
            if bc_p is None or feat_p is None:
                continue

            bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
            feats = pd.read_csv(feat_p, header=None, sep="\t")
            var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
            var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values

            mat = sio.mmread(mp)  # genes x cells

            batch_prefix_num: int | None = None
            try:
                batch_prefix_num = int(prefix.split("_")[0])
            except ValueError:
                batch_prefix_num = None

            samp_records = batch_to_samples.get(batch_prefix_num, []) if batch_prefix_num is not None else []

            if df_cells is not None and samp_records:
                samp_ids = [r["sample_ID"] for r in samp_records]
                sub_cells = df_cells[df_cells["sample_ID"].isin(samp_ids)]
                if sub_cells.empty:
                    continue

                bc_to_idx = {b: i for i, b in enumerate(bcs)}
                matched_indices: list[int] = []
                matched_cell_ids: list[str] = []
                matched_types: list[str] = []
                matched_patients: list[str] = []
                matched_samples: list[str] = []

                sample_map = {r["sample_ID"]: r for r in samp_records}

                for _, crow in sub_cells.iterrows():
                    b = str(crow["barcode"])
                    if b in bc_to_idx:
                        matched_indices.append(bc_to_idx[b])
                        matched_cell_ids.append(f"{prefix}_{b}")
                        matched_types.append(str(crow["cell_type_author"]))
                        s_info = sample_map.get(crow["sample_ID"], {})
                        matched_patients.append(str(s_info.get("patient_ID", "Unknown")))
                        matched_samples.append(str(s_info.get("sample_ID", prefix)))

                if not matched_indices:
                    continue

                sub_mat = mat.tocsc()[:, matched_indices].T.tocsr()
                obs_df = pd.DataFrame(
                    {
                        "patient": matched_patients,
                        "sample": matched_samples,
                        "cell_type_author": matched_types,
                        "batch": prefix,
                    },
                    index=matched_cell_ids,
                )
                sub_adata = ad.AnnData(X=sub_mat.astype(np.float32), obs=obs_df, var=var_df)
                adatas.append(sub_adata)
            else:
                # Fallback: transpose and filter out zero-count cells to avoid millions of empty droplets
                sub_mat = _to_csr(mat)
                counts = np.asarray(sub_mat.sum(axis=1)).ravel()
                valid_mask = counts > 0
                if not np.any(valid_mask):
                    continue
                sub_mat = sub_mat[valid_mask]
                valid_bcs = [f"{prefix}_{b}" for b in np.array(bcs)[valid_mask]]

                parts = prefix.split("_")
                patient = parts[1] if len(parts) > 1 else "Unknown"
                obs_df = pd.DataFrame(
                    {
                        "patient": patient,
                        "sample": prefix,
                        "batch": prefix,
                        "cell_type_author": "Unknown",
                    },
                    index=valid_bcs,
                )
                sub_adata = ad.AnnData(X=sub_mat.astype(np.float32), obs=obs_df, var=var_df)
                adatas.append(sub_adata)

        if not adatas:
            return Failure(f"No valid AnnData objects could be created from {extract_dir}")

        combined = ad.concat(adatas, axis=0, join="outer")
        combined = _standardize_obs(
            combined,
            dataset="Leader2021",
            organ="Lung",
            cancer_type="Non-Small Cell Lung",
            cancer_code="NSCLC",
            patient_col="patient",
            sample_col="sample",
            cell_type_col="cell_type_author",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Leader et al. dataset: {exc}")


# =========================================================================
# 6. Kim et al. 2020 (GSE131907) - Lung Adenocarcinoma
# =========================================================================
def load_kim_luad(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Kim et al. 2020 Lung Adenocarcinoma atlas (GSE131907, ~40k cells)."""
    expected = [
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    ]
    if auto_download:
        match download_geo_supplementary("GSE131907", raw_dir, expected_files=expected):
            case Failure(err):
                pass

    mat_p = raw_dir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    annot_p = raw_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"

    if not mat_p.exists():
        return Failure(f"Kim UMI matrix not found at {mat_p}")

    try:
        df = pl.read_csv(mat_p, separator="\t")
        gene_col = df.columns[0]
        genes = df[gene_col].to_list()
        cell_ids = [c for c in df.columns if c != gene_col]

        mat = sp.csr_matrix(df.select(cell_ids).to_numpy().T.astype(np.float32))
        adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=cell_ids), var=pd.DataFrame(index=genes))
        adata.var["gene_name"] = genes

        if annot_p.exists():
            ann_df = pd.read_csv(annot_p, sep="\t")
            if "Index" in ann_df.columns:
                ann_df = ann_df.set_index("Index")
            elif "Cell" in ann_df.columns:
                ann_df = ann_df.set_index("Cell")
            adata.obs = ann_df.reindex(adata.obs_names)

        adata = _standardize_obs(
            adata,
            dataset="Kim2020",
            organ="Lung",
            cancer_type="Lung Adenocarcinoma",
            cancer_code="LUAD",
            patient_col="Sample_Origin" if "Sample_Origin" in adata.obs.columns else "Sample",
            sample_col="Sample" if "Sample" in adata.obs.columns else None,
            cell_type_col="Cell_type" if "Cell_type" in adata.obs.columns else None,
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Kim et al. dataset: {exc}")


def _download_becker_files(extract_dir: Path) -> Result[Path, str]:
    """Download Becker et al. 2022 RNA expression files directly from GEO GSM suppl."""
    extract_dir.mkdir(parents=True, exist_ok=True)
    fl_path = extract_dir.parent / "filelist.txt"
    if not fl_path.exists():
        fl_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE201nnn/GSE201349/suppl/filelist.txt"
        match download_single_file(fl_url, fl_path):
            case Failure(err):
                return Failure(f"Failed to download Becker filelist: {err}")

    rna_files: list[str] = []
    with open(fl_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                name = parts[1]
                if name.endswith("_matrix.mtx.gz") or name.endswith("_barcodes.tsv.gz") or name.endswith("_features.tsv.gz"):
                    rna_files.append(name)

    logger.info("Checking %d Becker RNA matrix files from NCBI GEO...", len(rna_files))
    for fname in rna_files:
        dest = extract_dir / fname
        if dest.exists() and dest.stat().st_size > 0:
            continue
        gsm = fname.split("_")[0]
        bucket = f"{gsm[:7]}nnn"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{bucket}/{gsm}/suppl/{fname}"
        download_single_file(url, dest)

    return Success(extract_dir)


# =========================================================================
# 7. Becker et al. 2022 (GSE201349) - Colorectal Cancer
# =========================================================================
def load_becker_coad(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Becker et al. 2022 Colorectal Cancer continuum (GSE201349, ~30k cells)."""
    tar_path = raw_dir / "GSE201349_RAW.tar"
    extract_dir = raw_dir / "extracted_becker"

    if _is_dir_empty(extract_dir):
        if tar_path.exists():
            match unpack_tar(tar_path, extract_dir):
                case Failure(err):
                    pass
        if _is_dir_empty(extract_dir) and auto_download:
            _download_becker_files(extract_dir)

    try:
        files = sorted(list(extract_dir.glob("*matrix.mtx*")))
        adatas: list[ad.AnnData] = []
        for mtx_p in files:
            prefix = mtx_p.name.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
            bc_p = _resolve_10x_file(extract_dir, prefix, ["_barcodes.tsv.gz", "_barcodes.tsv"])
            feat_p = _resolve_10x_file(extract_dir, prefix, ["_features.tsv.gz", "_features.tsv", "_genes.tsv.gz", "_genes.tsv"])
            if bc_p is None or feat_p is None:
                continue

            mat = _to_csr(sio.mmread(mtx_p))
            bcs = [f"{prefix}_{b}" for b in pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()]
            feats = pd.read_csv(feat_p, header=None, sep="\t")

            var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
            var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values

            sub_adata = ad.AnnData(X=mat.astype(np.float32), obs=pd.DataFrame(index=bcs), var=var_df)
            sub_adata.obs["sample"] = prefix
            sub_adata.obs["geo_id"] = prefix.split("_")[0]
            adatas.append(sub_adata)

        if not adatas:
            return Failure(f"No valid single-cell matrices extracted from {extract_dir}")

        combined = ad.concat(adatas, axis=0, join="outer")
        combined = _standardize_obs(
            combined,
            dataset="Becker2022",
            organ="Colon",
            cancer_type="Colorectal Cancer",
            cancer_code="COAD",
            sample_col="sample",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Becker et al. dataset: {exc}")


# =========================================================================
# 8. Khaliq et al. 2022 (GSE200997) - Colorectal Cancer
# =========================================================================
def load_khaliq_cc(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Khaliq et al. 2022 Colorectal Cancer classification (GSE200997, ~25k cells)."""
    expected = [
        "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz",
        "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz",
    ]
    if auto_download:
        match download_geo_supplementary("GSE200997", raw_dir, expected_files=expected):
            case Failure(err):
                pass

    mat_p = raw_dir / "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"
    annot_p = raw_dir / "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"

    if not mat_p.exists():
        return Failure(f"Khaliq raw count matrix not found at {mat_p}")

    try:
        df = pl.read_csv(mat_p)
        gene_col = df.columns[0]
        genes = df[gene_col].to_list()
        cell_ids = [c for c in df.columns if c != gene_col]

        mat = sp.csr_matrix(df.select(cell_ids).to_numpy().T.astype(np.float32))
        adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=cell_ids), var=pd.DataFrame(index=genes))
        adata.var["gene_name"] = genes

        if annot_p.exists():
            ann_df = pd.read_csv(annot_p, index_col=0)
            adata.obs = ann_df.reindex(adata.obs_names)

        adata = _standardize_obs(
            adata,
            dataset="Khaliq2022",
            organ="Colon",
            cancer_type="Colorectal Cancer",
            cancer_code="COAD",
            patient_col="patient" if "patient" in adata.obs.columns else "samples",
            cell_type_col="CellType" if "CellType" in adata.obs.columns else None,
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Khaliq et al. dataset: {exc}")


# =========================================================================
# 9. Borcherding et al. 2021 (GSE121638) - Clear Cell Renal Cell Carcinoma
# =========================================================================
def load_borcherding_ccrcc(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Borcherding et al. 2021 ccRCC immune environment (GSE121638, ~25k cells)."""
    tar_path = raw_dir / "GSE121638_RAW.tar"
    if auto_download and not tar_path.exists():
        match download_geo_supplementary("GSE121638", raw_dir, expected_files=["GSE121638_RAW.tar"]):
            case Failure(err):
                pass

    extract_dir = raw_dir / "extracted_borcherding"
    if _is_dir_empty(extract_dir) and tar_path.exists():
        match unpack_tar(tar_path, extract_dir):
            case Failure(err):
                pass

    try:
        mtx_files = sorted(list(extract_dir.glob("*matrix.mtx*")))
        adatas: list[ad.AnnData] = []
        for mp in mtx_files:
            prefix = mp.name.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
            gene_p = _resolve_10x_file(extract_dir, prefix, ["_genes.tsv.gz", "_genes.tsv", "_features.tsv.gz", "_features.tsv"])
            bc_p = _resolve_10x_file(extract_dir, prefix, ["_barcodes.tsv.gz", "_barcodes.tsv"])
            if gene_p is None or bc_p is None:
                continue

            mat = _to_csr(sio.mmread(mp))
            bcs = [f"{prefix}_{b}" for b in pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()]
            genes_df = pd.read_csv(gene_p, header=None, sep="\t")

            var_df = pd.DataFrame(index=genes_df[0].astype(str).tolist())
            var_df["gene_name"] = (genes_df[1] if len(genes_df.columns) > 1 else genes_df[0]).astype(str).values

            sub_adata = ad.AnnData(X=mat.astype(np.float32), obs=pd.DataFrame(index=bcs), var=var_df)
            sub_adata.obs["sample"] = prefix
            parts = prefix.split("_")
            sub_adata.obs["geo_id"] = parts[0] if len(parts) > 0 else "Unknown"
            sub_adata.obs["patient"] = parts[1] if len(parts) > 1 else "Unknown"
            adatas.append(sub_adata)

        if not adatas:
            return Failure(f"No valid single-cell matrices extracted from {extract_dir}")

        combined = ad.concat(adatas, axis=0, join="outer")
        combined = _standardize_obs(
            combined,
            dataset="Borcherding2021",
            organ="Kidney",
            cancer_type="Clear Cell Renal Cell",
            cancer_code="KIRC",
            patient_col="patient",
            sample_col="sample",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Borcherding et al. dataset: {exc}")


# =========================================================================
# 10. Sharma et al. 2020 (GSE156625) - Hepatocellular Carcinoma
# =========================================================================
def load_sharma_hcc(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Sharma et al. 2020 HCC onco-fetal atlas (GSE156625, ~15k cells)."""
    expected = [
        "GSE156625_HCCmatrix.mtx.gz",
        "GSE156625_HCCbarcodes.tsv.gz",
        "GSE156625_HCCgenes.tsv.gz",
    ]
    if auto_download:
        match download_geo_supplementary("GSE156625", raw_dir, expected_files=expected):
            case Failure(err):
                pass

    mat_p = raw_dir / "GSE156625_HCCmatrix.mtx.gz"
    bc_p = raw_dir / "GSE156625_HCCbarcodes.tsv.gz"
    gene_p = raw_dir / "GSE156625_HCCgenes.tsv.gz"

    if not mat_p.exists():
        return Failure(f"Sharma HCC matrix not found at {mat_p}")

    try:
        mat = _to_csr(sio.mmread(mat_p))
        bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
        genes_df = pd.read_csv(gene_p, header=None, sep="\t")

        var_df = pd.DataFrame(index=genes_df[0].astype(str).tolist())
        var_df["gene_name"] = (genes_df[1] if len(genes_df.columns) > 1 else genes_df[0]).astype(str).values

        adata = ad.AnnData(X=mat.astype(np.float32), obs=pd.DataFrame(index=bcs), var=var_df)
        adata = _standardize_obs(
            adata,
            dataset="Sharma2020",
            organ="Liver",
            cancer_type="Hepatocellular Carcinoma",
            cancer_code="LIHC",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Sharma et al. dataset: {exc}")


# =========================================================================
# 11. Lu et al. 2022 (GSE149614) - Hepatocellular Carcinoma
# =========================================================================
def load_lu_hcc(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Lu et al. 2022 Liver Cancer multicellular atlas (GSE149614, ~18k cells)."""
    expected = [
        "GSE149614_HCC.scRNAseq.S71915.count.txt.gz",
        "GSE149614_HCC.metadata.updated.txt.gz",
    ]
    if auto_download:
        match download_geo_supplementary("GSE149614", raw_dir, expected_files=expected):
            case Failure(err):
                pass

    count_p = raw_dir / "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
    meta_p = raw_dir / "GSE149614_HCC.metadata.updated.txt.gz"

    if not count_p.exists():
        return Failure(f"Lu count matrix not found at {count_p}")

    try:
        # Read gzipped text count matrix safely
        with gzip.open(count_p, "rt", encoding="utf-8") as f_in:
            first_line = f_in.readline().rstrip("\r\n").split("\t")
            # If first column is missing gene name, prepend header
            has_gene_col = first_line[0].lower() in ("gene", "gene_id", "symbol")

        if has_gene_col:
            df = pl.read_csv(count_p, separator="\t")
        else:
            # Handle ragged missing header on column 0
            df_pd = pd.read_csv(count_p, sep="\t", index_col=0)
            genes = df_pd.index.astype(str).tolist()
            mat = sp.csr_matrix(df_pd.values.T.astype(np.float32))
            adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=df_pd.columns), var=pd.DataFrame(index=genes))
            adata.var["gene_name"] = genes
            if meta_p.exists():
                m_df = pd.read_csv(meta_p, sep="\t", index_col=0)
                adata.obs = m_df.reindex(adata.obs_names)
            adata = _standardize_obs(
                adata,
                dataset="Lu2022",
                organ="Liver",
                cancer_type="Hepatocellular Carcinoma",
                cancer_code="LIHC",
                patient_col="patient" if "patient" in adata.obs.columns else None,
                cell_type_col="celltype" if "celltype" in adata.obs.columns else None,
            )
            adata = tag_expression_metadata(adata)
            return Success(_apply_subset_and_subsample(adata, subset, subsample_n))

        gene_col = df.columns[0]
        genes = df[gene_col].to_list()
        cell_ids = [c for c in df.columns if c != gene_col]
        mat = sp.csr_matrix(df.select(cell_ids).to_numpy().T.astype(np.float32))
        adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=cell_ids), var=pd.DataFrame(index=genes))
        adata.var["gene_name"] = genes

        if meta_p.exists():
            m_df = pd.read_csv(meta_p, sep="\t", index_col=0)
            adata.obs = m_df.reindex(adata.obs_names)

        adata = _standardize_obs(
            adata,
            dataset="Lu2022",
            organ="Liver",
            cancer_type="Hepatocellular Carcinoma",
            cancer_code="LIHC",
            patient_col="patient" if "patient" in adata.obs.columns else None,
            cell_type_col="celltype" if "celltype" in adata.obs.columns else None,
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Lu et al. dataset: {exc}")


# =========================================================================
# 12. Pu et al. 2021 (GSE184362) - Papillary Thyroid Carcinoma
# =========================================================================
def load_pu_ptc(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Pu et al. 2021 Papillary Thyroid Carcinoma atlas (GSE184362, ~20k cells)."""
    tar_path = raw_dir / "GSE184362_RAW.tar"
    if auto_download and not tar_path.exists():
        match download_geo_supplementary("GSE184362", raw_dir, expected_files=["GSE184362_RAW.tar"]):
            case Failure(err):
                pass

    extract_dir = raw_dir / "extracted_pu"
    if _is_dir_empty(extract_dir) and tar_path.exists():
        match unpack_tar(tar_path, extract_dir):
            case Failure(err):
                pass

    try:
        mtx_files = sorted(list(extract_dir.glob("*matrix.mtx*")))
        adatas: list[ad.AnnData] = []
        for mp in mtx_files:
            prefix = mp.name.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
            bc_p = _resolve_10x_file(extract_dir, prefix, ["_barcodes.tsv.gz", "_barcodes.tsv"])
            feat_p = _resolve_10x_file(extract_dir, prefix, ["_features.tsv.gz", "_features.tsv", "_genes.tsv.gz", "_genes.tsv"])
            if bc_p is None or feat_p is None:
                continue

            mat = _to_csr(sio.mmread(mp))
            bcs = [f"{prefix}_{b}" for b in pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()]
            feats = pd.read_csv(feat_p, header=None, sep="\t")

            var_df = pd.DataFrame(index=feats[0].astype(str).tolist())
            var_df["gene_name"] = (feats[1] if len(feats.columns) > 1 else feats[0]).astype(str).values

            sub_adata = ad.AnnData(X=mat.astype(np.float32), obs=pd.DataFrame(index=bcs), var=var_df)
            parts = prefix.split("_")
            sub_adata.obs["geo_id"] = parts[0] if len(parts) > 0 else "Unknown"
            sub_adata.obs["patient"] = parts[1] if len(parts) > 1 else "Unknown"
            sub_adata.obs["biopsy_site"] = parts[2] if len(parts) > 2 else "Unknown"
            sub_adata.obs["sample"] = prefix
            adatas.append(sub_adata)

        if not adatas:
            return Failure(f"No valid single-cell matrices extracted from {extract_dir}")

        combined = ad.concat(adatas, axis=0, join="outer")
        combined = _standardize_obs(
            combined,
            dataset="Pu2021",
            organ="Thyroid",
            cancer_type="Papillary Thyroid Carcinoma",
            cancer_code="THCA",
            patient_col="patient",
            sample_col="sample",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Pu et al. dataset: {exc}")


# =========================================================================
# 13. Durante et al. 2020 (GSE139829) - Uveal Melanoma
# =========================================================================
def load_durante_uvm(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Durante et al. 2020 Uveal Melanoma atlas (GSE139829, ~10k cells)."""
    tar_path = raw_dir / "GSE139829_RAW.tar"
    if auto_download and not tar_path.exists():
        match download_geo_supplementary("GSE139829", raw_dir, expected_files=["GSE139829_RAW.tar"]):
            case Failure(err):
                pass

    extract_dir = raw_dir / "extracted_durante"
    if _is_dir_empty(extract_dir) and tar_path.exists():
        match unpack_tar(tar_path, extract_dir):
            case Failure(err):
                pass

    try:
        mtx_files = sorted(list(extract_dir.glob("*matrix.mtx*")))
        adatas: list[ad.AnnData] = []
        for mp in mtx_files:
            prefix = mp.name.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
            bc_p = _resolve_10x_file(extract_dir, prefix, ["_barcodes.tsv.gz", "_barcodes.tsv"])
            gene_p = _resolve_10x_file(extract_dir, prefix, ["_genes.tsv.gz", "_genes.tsv", "_features.tsv.gz", "_features.tsv"])
            if bc_p is None or gene_p is None:
                continue

            mat = _to_csr(sio.mmread(mp))
            bcs = [f"{prefix}_{b}" for b in pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()]
            genes = pd.read_csv(gene_p, header=None, sep="\t")

            var_df = pd.DataFrame(index=genes[0].astype(str).tolist())
            var_df["gene_name"] = (genes[1] if len(genes.columns) > 1 else genes[0]).astype(str).values

            sub_adata = ad.AnnData(X=mat.astype(np.float32), obs=pd.DataFrame(index=bcs), var=var_df)
            sub_adata.obs["sample"] = prefix
            sub_adata.obs["geo_id"] = prefix.split("_")[0]
            adatas.append(sub_adata)

        if not adatas:
            return Failure(f"No valid single-cell matrices extracted from {extract_dir}")

        combined = ad.concat(adatas, axis=0, join="outer")
        combined = _standardize_obs(
            combined,
            dataset="Durante2020",
            organ="Eye",
            cancer_type="Uveal Melanoma",
            cancer_code="UVM",
            sample_col="sample",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Durante et al. dataset: {exc}")


# =========================================================================
# 14. Biermann et al. 2022 (GSE200218) - Melanoma Brain Metastasis
# =========================================================================
def load_biermann_brainmet(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Biermann et al. 2022 Melanoma Brain Met atlas (GSE200218, ~12k cells)."""
    expected = [
        "GSE200218_sc_sn_counts.mtx.gz",
        "GSE200218_sc_sn_metadata.csv.gz",
        "GSE200218_sc_sn_gene_names.csv.gz",
    ]
    if auto_download:
        match download_geo_supplementary("GSE200218", raw_dir, expected_files=expected):
            case Failure(err):
                pass

    cnt_p = raw_dir / "GSE200218_sc_sn_counts.mtx.gz"
    meta_p = raw_dir / "GSE200218_sc_sn_metadata.csv.gz"
    gene_p = raw_dir / "GSE200218_sc_sn_gene_names.csv.gz"

    if not cnt_p.exists():
        return Failure(f"Biermann counts matrix not found at {cnt_p}")

    try:
        mat = _to_csr(sio.mmread(cnt_p))
        genes_df = pd.read_csv(gene_p, index_col=0)
        # Exclude Excel date conversion artifacts
        mask = (genes_df.index != "1-Mar") & (genes_df.index != "2-Mar")
        genes_clean = genes_df.loc[mask]
        mat_clean = mat[:, mask]

        meta_df = pd.read_csv(meta_p, index_col=0) if meta_p.exists() else pd.DataFrame(index=range(mat_clean.shape[0]))
        var_df = pd.DataFrame(index=genes_clean.index.astype(str))
        var_df["gene_name"] = genes_clean.index.astype(str)

        adata = ad.AnnData(X=mat_clean.astype(np.float32), obs=meta_df, var=var_df)
        adata = _standardize_obs(
            adata,
            dataset="Biermann2022",
            organ="Brain",
            cancer_type="Melanoma Brain Metastasis",
            cancer_code="SKCM",
            patient_col="patient" if "patient" in adata.obs.columns else None,
            cell_type_col="cell_type_main" if "cell_type_main" in adata.obs.columns else "cell_type",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Biermann et al. dataset: {exc}")


# =========================================================================
# 15. Vazquez et al. 2022 (GSE180661) - Ovarian Cancer
# =========================================================================
def load_vazquez_ov(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Vazquez et al. 2022 Ovarian Cancer atlas (GSE180661, ~15k cells)."""
    expected = [
        "GSE180661_matrix.h5",
        "GSE180661_GEO_cells.tsv.gz",
    ]
    if auto_download:
        match download_geo_supplementary("GSE180661", raw_dir, expected_files=expected):
            case Failure(err):
                pass

    h5_p = raw_dir / "GSE180661_matrix.h5"
    meta_p = raw_dir / "GSE180661_GEO_cells.tsv.gz"

    if not h5_p.exists():
        return Failure(f"Vazquez matrix H5 not found at {h5_p}")

    try:
        try:
            adata = ad.read_h5ad(h5_p)
        except Exception:
            import scanpy as sc
            adata = sc.read_10x_h5(h5_p)

        adata.obs_names_make_unique()
        if meta_p.exists():
            meta_df = pd.read_csv(meta_p, sep="\t", index_col=0)
            adata.obs = meta_df.reindex(adata.obs_names)

        adata = _standardize_obs(
            adata,
            dataset="Vazquez2022",
            organ="Ovary",
            cancer_type="Ovarian Cancer",
            cancer_code="OV",
            patient_col="patient_id" if "patient_id" in adata.obs.columns else "patient",
            cell_type_col="cell_type" if "cell_type" in adata.obs.columns else None,
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Vazquez et al. dataset: {exc}")


# =========================================================================
# 16. Zhang et al. 2021 (GSE169246) - Triple-Negative Breast Cancer
# =========================================================================
def load_zhang_tnbc(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Zhang et al. 2021 TNBC immunotherapy atlas (GSE169246, ~18k cells)."""
    expected = [
        "GSE169246_TNBC_RNA.counts.mtx.gz",
        "GSE169246_TNBC_RNA.barcode.tsv.gz",
        "GSE169246_TNBC_RNA.feature.tsv.gz",
    ]
    if auto_download:
        match download_geo_supplementary("GSE169246", raw_dir, expected_files=expected):
            case Failure(err):
                pass

    cnt_p = raw_dir / "GSE169246_TNBC_RNA.counts.mtx.gz"
    bc_p = raw_dir / "GSE169246_TNBC_RNA.barcode.tsv.gz"
    feat_p = raw_dir / "GSE169246_TNBC_RNA.feature.tsv.gz"

    if not cnt_p.exists():
        return Failure(f"Zhang TNBC counts matrix not found at {cnt_p}")

    try:
        mat = _to_csr(sio.mmread(cnt_p))
        bcs = pd.read_csv(bc_p, header=None, sep="\t")[0].astype(str).tolist()
        feats_df = pd.read_csv(feat_p, header=None, sep="\t", index_col=0)
        genes = feats_df.index.astype(str).tolist()

        var_df = pd.DataFrame(index=genes)
        var_df["gene_name"] = feats_df[1].astype(str).values if len(feats_df.columns) > 0 else genes

        adata = ad.AnnData(X=mat.astype(np.float32), obs=pd.DataFrame(index=bcs), var=var_df)
        adata = _standardize_obs(
            adata,
            dataset="Zhang2021",
            organ="Breast",
            cancer_type="Triple-Negative Breast",
            cancer_code="BRCA",
        )
        adata = tag_expression_metadata(adata)
        return Success(_apply_subset_and_subsample(adata, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Zhang 2021 dataset: {exc}")


# =========================================================================
# 17. Zhang et al. 2022 (GSE215120) - Pan-Cancer Myeloid
# =========================================================================
def load_zhang_myeloid(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Zhang et al. 2022 Pan-Cancer Myeloid atlas (GSE215120, ~50k cells)."""
    tar_path = raw_dir / "GSE215120_RAW.tar"
    if auto_download and not tar_path.exists():
        match download_geo_supplementary("GSE215120", raw_dir, expected_files=["GSE215120_RAW.tar"]):
            case Failure(err):
                pass

    extract_dir = raw_dir / "extracted_zhang2022"
    if _is_dir_empty(extract_dir) and tar_path.exists():
        match unpack_tar(tar_path, extract_dir):
            case Failure(err):
                pass

    try:
        import scanpy as sc
        h5_files = sorted(list(extract_dir.glob("*.h5*")))
        if not h5_files:
            return Failure(f"No 10x H5 files found in {extract_dir}")

        adatas: list[ad.AnnData] = []
        for hp in h5_files:
            sub_adata = sc.read_10x_h5(hp)
            if "gene_ids" in sub_adata.var.columns:
                sub_adata.var["gene_name"] = sub_adata.var_names
                sub_adata.var_names = sub_adata.var["gene_ids"].astype(str)

            sub_adata.obs_names = [f"{hp.stem}_{b}" for b in sub_adata.obs_names]
            parts = hp.stem.split("_")
            sub_adata.obs["geo_id"] = parts[0] if len(parts) > 0 else "Unknown"
            sub_adata.obs["patient_id"] = parts[1] if len(parts) > 1 else "Unknown"
            adatas.append(sub_adata)

        combined = ad.concat(adatas, axis=0, join="outer")
        combined = _standardize_obs(
            combined,
            dataset="Zhang2022",
            organ="Multi-organ",
            cancer_type="Pan-Cancer Myeloid",
            cancer_code="Pan-Cancer",
            patient_col="patient_id",
            sample_col="geo_id",
        )
        combined = tag_expression_metadata(combined)
        return Success(_apply_subset_and_subsample(combined, subset, subsample_n))
    except Exception as exc:
        return Failure(f"Failed to parse Zhang 2022 dataset: {exc}")
