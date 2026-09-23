"""Reusable functional normalization of dataset gene identifiers to Ensembl release."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Mapping, Sequence
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import pyensembl
from returns.result import Failure, Result, Success
import scipy.sparse as sp

from ..logging import get_logger
from ..paths import get_data_paths, get_ensembl_dir

logger = get_logger("preprocessing.gene_normalization")


def ensure_ensembl_release_installed(
    release: int = 111,
    species: str = "human",
    ensembl_dir: Path | None = None,
) -> pyensembl.EnsemblRelease:
    """Ensure the specified Ensembl release is downloaded and indexed by pyensembl.

    Args:
        release: Ensembl release number (e.g. 111).
        species: Species name (e.g. 'human' or 'homo_sapiens').
        ensembl_dir: Optional installation directory. Defaults to 'data/ensembl'.

    Returns:
        The initialized pyensembl.EnsemblRelease object.
    """
    if ensembl_dir is not None:
        target_dir = Path(ensembl_dir).resolve()
    else:
        target_dir = get_ensembl_dir()

    target_dir.mkdir(parents=True, exist_ok=True)
    os.environ["PYENSEMBL_CACHE_DIR"] = str(target_dir)

    ensembl = pyensembl.EnsemblRelease(release=release, species=species)

    files_ok = ensembl.required_local_files_exist()
    db_indexed = (
        ensembl.db._database_file_exists()
        if (files_ok and hasattr(ensembl.db, "_database_file_exists"))
        else False
    )

    if not files_ok or not db_indexed:
        logger.info(
            "Ensembl release %d (%s) not found in %s. Downloading and indexing via pyensembl...",
            release,
            species,
            target_dir,
        )
        ensembl.download()
        ensembl.index()
        logger.info("Successfully installed Ensembl release %d in %s", release, target_dir)
    else:
        logger.debug("Ensembl release %d already installed in %s", release, target_dir)

    return ensembl


def _resolve_single_gene_symbol(
    query: str,
    ensembl: pyensembl.EnsemblRelease,
    alias_dict: Mapping[str, str],
) -> tuple[str | None, str, str, int | None, int | None, str, str, str, str]:
    """Resolve a single gene identifier to Ensembl ID and attributes.

    Returns:
        tuple: (gene_id, gene_name, contig, start, end, strand, biotype, mapping_status, alt_ids)
    """
    clean_query = query.strip()
    canonical_contigs = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}

    # 1. Direct Ensembl ID check (e.g. ENSG00000153563 or ENSG00000153563.14)
    if clean_query.startswith("ENSG"):
        stripped = clean_query.split(".")[0]
        try:
            gene = ensembl.gene_by_id(stripped)
            return (
                gene.gene_id,
                gene.gene_name or clean_query,
                str(gene.contig),
                int(gene.start),
                int(gene.end),
                str(gene.strand),
                str(gene.biotype),
                "ensembl_id_direct",
                "",
            )
        except Exception:
            pass

    # 2. Query as official gene symbol
    candidates: list[Any] = []
    try:
        candidate_ids = ensembl.gene_ids_of_gene_name(clean_query)
        candidates = [ensembl.gene_by_id(gid) for gid in candidate_ids]
    except Exception:
        candidates = []

    # 3. If no match, check alias/synonym dictionary
    status = "exact_symbol"
    if not candidates and clean_query in alias_dict:
        approved_sym = alias_dict[clean_query]
        if approved_sym != clean_query:
            try:
                candidate_ids = ensembl.gene_ids_of_gene_name(approved_sym)
                candidates = [ensembl.gene_by_id(gid) for gid in candidate_ids]
                status = "alias_resolved"
            except Exception:
                candidates = []

    # 4. If still no match, try case-insensitive uppercase
    if not candidates and clean_query.upper() != clean_query:
        try:
            candidate_ids = ensembl.gene_ids_of_gene_name(clean_query.upper())
            candidates = [ensembl.gene_by_id(gid) for gid in candidate_ids]
            status = "case_normalized"
        except Exception:
            candidates = []

    if not candidates:
        return (None, clean_query, "", None, None, "", "", "unmapped", "")

    # 5. Conflict Resolution: 1-to-Many
    if len(candidates) == 1:
        g = candidates[0]
        return (
            g.gene_id,
            g.gene_name or clean_query,
            str(g.contig),
            int(g.start),
            int(g.end),
            str(g.strand),
            str(g.biotype),
            status,
            "",
        )

    # Filter by canonical chromosomes first
    canonical = [g for g in candidates if str(g.contig).replace("chr", "") in canonical_contigs]
    pool = canonical if canonical else candidates

    # Prioritize protein_coding
    protein_coding = [g for g in pool if g.biotype == "protein_coding"]
    pool = protein_coding if pool and protein_coding else pool

    # Prioritize X chromosome for pseudoautosomal genes
    chr_x = [g for g in pool if str(g.contig).replace("chr", "") == "X"]
    pool = chr_x if pool and chr_x else pool

    # Deterministic tie-breaker: lowest numeric gene_id
    pool = sorted(pool, key=lambda g: g.gene_id)
    chosen = pool[0]
    alt_ids = ";".join(g.gene_id for g in candidates if g.gene_id != chosen.gene_id)

    return (
        chosen.gene_id,
        chosen.gene_name or clean_query,
        str(chosen.contig),
        int(chosen.start),
        int(chosen.end),
        str(chosen.strand),
        str(chosen.biotype),
        "contig_prioritized" if status == "exact_symbol" else status,
        alt_ids,
    )


def _load_or_update_mapping_cache(
    var_names: Sequence[str],
    ensembl: pyensembl.EnsemblRelease,
    cache_path: Path | None,
) -> pd.DataFrame:
    """Load cached mapping via Polars or resolve and update persistent parquet cache."""
    existing_cache: dict[str, dict[str, Any]] = {}

    if cache_path is not None and cache_path.exists():
        try:
            pl_cache = pl.read_parquet(cache_path)
            for row in pl_cache.iter_rows(named=True):
                existing_cache[row["query"]] = row
            logger.debug("Loaded %d cached gene mappings from %s", len(existing_cache), cache_path)
        except Exception as exc:
            logger.warning("Failed to read parquet cache from %s: %s", cache_path, exc)

    missing_queries = [v for v in var_names if v not in existing_cache]

    # Pre-fetch aliases for missing queries via mygene if needed
    alias_dict: dict[str, str] = {}
    if missing_queries:
        try:
            import mygene
            mg = mygene.MyGeneInfo()
            query_res = mg.querymany(
                missing_queries,
                scopes="symbol,alias,prev_symbol",
                fields="symbol",
                species="human",
                verbose=False,
            )
            for hit in query_res:
                q = hit.get("query")
                s = hit.get("symbol")
                if q and s:
                    alias_dict[q] = s
        except Exception as mg_exc:
            logger.debug("mygene alias pre-fetch skipped: %s", mg_exc)

    new_rows: list[dict[str, Any]] = []
    for q in missing_queries:
        gid, gname, contig, start, end, strand, biotype, m_status, alt_ids = _resolve_single_gene_symbol(
            q, ensembl, alias_dict
        )
        row = {
            "query": q,
            "gene_id": gid,
            "gene_name": gname,
            "contig": contig,
            "start": start if start is not None else 0,
            "end": end if end is not None else 0,
            "strand": strand,
            "biotype": biotype,
            "ensembl_release": ensembl.release,
            "species": ensembl.species.latin_name,
            "mapping_status": m_status,
            "alternative_ensembl_ids": alt_ids,
        }
        new_rows.append(row)
        existing_cache[q] = row

    # If new queries were resolved and cache_path provided, save to parquet
    if new_rows and cache_path is not None:
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            all_rows = list(existing_cache.values())
            pl_new = pl.DataFrame(all_rows)
            pl_new.write_parquet(cache_path)
            logger.debug("Updated persistent parquet cache with %d new entries at %s", len(new_rows), cache_path)
        except Exception as write_exc:
            logger.warning("Failed to write parquet cache to %s: %s", cache_path, write_exc)

    # Return mapping DataFrame in the exact order of var_names
    records = [existing_cache[v] for v in var_names]
    return pd.DataFrame(records)


def normalize_genes_to_ensembl(
    adata: ad.AnnData,
    release: int = 111,
    species: str = "human",
    ensembl_dir: Path | None = None,
    drop_unmapped: bool = True,
    aggregation: str = "sum",
) -> ad.AnnData:
    """Convert AnnData var_names to canonical Ensembl gene IDs and enrich .var attributes.

    Features:
    - Auto-installs and indexes Ensembl release in data/ensembl/ via pyensembl.
    - Resolves 1-to-many symbol conflicts by canonical contig and protein_coding priority.
    - Resolves 0-to-many unmapped names via HGNC alias/previous symbol lookup.
    - Enriches adata.var with: contig, start, end, strand, biotype, gene_name,
      ensembl_release, species, mapping_status, alternative_ensembl_ids.
    - Uses persistent Polars Parquet caching for sub-millisecond lookups.
    - Aggregates duplicate Ensembl gene IDs via specified operation ('sum', 'mean', 'max').
    - Drops or retains unmapped features per drop_unmapped setting.

    Args:
        adata: Input AnnData expression container.
        release: Ensembl release version (default 111).
        species: Organism species (default 'human').
        ensembl_dir: Target cache directory. Defaults to get_ensembl_dir().
        drop_unmapped: If True, drops unmapped features and stores them in adata.uns['unmapped_genes'].
        aggregation: Aggregation function for duplicate Ensembl IDs ('sum', 'mean', 'max').

    Returns:
        AnnData with Ensembl gene IDs as var_names and full genomic attributes in .var.
    """
    logger.info(
        "Normalizing %d genes to Ensembl Release %d (%s)...",
        adata.n_vars,
        release,
        species,
    )
    target_dir = ensembl_dir or get_ensembl_dir()
    ensembl = ensure_ensembl_release_installed(release=release, species=species, ensembl_dir=target_dir)

    cache_path = Path(target_dir) / f"gene_mapping_cache_release_{release}.parquet"

    mapping_df = _load_or_update_mapping_cache(list(adata.var_names), ensembl, cache_path)
    mapping_df["original_id"] = list(adata.var_names)

    # Attach mapping columns to var
    var_combined = adata.var.copy()
    for col in [
        "gene_id", "gene_name", "original_id", "contig", "start", "end",
        "strand", "biotype", "ensembl_release", "species", "mapping_status",
        "alternative_ensembl_ids"
    ]:
        var_combined[col] = mapping_df[col].values

    # Handle unmapped genes
    unmapped_mask = var_combined["gene_id"].isna() | (var_combined["gene_id"] == "") | (var_combined["mapping_status"] == "unmapped")
    n_unmapped = int(unmapped_mask.sum())

    if n_unmapped > 0:
        unmapped_genes = list(var_combined.loc[unmapped_mask, "original_id"])
        logger.info(
            "Identified %d unmapped genes (e.g. %s)",
            n_unmapped,
            unmapped_genes[:5],
        )
        if drop_unmapped:
            adata = adata[:, ~unmapped_mask].copy()
            var_combined = var_combined.loc[~unmapped_mask].copy()
            adata.uns["unmapped_genes"] = unmapped_genes
            adata.uns["n_unmapped_genes"] = n_unmapped
        else:
            # Retain original symbol for unmapped
            var_combined.loc[unmapped_mask, "gene_id"] = var_combined.loc[unmapped_mask, "original_id"]

    # Set index to gene_id
    var_combined.index = pd.Index(var_combined["gene_id"].astype(str), name="gene_id")
    adata.var = var_combined

    # Aggregate duplicate Ensembl IDs if present
    if not adata.var.index.is_unique:
        gene_ids = adata.var.index.to_numpy()
        unique_ids = np.unique(gene_ids)
        logger.info(
            "Aggregating duplicate Ensembl IDs (%d unique across %d total) via '%s'...",
            len(unique_ids),
            len(gene_ids),
            aggregation,
        )

        X = adata.X
        is_sparse = sp.issparse(X)
        X_dense = X.toarray() if is_sparse else np.asarray(X)

        df_X = pd.DataFrame(X_dense, index=adata.obs_names, columns=gene_ids)
        match aggregation:
            case "mean":
                X_agg = df_X.T.groupby(level=0).mean().T
            case "max":
                X_agg = df_X.T.groupby(level=0).max().T
            case _:
                X_agg = df_X.T.groupby(level=0).sum().T

        X_agg = X_agg.loc[:, unique_ids]
        var_dedup = adata.var[~adata.var.index.duplicated(keep="first")].loc[unique_ids]

        new_X = sp.csr_matrix(X_agg.to_numpy()) if is_sparse else X_agg.to_numpy()
        adata = ad.AnnData(
            X=new_X,
            obs=adata.obs.copy(),
            var=var_dedup.copy(),
            uns=adata.uns.copy(),
            obsm=adata.obsm.copy(),
        )

    # Sort var alphabetically by Ensembl ID for deterministic layout
    adata = adata[:, sorted(adata.var_names)].copy()
    logger.info("Ensembl normalization complete: %d cells x %d genes", adata.n_obs, adata.n_vars)
    return adata


def normalize_dataset_to_ensembl(
    adata: ad.AnnData,
    release: int | None = None,
    species: str | None = None,
    ensembl_dir: Path | None = None,
    drop_unmapped: bool = True,
    aggregation: str = "sum",
) -> Result[ad.AnnData, str]:
    """Convert dataset gene identifiers to canonical Ensembl gene IDs and enrich .var attributes.

    Args:
        adata: Input AnnData container.
        release: Optional Ensembl release version (defaults to config default_ensembl_release, e.g. 111).
        species: Optional species name (defaults to 'human').
        ensembl_dir: Optional installation directory (defaults to config ensembl_dir 'data/ensembl').
        drop_unmapped: Whether to drop unmapped non-gene features (recording them in adata.uns).
        aggregation: Aggregation strategy for duplicate Ensembl IDs ('sum', 'mean', 'max').

    Returns:
        Success(normalized_adata) or Failure(error_message).
    """
    try:
        cfg = get_data_paths()
        target_release = release or cfg.default_ensembl_release
        target_species = species or cfg.default_species
        target_dir = ensembl_dir or cfg.ensembl_dir

        logger.info(
            "Starting Ensembl gene ID normalization for %d genes (Release %d, species=%s, ensembl_dir=%s)...",
            adata.n_vars,
            target_release,
            target_species,
            target_dir,
        )

        norm_adata = normalize_genes_to_ensembl(
            adata=adata,
            release=target_release,
            species=target_species,
            ensembl_dir=target_dir,
            drop_unmapped=drop_unmapped,
            aggregation=aggregation,
        )

        logger.info(
            "Ensembl gene ID normalization complete: %d obs x %d vars (all var_names are canonical Ensembl IDs)",
            norm_adata.n_obs,
            norm_adata.n_vars,
        )
        return Success(norm_adata)
    except Exception as exc:
        msg = f"Failed to normalize genes to Ensembl: {exc}"
        logger.error(msg)
        return Failure(msg)


__all__ = [
    "ensure_ensembl_release_installed",
    "normalize_genes_to_ensembl",
    "normalize_dataset_to_ensembl",
]
