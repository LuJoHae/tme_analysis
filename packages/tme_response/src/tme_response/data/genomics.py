"""Pure functional genomic data extractors for TMB, somatic mutations, and CNA."""

from __future__ import annotations

from pathlib import Path
import polars as pl
from returns.result import Failure, Result, Success

# Standard non-synonymous mutation classifications in MAF files
NON_SYNONYMOUS_CLASSES = frozenset({
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Splice_Site",
    "Translation_Start_Site",
    "Nonstop_Mutation",
})

# Inactivating / loss-of-function variants
INACTIVATING_CLASSES = frozenset({
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "Translation_Start_Site",
})

KEY_ICB_GENES = frozenset({
    "B2M", "JAK1", "JAK2", "STK11", "KEAP1", "PTEN", "TP53", "EGFR", "KRAS", "BRAF", "NRAS"
})

ESTIMATED_EXOME_MB = 38.0


def parse_maf_mutations(maf_path: Path) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Parse MAF/mutation file to calculate TMB and extract key driver alterations.

    Returns:
        Result containing tuple of:
          1. tmb_df: [sample_id, tmb_per_mb, n_nonsynonymous]
          2. driver_df: [sample_id, gene, variant_classification, is_inactivating]
    """
    if not maf_path.exists():
        return Failure(f"Mutation file not found at: {maf_path}")

    try:
        # Read MAF with comment lines skipped
        df = pl.read_csv(
            maf_path,
            separator="\t",
            comment_prefix="#",
            infer_schema_length=10000,
            ignore_errors=True,
        )

        # Standardize column names to uppercase
        col_map = {c: c.upper() for c in df.columns}
        df = df.rename(col_map)

        sample_col = next((c for c in df.columns if any(k in c for k in ("TUMOR_SAMPLE_BARCODE", "SAMPLE_ID", "SAMPLEID"))), None)
        gene_col = next((c for c in df.columns if any(k in c for k in ("HUGO_SYMBOL", "GENE"))), None)
        class_col = next((c for c in df.columns if any(k in c for k in ("VARIANT_CLASSIFICATION", "MUTATION_TYPE"))), None)

        if not sample_col or not gene_col or not class_col:
            return Failure(f"MAF file at {maf_path} missing required columns (Sample, Gene, Classification).")

        maf_clean = df.select([
            pl.col(sample_col).cast(pl.String).alias("sample_id"),
            pl.col(gene_col).cast(pl.String).alias("gene"),
            pl.col(class_col).cast(pl.String).alias("variant_classification"),
        ])

        # Filter non-synonymous mutations for TMB
        non_syn = maf_clean.filter(pl.col("variant_classification").is_in(list(NON_SYNONYMOUS_CLASSES)))

        tmb_df = (
            non_syn.group_by("sample_id")
            .agg(pl.len().alias("n_nonsynonymous"))
            .with_columns(
                (pl.col("n_nonsynonymous") / ESTIMATED_EXOME_MB).alias("tmb_per_mb")
            )
        )

        # Extract driver mutations in key genes
        driver_df = (
            maf_clean.filter(pl.col("gene").is_in(list(KEY_ICB_GENES)))
            .with_columns(
                pl.col("variant_classification")
                .is_in(list(INACTIVATING_CLASSES))
                .alias("is_inactivating")
            )
        )

        return Success((tmb_df, driver_df))
    except Exception as exc:
        return Failure(f"Failed to parse MAF mutations from {maf_path}: {exc}")


def parse_cna_file(cna_path: Path) -> Result[pl.DataFrame, str]:
    """Parse copy number alterations table (genes as rows, samples as columns).

    Returns:
        Result containing cna_df: [sample_id, cdkn2a_loss, aneuploidy_score]
    """
    if not cna_path.exists():
        return Failure(f"CNA file not found at: {cna_path}")

    try:
        df = pl.read_csv(
            cna_path,
            separator="\t",
            comment_prefix="#",
            infer_schema_length=5000,
            ignore_errors=True,
        )

        # First column is typically Hugo_Symbol / Gene
        gene_col = df.columns[0]
        sample_cols = df.columns[1:]

        # Check for CDKN2A (9p21) deletion
        cdkn2a_row = df.filter(pl.col(gene_col).str.to_uppercase() == "CDKN2A")
        has_cdkn2a = not cdkn2a_row.is_empty()

        # Reshape to long format: sample_id, gene, cna_val
        cna_long = df.unpivot(index=[gene_col], on=sample_cols, variable_name="sample_id", value_name="cna_val")
        cna_long = cna_long.with_columns(pl.col("cna_val").cast(pl.Float64, strict=False))

        # Aneuploidy proxy: fraction of altered genes (cna_val != 0)
        aneuploidy_df = (
            cna_long.group_by("sample_id")
            .agg(
                (pl.col("cna_val").abs() >= 1.0).mean().alias("aneuploidy_score")
            )
        )

        if has_cdkn2a:
            cdkn2a_vals = cdkn2a_row.unpivot(index=[gene_col], on=sample_cols, variable_name="sample_id", value_name="cna_val")
            cdkn2a_df = cdkn2a_vals.select([
                pl.col("sample_id"),
                (pl.col("cna_val").cast(pl.Float64, strict=False) <= -1.5).alias("cdkn2a_loss")
            ])
            res_df = aneuploidy_df.join(cdkn2a_df, on="sample_id", how="left")
        else:
            res_df = aneuploidy_df.with_columns(pl.lit(False).alias("cdkn2a_loss"))

        return Success(res_df)
    except Exception as exc:
        return Failure(f"Failed to parse CNA file from {cna_path}: {exc}")
