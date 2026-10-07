"""Immutable data models and configuration for single-cell deconvolution references."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence
import anndata as ad
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success


class DeconvolutionReferenceConfig(BaseModel):
    """Configuration for building deconvolution references from single-cell data."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    cell_state_key: str = "cell_state"
    cell_type_key: Maybe[str] = Nothing
    gene_symbol_key: Maybe[str] = Nothing
    malignant_key: Maybe[str] = Nothing
    malignant_label: str = "Malignant"
    auto_detect_malignant: bool = True
    filter_confounding: bool = True
    min_cells_per_state: int = 15
    collinearity_threshold: float = 0.85
    pseudo_min: float = 1e-8
    expression_layer: Maybe[str] = Nothing
    normalize_multinomial: bool = True


class DeconvolutionReferenceResult(BaseModel):
    """Immutable result of single-cell deconvolution reference construction."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    phi_state: pl.DataFrame
    phi_type: pl.DataFrame
    hierarchy_table: pl.DataFrame
    mrna_scaling: pl.DataFrame
    collinearity: pl.DataFrame
    gene_names: tuple[str, ...]
    cell_states: tuple[str, ...]
    cell_types: tuple[str, ...]
    malignant_states: tuple[str, ...]
    metadata: Mapping[str, Any] = Field(default_factory=dict)

    def to_numpy(self, which: str = "state") -> tuple[np.ndarray, list[str], list[str]]:
        """Extract a 2D numpy expression matrix (K x G), state/type labels, and gene names.

        Args:
            which: 'state' for fine cell states (default) or 'type' for broad cell types.

        Returns:
            (matrix, row_labels, gene_names)
        """
        df = self.phi_state if which == "state" else self.phi_type
        id_col = "cell_state" if which == "state" else "cell_type"
        labels = df[id_col].to_list()
        genes = [c for c in df.columns if c != id_col]
        mat = df.select(genes).to_numpy().astype(np.float64)
        return mat, labels, genes

    def to_instaprism(self) -> tuple[np.ndarray, list[str]]:
        """Extract linear reference matrix and state labels for InstaPrism deconvolution.

        Returns:
            (reference_matrix, state_labels) where matrix has shape (K_states, G_genes).
        """
        mat, labels, _ = self.to_numpy(which="state")
        return mat, labels

    def to_bayesprism(self) -> Result[Any, str]:
        """Convert reference into BayesPrism S4 equivalents (RefPhi and RefTumor)."""
        from .adapters import export_to_bayesprism
        return export_to_bayesprism(self)

    def to_anndata(self) -> ad.AnnData:
        """Convert the deconvolution reference into an AnnData object."""
        mat, labels, genes = self.to_numpy(which="state")
        obs_df = self.hierarchy_table.to_pandas().set_index("cell_state").loc[labels]
        var_df = pl.DataFrame({"gene_name": genes}).to_pandas().set_index("gene_name")
        adata = ad.AnnData(X=mat, obs=obs_df, var=var_df)
        adata.uns["mrna_scaling"] = self.mrna_scaling.to_dicts()
        adata.uns["collinearity"] = self.collinearity.to_dicts()
        adata.uns["cell_types"] = list(self.cell_types)
        adata.uns["malignant_states"] = list(self.malignant_states)
        return adata

    def export_parquet(self, out_dir: Path) -> Result[Path, str]:
        """Export all reference tables to standardized Parquet files in the specified directory.

        Files written:
        - curated_reference_phi.parquet
        - reference_phi_cell_type.parquet
        - cell_hierarchy.parquet
        - mrna_scaling.parquet
        - reference_collinearity.parquet
        """
        try:
            out_dir.mkdir(parents=True, exist_ok=True)
            self.phi_state.write_parquet(out_dir / "curated_reference_phi.parquet")
            self.phi_state.write_parquet(out_dir / "reference_phi.parquet")
            self.phi_type.write_parquet(out_dir / "reference_phi_cell_type.parquet")
            self.hierarchy_table.write_parquet(out_dir / "cell_hierarchy.parquet")
            self.mrna_scaling.write_parquet(out_dir / "mrna_scaling.parquet")
            self.collinearity.write_parquet(out_dir / "reference_collinearity.parquet")
            return Success(out_dir)
        except Exception as exc:
            return Failure(f"Failed to export reference parquets to '{out_dir}': {exc}")
