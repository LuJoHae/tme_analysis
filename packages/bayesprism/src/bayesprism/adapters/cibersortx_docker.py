"""
Adapter for CIBERSORTx Docker Container Deconvolution.

Provides functional integration of the official containerized CIBERSORTx suite
(Newman et al., Nature Biotechnology 2019):
- Packaging of bulk mixture and reference signature matrices into Stanford CIBERSORT TSV format
- Declarative generation of 'docker run' commands for 'cibersortx/fractions'
- Configurable execution with B-mode and S-mode batch correction, quantile normalization, and permutations
- Subprocess execution with Docker volume mounting and credential resolution
- Robust Polars-based parsing of inferred cellular fractions, empirical p-values, RMSE, and correlation
- Dependency injection / mock execution support for robust offline testing in CI/CD

Adheres to strict functional programming standards, Pydantic (frozen=True), and Returns (Result/Maybe).
"""

from __future__ import annotations

import io
import os
from pathlib import Path
import subprocess
from typing import Callable, Sequence
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Some, Nothing
from returns.result import Result, Success, Failure


class CibersortXDockerConfig(BaseModel):
    """Configuration options for CIBERSORTx Docker container."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    username: str = ""
    token: str = ""
    image_name: str = "cibersortx/fractions"
    rmbatch_bmode: bool = False
    rmbatch_smode: bool = False
    qn: bool = False
    n_perm: int = 0
    absolute: bool = False
    abs_method: str = "sigscore"
    docker_binary: str = "docker"
    timeout_seconds: int = 3600


class CibersortXResult(BaseModel):
    """Container holding inferred cell proportions and diagnostics from CIBERSORTx."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    proportions: pl.DataFrame
    cell_types: tuple[str, ...]
    sample_names: tuple[str, ...]
    p_values: Maybe[pl.Series] = Field(default=Nothing)
    correlation: pl.Series
    rmse: pl.Series
    absolute_scores: Maybe[pl.Series] = Field(default=Nothing)
    raw_output_text: Maybe[str] = Field(default=Nothing)


# =========================================================================
# Pure Core Functions
# =========================================================================


def _parse_env_file(env_path: Path) -> dict[str, str]:
    """Purely parse key-value pairs from a .env file without external dependencies."""
    if not env_path.is_file():
        return {}
    res: dict[str, str] = {}
    try:
        for line in env_path.read_text(encoding="utf-8").splitlines():
            sline = line.strip()
            if not sline or sline.startswith("#") or "=" not in sline:
                continue
            k, v = sline.split("=", 1)
            k = k.strip()
            v = v.strip().strip("'\"")
            if k:
                res[k] = v
    except Exception:
        pass
    return res


def resolve_cibersortx_credentials(
    config: CibersortXDockerConfig,
    env_path: Path = Path(".env"),
) -> Result[tuple[str, str], str]:
    """
    Purely resolve CIBERSORTx username and token from configuration, environment, or .env.

    Checked sources in order:
    1. Direct fields in CibersortXDockerConfig
    2. Environment variables: CIBERSORTX_USERNAME, CIBERSORTX_TOKEN
    3. Workspace .env file
    """
    user = (
        config.username.strip()
        if config.username.strip()
        else os.environ.get("CIBERSORTX_USERNAME", "").strip()
    )
    tok = (
        config.token.strip()
        if config.token.strip()
        else os.environ.get("CIBERSORTX_TOKEN", "").strip()
    )

    if not user or not tok:
        env_dict = _parse_env_file(env_path)
        if not user:
            user = env_dict.get("CIBERSORTX_USERNAME", "").strip()
        if not tok:
            tok = env_dict.get("CIBERSORTX_TOKEN", "").strip()

    if not user:
        return Failure(
            "CIBERSORTx username is required. Specify it in CibersortXDockerConfig, set the CIBERSORTX_USERNAME environment variable, or add CIBERSORTX_USERNAME to .env."
        )
    if not tok:
        return Failure(
            "CIBERSORTx token is required. Specify it in CibersortXDockerConfig, set the CIBERSORTX_TOKEN environment variable, or add CIBERSORTX_TOKEN to .env."
        )

    return Success((user, tok))


def format_cibersort_matrix(
    matrix: np.ndarray,
    row_names: Sequence[str],
    col_names: Sequence[str],
) -> Result[str, str]:
    """
    Format a 2D matrix (genes x samples or genes x cell types) into CIBERSORT tab-delimited text.

    First column header is 'GeneSymbol'.
    """
    if matrix.ndim != 2:
        return Failure(f"Expected 2D matrix, got shape {matrix.shape}")

    n_rows, n_cols = matrix.shape
    if n_rows != len(row_names):
        return Failure(
            f"Row count mismatch: matrix has {n_rows} rows, got {len(row_names)} row names."
        )
    if n_cols != len(col_names):
        return Failure(
            f"Column count mismatch: matrix has {n_cols} columns, got {len(col_names)} column names."
        )
    if n_rows == 0 or n_cols == 0:
        return Failure("Cannot format empty matrix.")

    # Deduplicate row names if needed or verify uniqueness
    if len(set(row_names)) != len(row_names):
        return Failure(
            "Gene names must be unique for CIBERSORTx tab-delimited formatting."
        )

    buf = io.StringIO()
    buf.write("GeneSymbol\t" + "\t".join(col_names) + "\n")

    for i in range(n_rows):
        row_vals = "\t".join(f"{matrix[i, j]:.6g}" for j in range(n_cols))
        buf.write(f"{row_names[i]}\t{row_vals}\n")

    return Success(buf.getvalue())


def build_cibersortx_docker_command(
    config: CibersortXDockerConfig,
    username: str,
    token: str,
    input_dir: str,
    output_dir: str,
    sig_filename: str,
    mix_filename: str,
) -> list[str]:
    """Construct declarative CLI command argument list for running CIBERSORTx container."""
    sig_basename = Path(sig_filename).name
    mix_basename = Path(mix_filename).name

    cmd = [
        config.docker_binary,
        "run",
        "--rm",
        "--platform",
        "linux/amd64",
        "-v",
        f"{input_dir}:/src/data",
        "-v",
        f"{output_dir}:/src/outdir",
        config.image_name,
        "--username",
        username,
        "--token",
        token,
        "--sigmatrix",
        sig_basename,
        "--mixture",
        mix_basename,
        "--perm",
        str(config.n_perm),
        "--rmbatchBmode",
        "TRUE" if config.rmbatch_bmode else "FALSE",
        "--rmbatchSmode",
        "TRUE" if config.rmbatch_smode else "FALSE",
        "--QN",
        "TRUE" if config.qn else "FALSE",
        "--verbose",
        "TRUE",
    ]

    if config.absolute:
        cmd.extend(["--absolute", "TRUE", "--abs_method", config.abs_method])

    return cmd


def parse_cibersortx_results(content: str) -> Result[CibersortXResult, str]:
    """
    Purely parse the tab-delimited output table from CIBERSORTx into a CibersortXResult.

    Extracts proportions, sample names, cell types, empirical p-values, correlation, and RMSE.
    """
    cleaned_content = content.strip()
    if not cleaned_content:
        return Failure("CIBERSORTx output is empty.")

    try:
        df = pl.read_csv(io.StringIO(cleaned_content), separator="\t")
    except Exception as e:
        return Failure(f"Failed to parse CIBERSORTx output table: {e}")

    cols = df.columns
    if len(cols) < 2:
        return Failure(
            f"Expected at least 2 columns in CIBERSORTx results, got {cols}"
        )

    sample_col = cols[0]
    sample_names = tuple(df[sample_col].cast(pl.String).to_list())

    # Identify diagnostic columns (case-insensitive substring matches)
    p_val_col = next(
        (c for c in cols if "p-value" in c.lower() or "p.value" in c.lower()),
        None,
    )
    corr_col = next((c for c in cols if "correlation" in c.lower()), None)
    rmse_col = next((c for c in cols if "rmse" in c.lower()), None)
    abs_col = next((c for c in cols if "absolute" in c.lower()), None)

    diagnostic_cols = {c for c in [sample_col, p_val_col, corr_col, rmse_col, abs_col] if c is not None}
    cell_type_cols = tuple(c for c in cols if c not in diagnostic_cols)

    if not cell_type_cols:
        return Failure(
            "No cell type fraction columns identified in CIBERSORTx output."
        )

    # Cast cell type columns to Float64
    cast_exprs = [pl.col(sample_col).cast(pl.String)] + [
        pl.col(c).cast(pl.Float64) for c in cell_type_cols
    ]
    proportions_df = df.select(cast_exprs)

    p_values: Maybe[pl.Series] = (
        Some(df[p_val_col].cast(pl.Float64).alias("p_value"))
        if p_val_col
        else Nothing
    )
    correlation: pl.Series = (
        df[corr_col].cast(pl.Float64).alias("correlation")
        if corr_col
        else pl.Series("correlation", [0.0] * df.height, dtype=pl.Float64)
    )
    rmse: pl.Series = (
        df[rmse_col].cast(pl.Float64).alias("rmse")
        if rmse_col
        else pl.Series("rmse", [0.0] * df.height, dtype=pl.Float64)
    )
    absolute_scores: Maybe[pl.Series] = (
        Some(df[abs_col].cast(pl.Float64).alias("absolute_score"))
        if abs_col
        else Nothing
    )

    return Success(
        CibersortXResult(
            proportions=proportions_df,
            cell_types=cell_type_cols,
            sample_names=sample_names,
            p_values=p_values,
            correlation=correlation,
            rmse=rmse,
            absolute_scores=absolute_scores,
            raw_output_text=Some(cleaned_content),
        )
    )


# =========================================================================
# Imperative Shell Functions
# =========================================================================


def is_docker_available(docker_binary: str = "docker") -> bool:
    """Check if the Docker binary and daemon are reachable."""
    try:
        proc = subprocess.run(
            [docker_binary, "info"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=5,
            check=False,
        )
        return proc.returncode == 0
    except (FileNotFoundError, subprocess.SubprocessError, PermissionError):
        return False


def run_cibersortx_docker(
    mixture: np.ndarray,
    signature: np.ndarray,
    gene_names: Sequence[str],
    cell_types: Sequence[str],
    sample_names: Sequence[str],
    config: CibersortXDockerConfig,
    mock_runner: Maybe[
        Callable[[list[str], Path, Path], Result[str, str]]
    ] = Nothing,
) -> Result[CibersortXResult, str]:
    """
    Run CIBERSORTx deconvolution using the official Docker container.

    Parameters:
        mixture: 2D array of shape (genes x samples)
        signature: 2D array of shape (genes x cell_types)
        gene_names: Gene symbols corresponding to rows
        cell_types: Cell type labels corresponding to signature columns
        sample_names: Sample labels corresponding to mixture columns
        config: CibersortXDockerConfig container settings
        mock_runner: Optional injected callable for offline testing: (cmd, in_dir, out_dir) -> Result[content, err]

    Returns:
        Result[CibersortXResult, str] containing inferred cell fractions and diagnostics.
    """
    import tempfile

    # 1. Resolve credentials
    creds_res = resolve_cibersortx_credentials(config)
    match creds_res:
        case Failure(err):
            return Failure(err)
        case Success((username, token)):
            pass

    # 2. Format input matrices
    sig_res = format_cibersort_matrix(signature, gene_names, cell_types)
    match sig_res:
        case Failure(err):
            return Failure(f"Signature matrix formatting error: {err}")
        case Success(sig_text):
            pass

    mix_res = format_cibersort_matrix(mixture, gene_names, sample_names)
    match mix_res:
        case Failure(err):
            return Failure(f"Mixture matrix formatting error: {err}")
        case Success(mix_text):
            pass

    # 3. Create temporary workspace directories
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmpdir_raw:
        tmpdir = Path(tmpdir_raw)
        in_dir = tmpdir / "input"
        out_dir = tmpdir / "output"
        in_dir.mkdir(parents=True, exist_ok=True)
        out_dir.mkdir(parents=True, exist_ok=True)

        sig_file = in_dir / "sigmatrix.txt"
        mix_file = in_dir / "mixture.txt"
        sig_file.write_text(sig_text)
        mix_file.write_text(mix_text)

        cmd = build_cibersortx_docker_command(
            config=config,
            username=username,
            token=token,
            input_dir=str(in_dir),
            output_dir=str(out_dir),
            sig_filename="sigmatrix.txt",
            mix_filename="mixture.txt",
        )

        # 4. Injected mock runner branch (for testing without Docker)
        match mock_runner:
            case Some(runner):
                run_res = runner(cmd, in_dir, out_dir)
                match run_res:
                    case Failure(err):
                        return Failure(err)
                    case Success(output_text):
                        return parse_cibersortx_results(output_text)
            case Nothing:
                pass

        # 5. Live Docker execution
        if not is_docker_available(config.docker_binary):
            return Failure(
                f"Docker binary '{config.docker_binary}' or Docker daemon is unavailable."
            )

        try:
            proc = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=config.timeout_seconds,
                check=False,
                text=True,
            )
        except subprocess.TimeoutExpired:
            return Failure(
                f"CIBERSORTx container timed out after {config.timeout_seconds} seconds."
            )
        except Exception as e:
            return Failure(f"Failed to invoke Docker container: {e}")

        if proc.returncode != 0:
            return Failure(
                f"CIBERSORTx container exited with code {proc.returncode}.\n"
                f"Stderr: {proc.stderr.strip() if proc.stderr else '(empty)'}\n"
                f"Stdout: {proc.stdout.strip() if proc.stdout else '(empty)'}"
            )

        # 6. Locate output file (prioritize *Results*.txt, then *Adjusted*.txt, then any non-empty .txt/.csv/.tsv)
        all_output_files = [
            p
            for p in out_dir.rglob("*")
            if p.is_file()
            and p.suffix.lower() in {".txt", ".tsv", ".csv"}
            and p.stat().st_size > 0
        ]
        results_candidates = [
            p for p in all_output_files if "results" in p.name.lower()
        ]
        if not results_candidates:
            results_candidates = all_output_files

        if not results_candidates:
            return Failure(
                f"No output results file found in CIBERSORTx output directory ({out_dir}).\n"
                f"Container Exit Code: {proc.returncode}\n"
                f"Stdout:\n{proc.stdout.strip() if proc.stdout else '(empty)'}\n"
                f"Stderr:\n{proc.stderr.strip() if proc.stderr else '(empty)'}"
            )

        # Pick the most recently modified results file
        results_file = sorted(results_candidates, key=lambda p: p.stat().st_mtime)[-1]
        output_content = results_file.read_text(encoding="utf-8", errors="replace")

        return parse_cibersortx_results(output_content)
