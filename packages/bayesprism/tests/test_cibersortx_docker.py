"""Unit tests for CIBERSORTx Docker container adapter."""

from pathlib import Path
import numpy as np
import polars as pl
import pytest
from returns.maybe import Some, Nothing
from returns.result import Success, Failure

from bayesprism.adapters.cibersortx_docker import (
    CibersortXDockerConfig,
    CibersortXResult,
    format_cibersort_matrix,
    build_cibersortx_docker_command,
    resolve_cibersortx_credentials,
    parse_cibersortx_results,
    run_cibersortx_docker,
    is_docker_available,
)


def test_cibersortx_config_frozen() -> None:
    """Verify CibersortXDockerConfig is immutable and has expected defaults."""
    cfg = CibersortXDockerConfig(
        username="user@example.com",
        token="tok_12345",
        rmbatch_bmode=True,
    )
    assert cfg.username == "user@example.com"
    assert cfg.token == "tok_12345"
    assert cfg.rmbatch_bmode is True
    assert cfg.rmbatch_smode is False
    assert cfg.image_name == "cibersortx/fractions"

    with pytest.raises(Exception):
        cfg.rmbatch_bmode = False  # type: ignore


def test_format_cibersort_matrix_success() -> None:
    """Verify proper tab-separated formatting with GeneSymbol header."""
    genes = ("CD3D", "MS4A1", "CD8A")
    cell_types = ("T_cells", "B_cells")
    mat = np.array([
        [10.5, 0.2],
        [0.1, 25.8],
        [8.4, 0.0],
    ])

    match format_cibersort_matrix(mat, genes, cell_types):
        case Failure(err):
            raise AssertionError(f"Formatting failed: {err}")
        case Success(text):
            lines = text.strip().split("\n")
            assert len(lines) == 4
            assert lines[0] == "GeneSymbol\tT_cells\tB_cells"
            assert lines[1].startswith("CD3D\t10.5\t0.2")
            assert lines[2].startswith("MS4A1\t0.1\t25.8")
            assert lines[3].startswith("CD8A\t8.4\t0")


def test_format_cibersort_matrix_validation() -> None:
    """Verify validation on dimension mismatches and duplicate genes."""
    mat = np.ones((3, 2))

    # Row mismatch
    match format_cibersort_matrix(mat, ("G1", "G2"), ("C1", "C2")):
        case Failure(err):
            assert "Row count mismatch" in err
        case Success(_):
            raise AssertionError("Should fail on row count mismatch")

    # Column mismatch
    match format_cibersort_matrix(mat, ("G1", "G2", "G3"), ("C1",)):
        case Failure(err):
            assert "Column count mismatch" in err
        case Success(_):
            raise AssertionError("Should fail on column count mismatch")

    # Duplicate gene names
    match format_cibersort_matrix(mat, ("G1", "G1", "G3"), ("C1", "C2")):
        case Failure(err):
            assert "unique" in err.lower()
        case Success(_):
            raise AssertionError("Should fail on duplicate genes")


def test_resolve_cibersortx_credentials(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify credential resolution from config and environment variables."""
    # Direct config
    cfg1 = CibersortXDockerConfig(username="u1@test.org", token="tok1")
    match resolve_cibersortx_credentials(cfg1):
        case Success((u, t)):
            assert u == "u1@test.org"
            assert t == "tok1"
        case Failure(err):
            raise AssertionError(err)

    # Empty config, missing env vars, missing .env -> Failure
    cfg_empty = CibersortXDockerConfig()
    monkeypatch.delenv("CIBERSORTX_USERNAME", raising=False)
    monkeypatch.delenv("CIBERSORTX_TOKEN", raising=False)
    dummy_empty_env = Path("/tmp/non_existent_env_file_xyz")
    match resolve_cibersortx_credentials(cfg_empty, env_path=dummy_empty_env):
        case Failure(err):
            assert "username is required" in err.lower()
        case Success(_):
            raise AssertionError("Should fail when credentials absent")

    # Fallback to env vars
    monkeypatch.setenv("CIBERSORTX_USERNAME", "env_user@test.org")
    monkeypatch.setenv("CIBERSORTX_TOKEN", "env_tok_999")
    match resolve_cibersortx_credentials(cfg_empty, env_path=dummy_empty_env):
        case Success((u, t)):
            assert u == "env_user@test.org"
            assert t == "env_tok_999"
        case Failure(err):
            raise AssertionError(err)


def test_build_cibersortx_docker_command() -> None:
    """Verify synthesis of Docker CLI arguments."""
    cfg = CibersortXDockerConfig(
        rmbatch_bmode=True,
        rmbatch_smode=False,
        qn=False,
        n_perm=100,
        absolute=True,
    )
    cmd = build_cibersortx_docker_command(
        config=cfg,
        username="u@test.org",
        token="tok123",
        input_dir="/tmp/input",
        output_dir="/tmp/output",
        sig_filename="sig.txt",
        mix_filename="mix.txt",
    )

    assert cmd[0] == "docker"
    assert cmd[1] == "run"
    assert "-v" in cmd
    assert "/tmp/input:/src/data" in cmd
    assert "/tmp/output:/src/outdir" in cmd
    assert "--username" in cmd and cmd[cmd.index("--username") + 1] == "u@test.org"
    assert "--sigmatrix" in cmd and cmd[cmd.index("--sigmatrix") + 1] == "sig.txt"
    assert "--mixture" in cmd and cmd[cmd.index("--mixture") + 1] == "mix.txt"
    assert "--perm" in cmd and cmd[cmd.index("--perm") + 1] == "100"
    assert "--rmbatchBmode" in cmd and cmd[cmd.index("--rmbatchBmode") + 1] == "TRUE"
    assert "--rmbatchSmode" in cmd and cmd[cmd.index("--rmbatchSmode") + 1] == "FALSE"
    assert "--QN" in cmd and cmd[cmd.index("--QN") + 1] == "FALSE"
    assert "--absolute" in cmd and cmd[cmd.index("--absolute") + 1] == "TRUE"


def test_parse_cibersortx_results() -> None:
    """Verify parsing of CIBERSORTx output TSV with cell types and diagnostics."""
    mock_tsv = (
        "Mixture\tB.cells\tT.cells\tNK.cells\tP-value\tCorrelation\tRMSE\tAbsolute score (sig.score)\n"
        "Sample_01\t0.20\t0.50\t0.30\t0.01\t0.88\t0.32\t1.45\n"
        "Sample_02\t0.10\t0.70\t0.20\t0.00\t0.94\t0.21\t1.80\n"
    )

    match parse_cibersortx_results(mock_tsv):
        case Failure(err):
            raise AssertionError(f"Parsing failed: {err}")
        case Success(res):
            assert isinstance(res, CibersortXResult)
            assert res.sample_names == ("Sample_01", "Sample_02")
            assert res.cell_types == ("B.cells", "T.cells", "NK.cells")

            # Check proportions DataFrame
            assert "Mixture" in res.proportions.columns
            assert "B.cells" in res.proportions.columns
            b_vals = res.proportions["B.cells"].to_list()
            assert b_vals == [0.20, 0.10]

            # Check diagnostics
            assert res.p_values != Nothing
            assert res.p_values.unwrap().to_list() == [0.01, 0.00]
            assert res.correlation.to_list() == [0.88, 0.94]
            assert res.rmse.to_list() == [0.32, 0.21]
            assert res.absolute_scores != Nothing
            assert res.absolute_scores.unwrap().to_list() == [1.45, 1.80]


def test_run_cibersortx_docker_mock() -> None:
    """Test full workflow using an injected mock runner simulating Docker execution."""
    rng = np.random.default_rng(42)
    G = 10
    K = 3
    N = 2
    genes = tuple(f"Gene_{i:02d}" for i in range(G))
    cell_types = ("T_cells", "B_cells", "Myeloid")
    sample_names = ("Patient_A", "Patient_B")

    sig = rng.gamma(2.0, 1.0, (G, K))
    mix = rng.gamma(2.0, 1.0, (G, N))

    cfg = CibersortXDockerConfig(
        username="test@lab.org",
        token="tok_secret",
        n_perm=10,
    )

    def mock_container_runner(
        cmd: list[str],
        in_dir: Path,
        out_dir: Path,
    ) -> Result[str, str]:
        # Assert input files were generated correctly
        sig_file = in_dir / "sigmatrix.txt"
        mix_file = in_dir / "mixture.txt"
        assert sig_file.exists()
        assert mix_file.exists()
        assert "GeneSymbol" in sig_file.read_text()
        assert "GeneSymbol" in mix_file.read_text()

        # Simulate CIBERSORTx writing results
        simulated_output = (
            "Mixture\tT_cells\tB_cells\tMyeloid\tP-value\tCorrelation\tRMSE\n"
            "Patient_A\t0.50\t0.30\t0.20\t0.02\t0.91\t0.25\n"
            "Patient_B\t0.20\t0.70\t0.10\t0.00\t0.95\t0.18\n"
        )
        return Success(simulated_output)

    match run_cibersortx_docker(
        mixture=mix,
        signature=sig,
        gene_names=genes,
        cell_types=cell_types,
        sample_names=sample_names,
        config=cfg,
        mock_runner=Some(mock_container_runner),
    ):
        case Failure(err):
            raise AssertionError(f"run_cibersortx_docker failed: {err}")
        case Success(res):
            assert res.sample_names == sample_names
            assert res.cell_types == cell_types
            assert res.proportions.shape == (2, 4)  # Mixture + 3 cell types
            assert res.correlation.to_list() == [0.91, 0.95]


def test_docker_missing_graceful_failure() -> None:
    """Verify that an unavailable Docker binary produces an explicit Failure."""
    cfg = CibersortXDockerConfig(
        username="u@test.org",
        token="tok",
        docker_binary="non_existent_docker_binary_xyz",
    )
    sig = np.ones((5, 2))
    mix = np.ones((5, 2))
    genes = ("G1", "G2", "G3", "G4", "G5")
    cell_types = ("C1", "C2")
    sample_names = ("S1", "S2")

    match run_cibersortx_docker(
        mixture=mix,
        signature=sig,
        gene_names=genes,
        cell_types=cell_types,
        sample_names=sample_names,
        config=cfg,
    ):
        case Failure(err):
            assert "unavailable" in err.lower() or "docker" in err.lower()
        case Success(_):
            raise AssertionError("Should have failed gracefully when Docker is missing")
