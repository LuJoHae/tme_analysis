#!/usr/bin/env bash
set -euo pipefail

JOB_NAME="jupyter_tme"
PORT="${JUPYTER_PORT:-8888}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-em25}"
SLURM_PARTITION="${SLURM_PARTITION:-normal}"
SCRATCH_BASE="/ritom/scratch/cscs/${USER}"
UV_BIN="${REMOTE_UV:-/users/${USER}/.local/bin/uv}"
WORK_DIR="${REMOTE_DIR:-$SCRATCH_BASE/code}"
DATA_DIR="${TME_DATA_DIR:-$SCRATCH_BASE/data-test}"
TMP_DIR="${TMPDIR:-$SCRATCH_BASE/tmp}"
LOG_DIR="$SCRATCH_BASE/slurm_logs"
LOG_FILE="$LOG_DIR/jupyter_tme.log"

mkdir -p "$LOG_DIR" "$TMP_DIR"

case "${1:-status}" in
    start)
        ACTIVE_JOB=$(squeue -u "$USER" -n "$JOB_NAME" -h -o "%i" 2>/dev/null || true)
        if [ -n "$ACTIVE_JOB" ]; then
            echo "Jupyter server is already running in Slurm job $ACTIVE_JOB."
        else
            echo "Submitting JupyterLab Slurm job '$JOB_NAME' on port $PORT..."
            sbatch -A "$SLURM_ACCOUNT" -p "$SLURM_PARTITION" -J "$JOB_NAME" -n 1 -c 8 --mem=32G -t 04:00:00 \
                -o "$LOG_FILE" \
                --wrap="cd '$WORK_DIR' && TME_DATA_DIR='$DATA_DIR' TME_DISABLE_DOWNLOADS=1 TMPDIR='$TMP_DIR' UV_CACHE_DIR='$SCRATCH_BASE/.cache/uv' '$UV_BIN' run jupyter lab --no-browser --port=$PORT --ip=0.0.0.0 --notebook-dir='$WORK_DIR'"
            sleep 3
        fi
        "$0" status
        ;;

    status)
        ACTIVE_INFO=$(squeue -u "$USER" -n "$JOB_NAME" -h -o "%i %T %N" 2>/dev/null || true)
        if [ -n "$ACTIVE_INFO" ]; then
            echo "=== Slurm Job Status ==="
            echo "Job '$JOB_NAME': $ACTIVE_INFO"
            NODE=$(echo "$ACTIVE_INFO" | awk '{print $3}')
            echo ""
            echo "=== Connection Instructions ==="
            echo "Run this command on your LOCAL machine to forward the port from compute node $NODE:"
            echo "  ssh -N -L $PORT:$NODE:$PORT eiger"
            echo ""
            if [ -f "$LOG_FILE" ]; then
                echo "=== Recent Jupyter Log (for token URL) ==="
                tail -n 20 "$LOG_FILE"
            fi
        else
            echo "Slurm job '$JOB_NAME' is NOT running."
        fi
        ;;

    stop)
        echo "Cancelling Slurm job '$JOB_NAME'..."
        scancel -u "$USER" -n "$JOB_NAME" 2>/dev/null || true
        echo "Jupyter Slurm job stopped."
        ;;

    logs)
        if [ -f "$LOG_FILE" ]; then
            tail -n 50 "$LOG_FILE"
        else
            echo "No log file found at $LOG_FILE."
        fi
        ;;

    *)
        echo "Usage: $0 {start|status|stop|logs}"
        exit 1
        ;;
esac
