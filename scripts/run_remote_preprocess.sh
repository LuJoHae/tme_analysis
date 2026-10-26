#!/usr/bin/env bash
# ==============================================================================
# Remote Slurm Control Wrapper for TME Dataset Preprocessing Pipeline
# Runs strictly via Slurm on CSCS Alps Eiger (internet downloads disabled).
#
# Commands:
#   bash scripts/run_remote_preprocess.sh start   # Submit Slurm batch job
#   bash scripts/run_remote_preprocess.sh status  # Inspect Slurm job state & progress
#   bash scripts/run_remote_preprocess.sh logs    # Display recent log output (-f to follow)
#   bash scripts/run_remote_preprocess.sh attach  # Show Slurm queue status
#   bash scripts/run_remote_preprocess.sh stop    # Cancel Slurm job
# ==============================================================================

set -euo pipefail

JOB_NAME="tme_preprocess"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-em25}"
SLURM_PARTITION="${SLURM_PARTITION:-normal}"
WORKERS="${WORKERS:-8}"
MAX_RAM_GB="${MAX_RAM_GB:-140}"
SCRATCH_BASE="/ritom/scratch/cscs/${USER}"
UV_BIN="${REMOTE_UV:-/users/${USER}/.local/bin/uv}"
if [ ! -x "$UV_BIN" ]; then
    UV_BIN="$(which uv 2>/dev/null || echo "uv")"
fi

WORK_DIR="${REMOTE_DIR:-$SCRATCH_BASE/code}"
DATA_DIR="${TME_DATA_DIR:-$SCRATCH_BASE/data-test}"
TMP_DIR="${TMPDIR:-$SCRATCH_BASE/tmp}"
LOG_DIR="$SCRATCH_BASE/slurm_logs"
OUTPUT_DIR="$WORK_DIR/output"
LOG_FILE="$LOG_DIR/tme_preprocess.log"
PROGRESS_FILE="$OUTPUT_DIR/preprocessing_progress.json"

mkdir -p "$LOG_DIR" "$OUTPUT_DIR" "$TMP_DIR"

case "${1:-status}" in
    start)
        echo "=========================================================================================="
        echo "ERROR: Internet dataset downloads are strictly disabled on CSCS Alps Eiger."
        echo "Download raw datasets on 'olm' first, transfer to /capstor/store/cscs/empa/em25/halu/data-test,"
        echo "and run 'make stage-in' + 'make preprocess-discovered-bg' to preprocess via Slurm."
        echo "=========================================================================================="
        exit 1
        ;;

    status)
        echo "=================================================================="
        echo "TME Slurm Preprocessing Status on $(hostname)"
        echo "=================================================================="
        if command -v squeue >/dev/null 2>&1; then
            squeue -u "$USER" -n "$JOB_NAME,preprocess_discovered" || true
        fi
        if [ -f "$PROGRESS_FILE" ]; then
            echo ""
            echo "--- Latest Checkpoint Progress ---"
            python3 -c "
import json
try:
    with open('$PROGRESS_FILE') as f:
        data = json.load(f)
    print(f'Last Updated:        {data.get(\"updated_at\", \"N/A\")}')
    print(f'Completed:           {data.get(\"completed_count\", 0)} / {data.get(\"total_target_datasets\", \"N/A\")}')
    print(f'Failed:              {data.get(\"failed_count\", 0)}')
except Exception as e:
    print(f'Could not parse checkpoint: {e}')
" 2>/dev/null || true
        fi
        if [ -f "$LOG_FILE" ]; then
            echo ""
            echo "--- Recent Log Output ---"
            tail -n 12 "$LOG_FILE"
        fi
        echo "=================================================================="
        ;;

    logs)
        if [ ! -f "$LOG_FILE" ]; then
            echo "No log file found at $LOG_FILE."
            exit 0
        fi
        if [ "${2:-}" = "-f" ] || [ "${2:-}" = "--follow" ]; then
            tail -f -n 50 "$LOG_FILE"
        else
            tail -n 50 "$LOG_FILE"
        fi
        ;;

    attach)
        squeue -u "$USER" -n "$JOB_NAME,preprocess_discovered" || true
        ;;

    stop)
        echo "Cancelling Slurm preprocessing job(s)..."
        scancel -u "$USER" -n "$JOB_NAME" 2>/dev/null || true
        scancel -u "$USER" -n "preprocess_discovered" 2>/dev/null || true
        echo "Slurm preprocessing job(s) cancelled."
        ;;

    *)
        echo "Usage: $0 {start|status|logs|attach|stop}"
        exit 1
        ;;
esac
