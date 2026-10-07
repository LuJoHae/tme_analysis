#!/usr/bin/env bash
# ==============================================================================
# Remote Daemon Control Wrapper for TME Dataset Preprocessing Pipeline
# Survives remote SSH connection drop via detached tmux session or nohup.
#
# Commands:
#   bash scripts/run_remote_preprocess.sh start   # Launch detached background daemon
#   bash scripts/run_remote_preprocess.sh status  # Inspect daemon state, RAM & progress
#   bash scripts/run_remote_preprocess.sh logs    # Display recent log output (-f to follow)
#   bash scripts/run_remote_preprocess.sh attach  # Attach directly to running tmux session
#   bash scripts/run_remote_preprocess.sh stop    # Gracefully terminate daemon and workers
# ==============================================================================

set -euo pipefail

SESSION_NAME="tme_preprocess"
WORKERS="${WORKERS:-8}"
MAX_RAM_GB="${MAX_RAM_GB:-140.0}"
UV_BIN="${REMOTE_UV:-/home/halu/.local/bin/uv}"
if [ ! -x "$UV_BIN" ]; then
    UV_BIN="$(which uv 2>/dev/null || echo "uv")"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$WORK_DIR/logs"
OUTPUT_DIR="$WORK_DIR/output"
LOG_FILE="$LOG_DIR/download_and_preprocess.log"
PID_FILE="$LOG_DIR/preprocess.pid"
PROGRESS_FILE="$OUTPUT_DIR/preprocessing_progress.json"

mkdir -p "$LOG_DIR" "$OUTPUT_DIR"

is_tmux_running() {
    tmux has-session -t "$SESSION_NAME" 2>/dev/null
}

is_pid_running() {
    if [ -f "$PID_FILE" ]; then
        local pid
        pid=$(cat "$PID_FILE" 2>/dev/null || true)
        if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then
            return 0
        fi
    fi
    return 1
}

case "${1:-status}" in
    start)
        echo "=================================================================="
        echo "Starting Remote Preprocessing Daemon on $(hostname)"
        echo "=================================================================="

        if is_tmux_running; then
            echo "[WARNING] Preprocessing daemon is ALREADY running in tmux session '$SESSION_NAME'."
            echo "Use '$0 status' or '$0 attach' to monitor execution."
            exit 0
        fi

        if is_pid_running; then
            echo "[WARNING] Process is already active with PID $(cat "$PID_FILE")."
            exit 0
        fi

        cd "$WORK_DIR"

        if command -v tmux >/dev/null 2>&1; then
            echo "--> Spawning detached tmux session '$SESSION_NAME'..."
            tmux new-session -d -s "$SESSION_NAME" \
                "cd '$WORK_DIR' && TMPDIR=/storage/halu/tmp '$UV_BIN' run python scripts/download_and_preprocess_all.py --workers $WORKERS --max-ram-gb $MAX_RAM_GB 2>&1 | tee -a '$LOG_FILE'"
            
            sleep 2
            # Save pane PID to pidfile
            tmux list-panes -t "$SESSION_NAME" -F "#{pane_pid}" > "$PID_FILE" 2>/dev/null || true
            echo "--> Daemon successfully launched in tmux session '$SESSION_NAME'."
        else
            echo "--> tmux not found. Spawning via nohup with SIGHUP detachment..."
            nohup "$UV_BIN" run python scripts/download_and_preprocess_all.py --workers "$WORKERS" --max-ram-gb "$MAX_RAM_GB" >> "$LOG_FILE" 2>&1 < /dev/null &
            echo $! > "$PID_FILE"
            echo "--> Daemon successfully launched in background (PID: $(cat "$PID_FILE"))."
        fi

        echo ""
        "$0" status
        ;;

    status)
        echo "=================================================================="
        echo "TME Preprocessing Daemon Status: $(hostname)"
        echo "=================================================================="

        ACTIVE=0
        if is_tmux_running; then
            echo "Session Status:      ACTIVE (tmux session: '$SESSION_NAME')"
            ACTIVE=1
        elif is_pid_running; then
            echo "Session Status:      ACTIVE (PID: $(cat "$PID_FILE"))"
            ACTIVE=1
        else
            echo "Session Status:      STOPPED / INACTIVE"
        fi

        # Process tree metrics
        if [ "$ACTIVE" -eq 1 ]; then
            echo ""
            echo "--- Active Processes ---"
            ps -u "$USER" -o pid,%cpu,%mem,rss,etime,command | grep -E "download_and_preprocess|tme_preprocess" | grep -v grep || true
        fi

        # Checkpoint JSON summary
        if [ -f "$PROGRESS_FILE" ]; then
            echo ""
            echo "--- Latest Checkpoint Progress ---"
            python3 -c "
import json
try:
    with open('$PROGRESS_FILE') as f:
        data = json.load(f)
    print(f'Last Updated:        {data.get(\"updated_at\", \"N/A\")}')
    print(f'Target Datasets:     {data.get(\"total_target_datasets\", \"N/A\")}')
    print(f'Completed:           {data.get(\"completed_count\", 0)} / {data.get(\"total_target_datasets\", \"N/A\")}')
    print(f'Failed:              {data.get(\"failed_count\", 0)}')
    print(f'Remaining:           {data.get(\"pending_count\", \"N/A\")}')
    print(f'Active Workers:      {data.get(\"active_workers\", 0)}')
    print(f'Process Tree RAM:    {data.get(\"current_rss_gb\", 0.0):.1f} GB (Safety Cap: {data.get(\"max_ram_cap_gb\", 140.0):.1f} GB)')
    failed = data.get('failed_datasets', [])
    if failed:
        print('Failed Datasets:')
        for item in failed[:5]:
            print(f'  - {item.get(\"id\")}: {item.get(\"error\", \"\")[:80]}...')
except Exception as e:
    print(f'Could not parse checkpoint: {e}')
" 2>/dev/null || true
        fi

        # Log preview
        if [ -f "$LOG_FILE" ]; then
            echo ""
            echo "--- Recent Log Output (Last 12 Lines) ---"
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
        if is_tmux_running; then
            echo "Attaching to tmux session '$SESSION_NAME'. (Detach using Ctrl+B then D)..."
            tmux attach-session -t "$SESSION_NAME"
        else
            echo "tmux session '$SESSION_NAME' is not running."
            exit 1
        fi
        ;;

    stop)
        echo "Stopping preprocessing daemon..."
        if is_tmux_running; then
            echo "--> Killing tmux session '$SESSION_NAME'..."
            tmux kill-session -t "$SESSION_NAME" 2>/dev/null || true
        fi

        if [ -f "$PID_FILE" ]; then
            PID=$(cat "$PID_FILE" 2>/dev/null || true)
            if [ -n "$PID" ] && kill -0 "$PID" 2>/dev/null; then
                echo "--> Sending SIGTERM to PID $PID..."
                kill -TERM "$PID" 2>/dev/null || true
                sleep 2
                if kill -0 "$PID" 2>/dev/null; then
                    echo "--> Sending SIGKILL to PID $PID..."
                    kill -9 "$PID" 2>/dev/null || true
                fi
            fi
            rm -f "$PID_FILE"
        fi

        # Clean any remaining worker processes matching the script name
        pkill -f "download_and_preprocess_all.py" 2>/dev/null || true
        echo "Preprocessing daemon stopped."
        ;;

    *)
        echo "Usage: $0 {start|status|logs|attach|stop}"
        exit 1
        ;;
esac
