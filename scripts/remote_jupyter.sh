#!/usr/bin/env bash
set -euo pipefail

SESSION_NAME="jupyter_tme"
PORT="${JUPYTER_PORT:-8888}"
UV_BIN="${REMOTE_UV:-/home/halu/.local/bin/uv}"
WORK_DIR="${REMOTE_DIR:-$HOME/python-venv/tme_analysis}"
LOG_FILE="$HOME/.jupyter_tme.log"

case "${1:-status}" in
    start)
        if tmux has-session -t "$SESSION_NAME" 2>/dev/null; then
            echo "Jupyter server is already running in tmux session '$SESSION_NAME'."
        else
            echo "Starting JupyterLab in tmux session '$SESSION_NAME' on port $PORT..."
            cd "$WORK_DIR"
            tmux new-session -d -s "$SESSION_NAME" \
                "cd '$WORK_DIR' && '$UV_BIN' run jupyter lab --no-browser --port=$PORT --ip=127.0.0.1 --notebook-dir='$WORK_DIR' 2>&1 | tee '$LOG_FILE'"
            sleep 3
        fi
        "$0" status
        ;;

    status)
        if tmux has-session -t "$SESSION_NAME" 2>/dev/null; then
            echo "=== Session Status ==="
            echo "tmux session '$SESSION_NAME': ACTIVE"
            echo ""
            echo "=== Active Jupyter Servers ==="
            cd "$WORK_DIR"
            "$UV_BIN" run jupyter server list 2>/dev/null || true
            echo ""
            echo "=== Connection Instructions ==="
            echo "Run this command on your LOCAL machine to forward the port:"
            echo "  ssh -N -L $PORT:127.0.0.1:$PORT olm"
            echo ""
            echo "Then open your browser to the URL (including token) shown above."
        else
            echo "tmux session '$SESSION_NAME' is NOT running."
            cd "$WORK_DIR"
            "$UV_BIN" run jupyter server list 2>/dev/null || true
        fi
        ;;

    stop)
        if tmux has-session -t "$SESSION_NAME" 2>/dev/null; then
            echo "Stopping tmux session '$SESSION_NAME'..."
            tmux kill-session -t "$SESSION_NAME" || true
            cd "$WORK_DIR"
            "$UV_BIN" run jupyter server stop "$PORT" 2>/dev/null || true
            echo "Jupyter server stopped."
        else
            echo "tmux session '$SESSION_NAME' is not running."
            cd "$WORK_DIR"
            "$UV_BIN" run jupyter server stop "$PORT" 2>/dev/null || true
        fi
        ;;

    logs)
        if [ -f "$LOG_FILE" ]; then
            tail -n 50 "$LOG_FILE"
        elif tmux has-session -t "$SESSION_NAME" 2>/dev/null; then
            tmux capture-pane -pt "$SESSION_NAME" -S -50
        else
            echo "No log file or running session found."
        fi
        ;;

    *)
        echo "Usage: $0 {start|status|stop|logs}"
        exit 1
        ;;
esac
