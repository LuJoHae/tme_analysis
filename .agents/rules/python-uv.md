# Execute Python with `uv`

This project uses `uv` for Python dependency management and script execution.

## Rule
Whenever you need to run a Python script, module, or tool (such as `pytest`, `snakemake`, etc.) locally or remotely, you must **always** execute it using `uv run`. 

Do not use the bare `python` command, and do not activate the virtual environment manually.

**Correct:**
```bash
uv run python script.py
uv run pytest
```

**Incorrect:**
```bash
python script.py
pytest
source .venv/bin/activate && python script.py
```
