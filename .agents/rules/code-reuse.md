# Code Reuse, Direct Upstreaming, and Anti-Duplication Policy

This repository is designed with modular packages in `packages/*` and thin analysis scripts in `scripts/` and `workflow/`. Writing redundant code, duplicate helper functions, or ad-hoc data processors is strictly prohibited.

---

## 1. The Direct Upstreaming Protocol

All domain logic must live in `packages/*`, not in scripts:
- **Scripts are Pure Orchestration**: Files in `scripts/` and `workflow/` are entry points. They should only handle command-line arguments, orchestrate high-level pipeline calls, and read/write files.
- **No Script-Local Domain Helpers**: Do not implement data parsing, biological conversions, feature engineering, statistical modeling, or custom transformations directly inside script files.
- **Dataset Loaders and Builders**: All dataset ingestion, custom downloaders (e.g., Google Drive, GEO), format converters, and AnnData matrix builders MUST be implemented as reusable functions within `packages/tme_datasets/src/tme_datasets/` (in `providers/` and `download/`), accompanied by unit tests in `packages/tme_datasets/tests/`. Do NOT write one-off standalone scripts in the root directory or `scratch/` for dataset processing.
- **Direct Upstreaming**: If a script needs new functionality or a modification to existing behavior:
  1. Add or generalize the function directly inside the relevant package under `packages/*` (see [.agents/rules/package-map.md](file:///Users/halu/Code/tme_analysis/.agents/rules/package-map.md)).
  2. Add or update corresponding unit tests in `packages/<pkg>/tests/` and verify they pass with `pytest`.
  3. Import the function into the script.

---

## 2. Mandatory Pre-Implementation Discovery Phase

Before creating any function, class, or data model, you MUST search the existing codebase:
1. **Search Workspace Packages**:
   - Use `grep_search` and `find_by_name` across `packages/` to inspect existing implementations, utilities, and abstractions.
   - Refer to [.agents/rules/package-map.md](file:///Users/halu/Code/tme_analysis/.agents/rules/package-map.md) to locate the package owning the target domain.
2. **Audit Imports**:
   - Inspect existing scripts and package modules doing similar tasks to identify existing entry points and conventions.
3. **Verify Before Coding**:
   - Creating a duplicate function that already exists in any package is considered a critical error.

---

## 3. Amend and Generalize Over Re-implementing

When an existing package function does 70–80% of what you need:
- **Generalize the Existing Function**: Extend the function with backward-compatible optional arguments, polymorphic dispatch via Typeclasses/Protocols, or extract a shared core subroutine.
- **Maintain Compatibility**: Ensure all existing callers and tests remain functional.
- **Do Not Fork Logic**: Never create parallel variants like `process_data_v2()`, `custom_convert()`, or `my_loader()`.

---

## 4. Mandatory Automated Testing for Package Upstreams

Any code added or amended in `packages/*` must adhere to repo quality standards:
- **Unit Tests Required**: Every new or amended function must be covered by unit tests in `packages/<pkg>/tests/`.
- **Run Tests**: Execute `uv run pytest packages/<pkg>` to ensure no regressions.
- **Functional Style**: All package code must adhere to [.agents/rules/code-style-guide.md](file:///Users/halu/Code/tme_analysis/.agents/rules/code-style-guide.md) (pure functions, frozen Pydantic/dataclasses, Polars, Result monad, `mypy --strict`).

---

## 5. Out-of-Domain Logic Guardrail

If a required capability does not cleanly fit into any of the 12 packages listed in [.agents/rules/package-map.md](file:///Users/halu/Code/tme_analysis/.agents/rules/package-map.md):
- **STOP and Consult the User**: Do NOT autonomously create a brand-new package or silently leave the helper functions trapped inside an analysis script.
- Present the requirement to the user and ask where the logic should reside (e.g. creating a new package or extending an existing one).

---

## 6. Transparency in Plans and Summaries

When producing an implementation plan (`implementation_plan.md`) or summarizing completed work:
- Explicitly list:
  - **Reused Package Functions**: Which existing functions/classes were leveraged.
  - **Amended Package Functions**: Which existing functions were modified/generalized and how backward compatibility was preserved.
  - **Newly Upstreamed Functions**: New functions added to `packages/*`, their target package, and test coverage.

---

## 7. Strict Dependency Direction (Zero-Dependency Foundation Invariant)

- **Foundation Isolation**: `packages/tme_datasets` is the foundational data layer (Layer 0) and must remain **100% free of internal repository dependencies**. Never add an import of a sibling package (`gene_utils`, `ml_pipelines`, etc.) into `tme_datasets`.
- **Strict Downward Dependency Flow**:
  - Scripts/Pipelines (Layer 2) → Specialized Packages (Layer 1) → Foundational Data Layer (`tme_datasets`, Layer 0).
  - Cross-imports between sibling packages in Layer 1 must not introduce circular dependencies.
- **Manifest Dependency Verification**: Whenever adding or modifying code in any package under `packages/*`, ensure all imported third-party packages are explicitly listed in that package's local `pyproject.toml`. Do not rely on ambient dependencies provided by other workspace members.
