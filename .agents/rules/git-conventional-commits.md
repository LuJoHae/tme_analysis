---
trigger: always_on
---

# Conventional Commits

This project strictly enforces the **[Conventional Commits v1.0.0](https://www.conventionalcommits.org/en/v1.0.0/)** specification for all Git commit messages. When creating commits or suggesting commit messages, you must follow this format.

## Format
```
<type>[optional scope]: <description>

[optional body]

[optional footer(s)]
```

## Allowed Types
- **feat**: A new feature
- **fix**: A bug fix
- **docs**: Documentation only changes
- **style**: Changes that do not affect the meaning of the code (white-space, formatting, etc.)
- **refactor**: A code change that neither fixes a bug nor adds a feature
- **perf**: A code change that improves performance
- **test**: Adding missing tests or correcting existing tests
- **build**: Changes that affect the build system or external dependencies (e.g., pip, uv, npm)
- **ci**: Changes to CI configuration files and scripts
- **chore**: Other changes that don't modify `src` or `test` files
- **revert**: Reverts a previous commit

## Rules
1. **Lowercase Types**: The `<type>` must be in all lowercase.
2. **Imperative Mood**: The `<description>` must be written in the imperative, present tense (e.g., "add feature" not "added feature").
3. **No Period**: Do not capitalize the first letter of the description, and do not end it with a period.
4. **Breaking Changes**: If the commit introduces a breaking change, append a `!` after the type/scope (e.g., `feat!: drop python 3.8 support`), or indicate `BREAKING CHANGE: <description>` in the footer.
5. **Small, Logically Sorted Commits**: Keep commits small, focused, and easy to review (target roughly **50–250 changed lines** and **one logical concern** per commit). Never bundle large multi-layer changes (e.g. 500+ lines across packages, scripts, and workflows) into a single commit. When a task touches many lines, split it into multiple logically ordered commits along architectural boundaries:
   1. `chore(deps)` / `build` (`pyproject.toml`, `uv.lock`)
   2. `feat` / `fix` / `refactor` in core packages (`packages/tme_datasets`, `packages/tme_response`, `packages/plotting_utils`)
   3. `test` (unit and property-based tests)
   4. `feat` / `refactor` in analysis/CLI scripts (`scripts/`)
   5. `chore(workflow)` (`workflow/`, `Makefile`)
   6. `docs` (`README`, `article/`, `.agents/rules/`)
6. **Regular & Frequent Cadence**: Commit frequently and incrementally as work progresses. Never postpone committing until the end of an entire multi-step task or turn. Proactively commit after completing each meaningful milestone (e.g., after creating a file, writing a unit test, updating metadata, or completing a specific refactor).
7. **Selective Staging**: Always selectively stage only the relevant files (`git add <file>`) rather than blanket adding (`git add .` or `git add -A`) to avoid accidentally committing temporary files, scratch scripts, or intermediate build artifacts.
