# Conventional Commits

This project strictly enforces the **Conventional Commits** specification for all Git commit messages. When creating commits or suggesting commit messages, you must follow this format.

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
5. **Small & Atomic Commits**: Keep commits small and focused. Do not bundle multiple or unrelated changes into massive commits. Break your work down into small, logically separated commits that address a single coherent concern (e.g., separate a `refactor` from a `feat`).
6. **Regular & Frequent Cadence**: Commit frequently and incrementally as work progresses. Never postpone committing until the end of an entire multi-step task or turn. Proactively commit after completing each meaningful milestone (e.g., after creating a file, writing a unit test, updating metadata, or completing a specific refactor).
7. **Selective Staging**: Always selectively stage only the relevant files (`git add <file>`) rather than blanket adding (`git add .` or `git add -A`) to avoid accidentally committing temporary files, scratch scripts, or intermediate build artifacts.
