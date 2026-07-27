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
5. **Atomic Commits**: Do not bundle unrelated changes into one massive commit. Break your work down into smaller, logically separated, and atomic commits that correspond to specific tasks or changes (e.g., make one commit for a `refactor` and a separate commit for a `feat`).
