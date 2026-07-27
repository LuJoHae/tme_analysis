---
description: This workflow acts as a system instruction or "Skill" for an AI agent. When asked to review code, the agent MUST follow these step-by-step procedures to strictly enforce the Functional Python Coding Rules.
---

# Agent Workflow: Functional Python Code Review

This workflow acts as a system instruction or "Skill" for an AI agent. When asked to review code, the agent MUST follow these step-by-step procedures to strictly enforce the Functional Python Coding Rules.

## 1. Trigger
This workflow is triggered when the user requests a code review, submits a pull request, or explicitly asks the agent to "check if this code adheres to the functional guidelines."

## 2. Execution Steps

### Step 1: Automated Tooling Verification (Static Analysis)
Before manually inspecting the code, the agent should run or request the results of the following commands in the terminal:
1. `mypy --strict <target_files>`: Ensure there are zero type hinting errors.
2. `flake8 --select=WPS <target_files>` (assuming `wemake-python-styleguide` is installed): Check for mutability, shadowing, and complex logic violations.
*If the tools report errors, immediately report them to the user and halt the review until they are fixed.*

### Step 2: Purity and Expression-Oriented Check
Scan the source code for banned imperative keywords and paradigms:
- **FAIL** if `global`, `nonlocal`, or `del` are used.
- **FAIL** if data structures (lists, dicts, etc.) are mutated in-place (e.g., `.append()`, `.pop()`, `.update()`).
- **WARN** if standard `if/else` statements are used for routing logic where a `match...case` or ternary expression (`x if y else z`) would be more declarative.

### Step 3: Immutability and State
Verify how data is modeled and passed:
- Ensure all data classes use `@dataclass(frozen=True)` or `pydantic.BaseModel` with `frozen=True`.
- Ensure dependencies are passed as arguments (Dependency Injection) rather than instantiated inside pure functions.

### Step 4: Library and Dependency Audit
Check the `import` statements at the top of the file:
- **FAIL** if `import pandas` is found. Recommend `import polars as pl`.
- **FAIL** if `import matplotlib` or `import seaborn` is found. Recommend `import altair as alt`.
- **FAIL** if `Optional` or `None` are used in domain logic. Recommend `Maybe` or `Option` from `returns`.

### Step 5: Error Handling & Monadic Flow
Analyze how the code handles failures:
- **FAIL** if `raise` is used to control business logic or expected failures.
- **PASS** only if errors are returned as values using the `Result` monad.
- Check if pipelines are properly chained. Recommend the use of `@do` notation or `.bind()` for Railway Oriented Programming if the code is manually unpacking `Result` types.

## 3. Reporting
After the review is complete, the agent must generate a final report using the following format:

```markdown
# Functional Code Review Report

## 🟢 Passed Checks
- [List areas where the code adhered well to the rules]

## 🔴 Violations
- [List any strict rule violations, such as side-effects, mutable state, or banned libraries]

## 🟡 Refactoring Suggestions
- [Provide code snippets demonstrating how to rewrite the imperative/impure code into pure, functional code]
```
