---
trigger: always_on
---

# Three-Branch Worktree Workflow (`<type>/<slug>` ↔ `review/<slug>` → `main`)

This repository uses a strict three-branch `git worktree` workflow to keep agent execution isolated, user review conflict-free, and `main` always stable.

## 1. Branch Roles & Ownership

1. **`main`** — Protected integration branch.
   - **Agents must NEVER** switch to, commit directly to, merge into, or reset `main` unless the user explicitly asks to land work on `main`.
2. **`<type>/<slug>`** (e.g. `feat/bagaev-mfp`, `fix/dropdown-sort`) — Agent working branch inside its dedicated `git worktree`.
   - Each agent works and commits **exclusively** on its own assigned working branch in its own worktree.
3. **`review/<slug>`** — User review and polish branch in the main workspace.
   - **Agents must NEVER** check out, switch to, commit to, reset, or force-push `review/<slug>`.

## 2. Mandatory Agent Workflow Loop

1. **Verify Isolation Before Editing**:
   - Confirm your current directory is a worktree checked out on `<type>/<slug>` (not `main`, not `review/*`, and not detached `HEAD` unless instructed).
2. **Pull User Review Changes Before Resuming**:
   - Whenever the user asks you to continue after reviewing (or mentions commits on `review/<slug>`), first merge the review branch into your current working branch:
     ```bash
     git merge review/<slug>
     ```
   - If merge conflicts occur, resolve them carefully preserving both the user's review edits and your intended behavior, run the relevant tests/checks, and commit the merge before continuing.
3. **Commit Incrementally in Small Conventional Commits**:
   - Follow [git-conventional-commits.md](file:///Users/halu/Code/tme_analysis/.agents/rules/git-conventional-commits.md) (`https://www.conventionalcommits.org/en/v1.0.0/`).
   - Never leave uncommitted changes at the end of a turn or bundle large multi-layer edits into a single commit.
4. **Clean Handoff for Review**:
   - Verify all relevant tests and linters pass after your final commit.
   - Report your exact working branch (`<type>/<slug>`), the commit list (`git log --oneline`), and the verification results so the user can run `greview <slug>` in the main workspace.
