---
name: update-template
description: Update this ingest repo to the latest copier template version, resolving known conflicts
version: 1.0.0
triggers:
  - "update template"
  - "copier update"
  - "update to latest template"
---

# Update to Latest Template

Update this ingest repo to the latest koza-ingest-template version using copier, with deterministic conflict resolution.

**Use when:** The template has been updated (new CI patterns, justfile improvements, dependency changes) and you want to pull those changes into this repo.

**Do NOT use when:**
- This repo has no `.copier-answers.yml` (use the `migrate-standalone-ingest` skill first)
- You want to create a new ingest (use the `create-koza-ingest` skill)

---

## Phase 1: Pre-flight

### 1.1 Verify copier-managed repo

Check for `.copier-answers.yml` at repo root. If missing, stop and tell the user this repo needs migration first.

### 1.2 Record current state

Capture these values before making any changes:

```bash
git branch --show-current        # Current branch
git status --porcelain            # Dirty files
grep '_commit:' .copier-answers.yml  # Current template commit
```

### 1.3 Capture local justfile

**Read the entire justfile and save these values in memory** — you will need them after copier update:

1. **TRANSFORMS variable**: The value after `TRANSFORMS :=` (e.g., `"gene disease phenotype"`)
2. **`download:` recipe**: The full recipe body (e.g., `uv run downloader download.yaml`)
3. **`run:` recipe**: The full dependency chain (e.g., `run: download postdownload transform-all postprocess`)
4. **`transform-all:` dependencies**: What it depends on (e.g., `transform-all: download` or `transform-all: preprocess`)
5. **Custom recipes**: Any of these if they exist — save the full recipe including body:
   - `preprocess:`
   - `postprocess:`
   - `postdownload:`
   - Any other non-standard recipes

### 1.4 Report and confirm

Tell the user:
```
Current branch: <branch>
Current template commit: <commit>
Dirty files: <list or none>
Custom justfile recipes: <list>
```

---

## Phase 2: Stash if Needed

If `git status` shows uncommitted tracked changes:

```bash
git stash push -m "pre-template-update stash"
```

Record that a stash was created. Untracked files are fine — they won't interfere with copier.

---

## Phase 3: Run Copier Update

```bash
copier update --trust --defaults
```

After running, check if the template commit changed:

```bash
grep '_commit:' .copier-answers.yml
```

If the commit is the same as before, the repo is already up to date. Skip to Phase 7 (restore stash if needed).

---

## Phase 4: Resolve Conflicts

### 4.1 Check for conflict markers

```bash
grep -rn "^<<<<<<<\|^=======\|^>>>>>>>" --include="*.yaml" --include="*.yml" --include="justfile" --include="*.toml" .
```

Even if no conflict markers are found, still verify the justfile and workflows (copier may have cleanly overwritten local customizations).

### 4.2 Resolve justfile

**Strategy: Local-first rebuild.** The justfile has repo-specific content (TRANSFORMS, download command, pipeline dependencies) that must be preserved.

1. **Check for conflict markers** in justfile
2. If conflicts exist, resolve each hunk:
   - Hunks in `download:` recipe → **keep "before updating"** (local version)
   - Hunks in `run:` recipe → **keep "before updating"** (local version)
   - Hunks in custom recipes (preprocess/postdownload/etc.) → **keep "before updating"**
   - Hunks in standard recipes (install/test/lint/format) → **keep "after updating"** (template version)

   To keep "before updating": `sed -i '' '/^=======$/,/^>>>>>>> after updating$/d; /^<<<<<<< before updating$/d' justfile`
   To keep "after updating": `sed -i '' '/^<<<<<<< before updating$/,/^=======$/d; /^>>>>>>> after updating$/d' justfile`

3. **Verify TRANSFORMS** is not empty (unless it was empty before). If copier overwrote it with `""`, restore from Phase 1 capture.

4. **Verify download recipe** matches Phase 1 capture. If copier replaced it with `uv run python scripts/download.py`, restore the local version.

5. **Verify run recipe** dependencies match Phase 1 capture.

6. **Verify custom recipes** (preprocess, postprocess, postdownload) still exist if they did before. Restore from Phase 1 capture if missing.

7. **Remove duplicate recipes**: Check for duplicate `run:` (or any other) recipe definitions:
   ```bash
   grep -c '^run:' justfile
   ```
   If count > 1, the template added a duplicate. Remove the **second** occurrence (the template-added one, typically `run: transform-all test` with comment "Run full pipeline: install, download, transform, test"). Remove the comment, group annotation, and recipe definition (3 lines).

### 4.3 Resolve workflow files

**test.yaml** — Take the template version (standardized CI with `just run`):

If conflict markers exist:
```bash
sed -i '' '/^<<<<<<< before updating$/,/^=======$/d; /^>>>>>>> after updating$/d' .github/workflows/test.yaml
```

**release.yaml** — Conditional:
- If the local release.yaml has `schedule:`, `workflow_dispatch:`, custom `env:`, or `softprops/action-gh-release` with a custom body → **keep local** (it has repo-specific release logic)
- If it closely matches the template (simple trigger, no custom body) → **take template**

For keeping local:
```bash
sed -i '' '/^=======$/,/^>>>>>>> after updating$/d; /^<<<<<<< before updating$/d' .github/workflows/release.yaml
```

For taking template:
```bash
sed -i '' '/^<<<<<<< before updating$/,/^=======$/d; /^>>>>>>> after updating$/d' .github/workflows/release.yaml
```

### 4.4 Resolve other files

- **pyproject.toml**: If conflicts, take template version. Then verify custom dependencies (e.g., `kghub-downloader`, `duckdb`, `cyvcf2`) are still listed. Restore any that were lost.
- **.copier-answers.yml**: Always take the copier result (no conflicts expected).
- **.gitignore**: Take template version (additive changes are safe).

---

## Phase 5: Fix Unrendered Jinja Tags

Copier sometimes fails to render Jinja tags in workflow files. Check and fix:

```bash
grep -rn '{% raw %}\|{% endraw %}' .github/workflows/
```

If found, replace patterns like:
```
{% raw %}${{ matrix.python-version }}{% endraw %}
```
with:
```
${{ matrix.python-version }}
```

Use sed or manual edit to fix each occurrence.

---

## Phase 6: Verification

Run these checks before committing:

```bash
# No conflict markers remain
grep -rn "^<<<<<<<\|^=======\|^>>>>>>>" . --include="*.yaml" --include="*.yml" --include="justfile" --include="*.toml"

# No unrendered Jinja tags
grep -rn '{% raw %}\|{% endraw %}' .github/workflows/

# Justfile parses correctly
just --list

# .copier-answers.yml has new commit
grep '_commit:' .copier-answers.yml
```

Also verify by inspection:
- TRANSFORMS variable is populated correctly
- `download:` recipe has the correct command for this repo
- `run:` recipe has the correct dependency chain
- Custom recipes (preprocess, postprocess, etc.) are present if they should be

---

## Phase 7: Commit and Push

**Always ask the user before committing.**

List the changed files:
```bash
git status
```

Stage specific files (never `git add .`):
```bash
git add .copier-answers.yml justfile .github/workflows/test.yaml .github/workflows/release.yaml
# Add any other changed files as appropriate
```

Commit:
```bash
git commit -m "chore: update to latest ingest template"
```

Push to current branch:
```bash
git push origin $(git branch --show-current)
```

If a stash was created in Phase 2:
```bash
git stash pop
```

If stash pop has conflicts, warn the user and help resolve them.

---

## Conflict Resolution Quick Reference

| File/Section | Keep Local | Take Template | Notes |
|-------------|-----------|---------------|-------|
| justfile `TRANSFORMS` | Yes | | Template has empty string |
| justfile `download:` | Yes | | Custom download commands per repo |
| justfile `run:` | Yes | | Custom pipeline dependency chains |
| justfile `preprocess:/postprocess:/postdownload:` | Yes | | Pipeline-specific custom recipes |
| justfile `transform-all:` dependencies | Yes | | May depend on preprocess instead of download |
| justfile `transform-all:` body | | Yes | Standardized iteration loop |
| justfile `install:/test:/lint:/format:/clean:` | | Yes | Keep standardized |
| justfile `setup:/_git-init/_git-add` | | Yes | Only used at project init |
| `.github/workflows/test.yaml` | | Yes | Standardized CI |
| `.github/workflows/release.yaml` (custom) | Yes | | Repo-specific release workflows |
| `.github/workflows/release.yaml` (simple) | | Yes | Matches template pattern |
| `pyproject.toml` | | Yes* | *Restore custom dependencies |
| `.copier-answers.yml` | | Yes | Always from copier |
| Duplicate `run:` recipe | Remove 2nd | | Copier merge artifact |

---

## Troubleshooting

### "Updating is only supported in git-tracked templates"
The `_src_path` in `.copier-answers.yml` points to a non-git source (e.g., old cookiecutter). This repo needs full migration first — use the `migrate-standalone-ingest` skill.

### copier update prompts for input
Use `--defaults` flag to skip prompts. If copier still prompts, check that `.copier-answers.yml` has all required fields.

### Tests fail after update
Common causes:
- **Koza 2.x validation errors** (`Unexpected keyword argument` for reader/writer `name`): The transform YAML configs need migration to Koza 2.x format (remove `name` fields from reader/writer sections).
- **Missing `dev` dependency group**: The `pyproject.toml` needs a `[dependency-groups]` section.
- **Download failures in CI**: The `just run` recipe tries the full pipeline including download. If the data source requires secrets or authentication, configure them in the GitHub repo's secrets settings.
