# Git Wrappers: Minimal Examples

These examples show how MAPI makes common Git workflows easier and safer for non-experts.

## Upload your work safely

```bash
mapi git upload
```

What this is meant to do:
- respect `.gitignore`
- avoid committing large environments/data
- push changes with safe defaults

## Pull updates without breaking local work

```bash
mapi git update
```

What this is meant to do:
- pull upstream changes
- avoid overwriting local data
- keep environments untouched

## Switch branches safely

```bash
mapi git switch dev
```

This avoids common beginner pitfalls (detached HEAD, accidental overwrites).
