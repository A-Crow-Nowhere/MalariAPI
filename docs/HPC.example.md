# HPC Passthrough: Minimal Examples

These examples show the idea: the **same command** can be run locally or submitted to HPC, with simple path tokens.

## Optional aliases (recommended)

```bash
alias m="mapi"
alias mr="mapi rivanna"
```

## Submit a job to HPC

```bash
mr submit --cmd "echo hello > [scratch]/hello.txt"
```

- `[scratch]` expands to the cluster’s scratch path automatically.

## Use your HPC home directory

```bash
mr submit --cmd "ls [home]"
```

- `[home]` expands to your remote home path automatically.

## Check status and pull results

```bash
mr status
mr pull
```

After `mr pull`, results appear locally under:

```
~/MalariAPI/scratch/
```

No manual `scp`, no guessing remote paths.
