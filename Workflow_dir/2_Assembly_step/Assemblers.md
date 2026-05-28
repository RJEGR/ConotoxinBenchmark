# Assemblers_v3.sh — Usage Documentation

LOSO-aware, mode-aware transcriptome assembly driver with optional TransRate scoring.

`Assemblers_v3.sh` runs one or more de novo assemblers over RNA-seq read pairs
described by a manifest, then optionally scores each resulting assembly with
[TransRate](https://hibberdlab.com/transrate/) in any of four modes. It builds
on `Assemblers_v2.sh` and folds in the read-staging and sanitization machinery
from `Accuracy_v2.sh`.

---

## Table of contents

- [What's new in v3](#whats-new-in-v3)
- [Requirements](#requirements)
- [Synopsis](#synopsis)
- [Options](#options)
- [TransRate modes](#transrate-modes)
- [Input formats](#input-formats)
  - [Manifest](#manifest)
  - [Config file](#config-file)
- [How it works](#how-it-works)
- [Output layout](#output-layout)
- [Examples](#examples)
- [Running under SLURM](#running-under-slurm)
- [Read sanitization](#read-sanitization)
- [Checkpointing](#checkpointing)
- [Error handling](#error-handling)
- [Migration from v2](#migration-from-v2)
- [Troubleshooting](#troubleshooting)

---

## What's new in v3

| Feature | v2 | v3 |
|---|---|---|
| TransRate step | Always run (reference-only) | Selectable via `-m` (`none`/`reference`/`reads`/`both`) |
| `stage_reads()` | Inline in main loop | Dedicated function (ported from `Accuracy_v2.sh`) |
| `sanitize_reads()` | — | seqkit `sana` + `pair` sanitization, with `-S` opt-out |
| Read-based scoring | Not possible | `-m reads` / `-m both` |
| Per-mode scratch dirs | — | `<assembly>_<mode>_dir`, `<assembly>_<mode>.log` |
| Meaning of `-m` | Memory in GB | **TransRate mode** (memory moved to `-C`) |

> **Breaking change:** the `-m` flag changed meaning between v2 and v3. See
> [Migration from v2](#migration-from-v2).

---

## Requirements

- **Bash** (the script uses `set -eo pipefail`, arrays, and `getopts`).
- **Assembler wrappers** — either on `PATH`, or as `*.sh` files inside a
  `run_assemblers_dir/` directory next to where the script is invoked. If that
  directory exists, the script makes the wrappers executable and prepends it to
  `PATH` automatically.
- **TransRate** — required unless running with `-m none`. The script sets up a
  site-specific environment (Ruby 2.2.0 and `transrate` under
  `/LUSTRE/apps/bioinformatica/...`); adjust `setup_transrate_env()` for other
  hosts.
- **seqkit** — required for read sanitization (`-m reads`/`both` without `-S`).
  Available from <https://bioinf.shenwei.me/seqkit/>. If it is missing the
  script warns and proceeds with unsanitized reads.

---

## Synopsis

```text
Assemblers_v3.sh -s <Manifest> -c <Config>
                 [-m none|reference|reads|both] [-t <CPUs>] [-C <MEM_GB>] [-S]
```

`-s` and `-c` are mandatory. All other options have defaults.

---

## Options

| Flag | Argument | Default | Description |
|---|---|---|---|
| `-s` | manifest path or `all` | — *(required)* | Manifest TSV. `all` processes every `*.txt` manifest in the current directory. |
| `-c` | config path | — *(required)* | Assembler config file (one `tool=command_template` per line). |
| `-m` | mode | `reference` | TransRate mode for Step 2. One of `none`, `reference`, `reads`, `both`. |
| `-t` | integer | `$SLURM_CPUS_ON_NODE`, else `20` | CPU threads. |
| `-C` | integer | `$SLURM_MEM_PER_NODE`/1024, else `100` | Memory in GB. Informational — TransRate has no memory cap, but `MEM` is exported for config templates / user hooks. |
| `-S` | *(none)* | sanitization **on** | Skip the seqkit `sana` + `pair` read-sanitization step for `reads`/`both` modes. |
| `-h` | *(none)* | — | Print help and exit. |

> **Note:** memory is `-C` (uppercase) in v3 because lowercase `-c` is the
> config-file flag.

---

## TransRate modes

The `-m` value selects whether and how Step 2 runs. For every mode except
`none`, the resulting `transrate` call is:

| Mode | Behaviour | TransRate call |
|---|---|---|
| `none` | Assemble only; **Step 2 skipped entirely**. | *(not run)* |
| `reference` | Reference-only metrics (assembly vs held-out ground truth). | `transrate --assembly <query> --reference <REFERENCE> --output <tr_path> --threads <CPU>` |
| `reads` | Read-only metrics (left/right vs assembly). | `transrate --left <fwd> --right <rev> --assembly <query> --output <tr_path> --threads <CPU>` |
| `both` | Reference **and** reads. | `transrate --left <fwd> --right <rev> --assembly <query> --reference <REFERENCE> --output <tr_path> --threads <CPU>` |

- `reference` reproduces the behaviour of `Assemblers_v2.sh`.
- `reads` and `both` require a sanitized (or `-S`) read pair; the same pair
  staged for the assemblers is reused — reads are never staged twice.
- `none` is useful when you intend to score later with `Accuracy_v2.sh`.

---

## Input formats

### Manifest

A tab-separated file with four columns and no required header:

```text
sample_id    fwd_fq              rev_fq              ref_fasta
sampleA      /path/A_R1.fq       /path/A_R2.fq       /path/test_M_superfamily.fasta
```

- Blank lines and lines beginning with `#` are ignored.
- All rows in a fold should share the same `ref_fasta`; if they differ, the
  script warns and uses **row 1's** reference.
- **LOSO manifests** from `Simulate_rnaseq_loso.sh` (e.g.
  `test_M_superfamily_200x_PE_samples.txt`) are supported directly. A trailing
  `_samples` in the basename is stripped to form the **fold stem**, e.g.
  `test_M_superfamily_200x_PE`.
- Single-row manifests trigger an I/O optimization: reads are **symlinked**
  rather than concatenated.

### Config file

One assembler per line, in `tool=command_template` form. Lines starting with
`#` and blank lines are ignored. Templates are `eval echo`-expanded with these
in-scope variables:

```text
$forward_fq   $reverse_fq   $OUTDIR   $CPU   $MEM   $FASTA_DIR   $REFERENCE
```

Example config line:

```text
trinity=run_trinity.sh $forward_fq $reverse_fq $OUTDIR $CPU $MEM
```

Each assembler is expected to write its final assembly to
`<FASTA_DIR>/<stem>_<tool>.fa`.

---

## How it works

For each manifest the script:

1. **Derives the fold stem** from the manifest name (strips `_samples`).
2. **Loads manifest rows** into `sample_id` / `fwd` / `rev` / `ref` arrays.
3. **Stages the read pair** as `<stem>_concat_PE1.fq` / `<stem>_concat_PE2.fq`
   — symlinked for single-row manifests, concatenated otherwise.
4. **Sanitizes the read pair** (only if `-m reads`/`both` and not `-S`) using
   seqkit `sana` → `pair`.
5. For each tool in the config:
   - **Step 1 — assemble.** Runs the expanded command template. Failed or
     empty assemblies are moved to `issues_dir/` and the batch continues.
   - **Step 2 — TransRate.** Skipped if `-m none`; otherwise runs
     `run_transrate()` in the selected mode. `contigs.csv` outputs are
     collected into `transrate_contigs_dir/`.
   - Writes a checkpoint and removes the assembler's working directory.
6. **Cleans up** the staged read pair (symlinks or temp/sanitized files).

---

## Output layout

Relative to the working directory:

```text
<stem>_FASTA_DIR/
    <stem>_<tool>.fa            # final assembly per tool
chkp_dir/
    1_<tool>_<stem>.chkp        # per-(tool, fold) checkpoint
transrate_tmp_dir/
    <assembly>_<mode>.log       # TransRate run log
    <stem>_seqkit.log           # seqkit sanitization log
transrate_contigs_dir/
    .../contigs.csv             # collected TransRate contig metrics
issues_dir/
    <stem>_<tool>_dir/          # diverted failed / empty assemblies
```

Per-mode scratch directories (`<assembly>_<mode>_dir`) are created under
`transrate_tmp_dir/` during a run and removed afterward, so reference/reads/both
runs of the same assembly never collide.

---

## Examples

Assemble and score against the LOSO ground-truth reference (v2-equivalent):

```bash
./Assemblers_v3.sh -s test_M_superfamily_200x_PE_samples.txt -c assemblers.cfg
```

Assemble only — no scoring (score later with `Accuracy_v2.sh`):

```bash
./Assemblers_v3.sh -s test_M_superfamily_200x_PE_samples.txt -c assemblers.cfg -m none
```

Read-based scoring with explicit resources:

```bash
./Assemblers_v3.sh -s manifest.txt -c assemblers.cfg -m reads -t 24 -C 100
```

Reference + reads, skipping read sanitization:

```bash
./Assemblers_v3.sh -s manifest.txt -c assemblers.cfg -m both -S
```

Batch-process every manifest in the directory:

```bash
./Assemblers_v3.sh -s all -c assemblers.cfg -m reference
```

---

## Running under SLURM

The script carries `#SBATCH` directives (`--mem=100GB`, `--ntasks-per-node=24`,
`-t 6-00:00:00`) and can be submitted directly:

```bash
sbatch Assemblers_v3.sh -s all -c assemblers.cfg -m both
```

When running inside a SLURM allocation:

- `-t` defaults to `$SLURM_CPUS_ON_NODE`.
- `-C` defaults to `$SLURM_MEM_PER_NODE / 1024` (SLURM reports MB).
- Explicit `-t` / `-C` always override the SLURM-derived values.

---

## Read sanitization

For `-m reads` and `-m both`, the staged pair is sanitized before scoring to
prevent TransRate's `Unmatched read IDs ... Use the -I option to ignore this`
failure:

1. **`seqkit sana`** — repairs or drops broken single-line FASTQ records.
2. **`seqkit pair`** — keeps only reads whose IDs appear in **both** R1 and R2.

Sanitization is **in-place** on the staged files. When the staged files are
symlinks (single-row LOSO manifests) the script deletes the symlink and writes
a fresh file — the original source reads on disk are **never mutated**.

If any sanitization step fails or produces empty output, the script warns and
falls back to the unsanitized reads rather than aborting. Pass `-S` to skip
sanitization entirely.

---

## Checkpointing

Each `(tool, fold)` pair gets a checkpoint file `chkp_dir/1_<tool>_<stem>.chkp`
created after a successful assembly. On re-run, any tool/fold with an existing
checkpoint is skipped, so interrupted batches resume cheaply.

> Checkpoints track **assembly completion**, not the TransRate mode. To force
> re-assembly, delete the relevant `chkp_dir/*.chkp` files.

---

## Error handling

- The script runs under `set -eo pipefail`.
- An assembler that exits non-zero, or produces a missing/empty
  `<FASTA_DIR>/<stem>_<tool>.fa`, has its working directory moved to
  `issues_dir/`; the batch continues with the next tool.
- A failed TransRate run logs to `transrate_tmp_dir/<assembly>_<mode>.log` and
  is reported, but does not abort the batch.
- Invalid `-m` values, a missing config file, or a missing manifest cause an
  immediate exit with a usage message.

---

## Migration from v2

| Concern | v2 | v3 |
|---|---|---|
| `-m` | Memory in GB | TransRate mode |
| Memory flag | `-m` | `-C` |
| Default Step 2 | Always reference-only | `-m reference` (same result) |

**Action required:** existing job scripts that passed `-m <GB>` must be updated.
A call like `-m 100` will now fail with `invalid -m '100'`. Replace it with
`-C 100`, and add an explicit `-m <mode>` if a non-default mode is wanted.

```bash
# v2
./Assemblers_v2.sh -s manifest.txt -c cfg -t 24 -m 100

# v3 equivalent
./Assemblers_v3.sh -s manifest.txt -c cfg -t 24 -C 100 -m reference
```

---

## Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| `invalid -m '100'` | Using v2's `-m` for memory. Use `-C 100`. |
| `seqkit not on PATH; skipping pair sanitization` | Install seqkit, or accept unsanitized reads, or pass `-S`. |
| TransRate `Unmatched read IDs` | Sanitization was skipped or failed; ensure seqkit is available and re-run without `-S`. |
| `no contigs.csv produced` | TransRate ran but emitted nothing; inspect `transrate_tmp_dir/<assembly>_<mode>.log`. |
| Assembly diverted to `issues_dir/` | Assembler exited non-zero or produced an empty FASTA; inspect the moved working directory. |
| `run_assemblers_dir not found` note | Wrappers must be on `PATH`; otherwise place them in `run_assemblers_dir/`. |
| `manifest has mixed reference columns` | Rows disagree on `ref_fasta`; row 1's reference is used. |
