# WGBS benchmark v4

One workflow: benchmark all six tools once, improve epykit, rerun only epykit,
and compare it with the frozen baseline. No tuning matrices or calibration
experiments are part of the active workflow.

## Run

Use the existing Python 3.10 and R installation on this machine. Local Python
vendor dependencies are loaded automatically; no PYTHONPATH setup is needed.

```bash
cd /scratch/wgbs_benchmark_v4

# Verify Python/R dependencies and inspect the full input/tool selection.
python3 benchmark.py baseline --dry-run

# 1. Run epykit, DSS, methylKit, BSmooth, dmrseq and DMRcate sequentially
#    on both full simulations, then all twelve real samples.
python3 benchmark.py baseline

# 2. Review results/baseline, then edit epykit in /scratch/epykit3.
#    Keep benchmark runners, scoring code and inputs fixed.

# 3. Rerun only epykit on exactly the same inputs, then compare automatically.
python3 benchmark.py epykit --name epykit_v2

# Repeat after another epykit change.
python3 benchmark.py epykit --name epykit_v3
```

`--epykit-source /path/to/epykit/src` selects another epykit checkout.
To rerun all tools, use `baseline --name baseline_v2`; later select it with
`epykit --baseline baseline_v2 --name epykit_v4`. Existing run directories
are never overwritten. A failed caller leaves logs and a failed run manifest;
use a fresh name after fixing the problem. There is no automatic resume of
partial runs. If only comparison failed after the callers completed, rebuild
it without running callers again:

```bash
python3 benchmark.py compare --run epykit_v2
```

Optional small execution check (still runs all six tools):

```bash
python3 benchmark.py baseline --smoke --name smoke
python3 benchmark.py epykit --baseline smoke --name smoke_epykit
```

Smoke runs use simulation pilots and the first 10,000 rows of each real sample.
They are implementation checks, not biological benchmarks. Full runs can take
hours; run them in your usual persistent terminal/session. Each caller has
inherited limits of 48 GiB sampled process-tree RSS and 24 hours.

## Where things live

| Path | Purpose |
|---|---|
| `benchmark.py` | The only user-facing entry point |
| `run_epykit.py`, `run_r_tool.R` | Existing v3 caller settings, preserved |
| `analysis/` | Workflow helpers and real-data gene/interval concordance |
| `src/wgbs_v3/` | Retained scoring, monitoring and provenance code; internal name preserved |
| `data/chr1_3_full_baseline/` | Full chr1–3 signal simulation: 5 vs 5 samples, 659 true regions |
| `data/calibration_v5_full_null_870001/` | Full chr1–3 null simulation: 5 vs 5 samples, no true regions |
| `data/*pilot/` | Small signal/null inputs for smoke checks |
| `data/real_cov/` | Copied GSE64177 coverage inputs: six NI/MTB donor pairs |
| `data/real/` | Created once at first full run: normalized real counts and hashes |
| `results/` | New v4 baseline and epykit revisions only |
| `archive/v3/` | Historical v3 results, workbench runs and supporting context |

Counts, truth and eligibility files are copied, not linked to v3. Simulation
sample sheets point at the v4 copies. Real-data normalization changes decimal
commas in percentages to dots and preserves counts. The redundant GEO tarball,
compressed raw copies and 48 GiB BLUEPRINT download collection are omitted.
They are unnecessary for running the already-generated inputs. V3 is untouched.
See `migration.json` for the copy inventory and selection decisions.

## Read the results

- `results/baseline/simulated/simulation_metrics.csv`: truth-based scores,
  runtime and memory for each simulation/tool.
- `results/baseline/real/concordance/`: real DMR counts, gene annotations,
  interval and gene agreement. Per-tool calls, logs and timings are in `real/TOOL/`.
- `results/epykit_v2/comparison/simulation_metrics.csv`: new epykit plus the
  original five competitors, with source-run and reuse labels.
- `results/epykit_v2/comparison/epykit_changes.csv`: numeric changes from the
  original epykit result (`candidate - baseline`; interpret each metric's direction).
- `results/epykit_v2/comparison/real/concordance/`: agreement recomputed using
  the new epykit calls and original competitor calls.
- `results/epykit_v2/comparison/real_performance.csv`: real runtime/memory,
  explicitly marking reused competitor measurements.

Every run records input and harness hashes and snapshots epykit source.
The workflow checks for source/input changes before accepting a completed run.
Do not edit epykit during a running benchmark. An epykit-only rerun requires
identical inputs and benchmark code; changing caller settings or scoring
requires a new all-tool baseline. Competitor timings remain historical
measurements, not fresh measurements of the machine's current load.

## Interpretation inherited from v3

The signal dataset predates the v4 simulator recorded in the null manifest.
These are fixed development inputs, already inspected in v3; compare each
input across revisions, not signal/null differences as a controlled effect.
No new simulation generation or independent confirmation is implied. Use new
matched simulations and held-out seeds for a later confirmatory study.

Real samples use hg19; simulations use hg38. Real data have no latent truth:
concordance is descriptive, not precision, FDR or biological validation.
The inherited real runners omit a donor-pair model despite the paired design.
The epykit adapter assumes GEO rows already represent strand-collapsed CpGs.
Caller statistical units also differ (CpGs, windows, native regions); in
particular, heuristic regions and epykit asymptotic region rankings should not
be treated as demonstrated region-level FDR control. Those scientific changes
are separate from this workflow cleanup.

## Whole-genome simulation inputs

`python3 prepare_genome.py` prepares **all 22 autosomes (chr1–chr22, GRCh38)**,
excluding X, Y, mitochondrial DNA and alternate contigs. It calibrates against
all 206 verified BLUEPRINT sample pairs at eligible CpGs in the existing
macrophage coordinate catalogue. This is genome-wide catalogue coverage,
not every theoretical CpG in the reference.

The matched signal and null datasets use the same current simulator, seed
20261003, 5 samples per group, nominal 20× mean depth and empirical library-size
variation/dropout. Signal regions have a 0.20 absolute population methylation
difference; the null has no differences. This pair shares the random streams
for the initial coverage draws. At 20×, the inherited methylated/unmethylated
count thinning can produce different realised final depths and eligibility
masks between signal and null; each dataset records and validates its own mask.
Spatial correlation remains disabled, as in the original default simulation.

Generation state: `data/genome_autosomes_preparation/status.json`.
Log: `data/genome_preparation.log`. Each dataset receives a `validation.json`
after all sample counts, truth regions and eligibility masks pass checks.
Completed inputs are `data/genome_autosomes_signal/` and
`data/genome_autosomes_null/`; calibration and source snapshots remain in
`data/genome_autosomes_preparation/`. The one-time generator reads existing
v3 BLUEPRINT sources and the shared reference; generated counts are independent.

Once generation and validation are complete:

```bash
python3 benchmark.py baseline --whole-genome --name genome_baseline
# After improving epykit:
python3 benchmark.py epykit --baseline genome_baseline --name genome_epykit_v2
```

The existing chr1–3 workflow remains the default. Whole-genome runs are larger
and still use the caller limits described above. The new simulations are
additional development inputs, not a claim of independent biological validation.
