# Epykit comprehensive performance review

Read `report.md` for the complete technical report. `evidence.xlsx` contains the reviewed comparisons and diagnostic data; `figures/` contains four standalone PNG/SVG figures. Tables distinguish original epykit, the accepted region correction and experimental BB-F.

The main six-tool comparisons cover the completed autosomal simulations and GSE64177 real data. `gse263850_snapshot.json` freezes a separate partial comparison of verified completed GSE263850 callers; pending tools are unavailable, not zero. This is a report of available evidence, not a rerun or implementation of the recommendations.

## Reproduce

From `/scratch/wgbs_benchmark_v4`:

```bash
python3 analysis/epykit_final_report_20261009/analyze.py
python3 analysis/epykit_final_report_20261009/candidate_ceiling.py
python3 analysis/epykit_final_report_20261009/profile_bb.py
python3 analysis/epykit_final_report_20261009/collect_latest_real.py
python3 analysis/epykit_final_report_20261009/build_report.py
python3 analysis/epykit_final_report_20261009/plot_report.py
python3 analysis/epykit_final_report_20261009/verify_report.py
```

`analyze.py` verifies arithmetic, original/corrected CpG byte equality, source hashes and full-family region correction. `candidate_ceiling.py` independently reconstructs reciprocal overlap matches from truth/eligibility and cached pre-q candidates. `profile_bb.py` measures bounded warm compiled kernels; it does not claim whole-genome CPU profiling or a statistically equivalent fixed-rho replacement. `collect_latest_real.py` verifies completion and independently reconstructs GSE263850 50% paper matches before freezing its evolving evidence.

`report.template.md` owns prose and figure/table placements. `build_report.py` fills the tables and creates the workbook and reviewed snapshot. `verify_report.py` checks deliverable coverage, recorded results, supporting diagnostics and artifacts. Re-running collection later updates the evidence cutoff; review the prose before refreshing a partial study.

## Interactive version

`app/` contains a supplementary editable report draft prepared with the Data shared runtime. It is not a delivered or verified app. The standard build failed with `spawnSync /usr/bin/node EPERM` while launching the canonical authored-file verifier. No protected runtime was changed or verification bypassed. Its authoring state is paused. The delivered Markdown, workbook and scientific figures are independent of this draft.
