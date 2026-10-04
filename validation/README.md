# Verification — 2026-10-02

- Five orchestration tests passed: all-tool baseline followed by epykit-only
  execution; reuse labels and baseline preservation; changed-input rejection;
  failed-run rejection; overwrite protection; source-change rejection.
- Python compilation and parsing of both R scripts passed.
- Baseline dry run passed, including Python and all six tools' R dependencies.
- Actual epykit execution passed on the signal pilot and on 10,000 rows per
  real sample; actual DSS execution passed on the same small real input.
- Real hg19 annotation and two-tool concordance passed using those actual
  epykit/DSS outputs, including an empty epykit DMR list.
- SHA-256 contents of all 88 selected input files (excluding intentionally
  rebased sample sheets) match v3. See `copied_data_hashes.json`.
- Historical results, workbench and development copies match original file
  counts, paths and sizes. Historical output contents were not fully rehashed.

Execution-check outputs are under `smoke/`. They are not a baseline and are
never selected by `benchmark.py`. The full all-tool benchmark has not run.
Warnings observed during small checks concern an epykit deprecated column,
R data.table copying, and unavailable timedatectl in the sandbox; all callers
and annotation exited successfully.
