# Default-path audit — 2026-10-07

The committed calibration work targets optional permutation inference. A later uncommitted patch targets the default chain-merge region caller. Neither has a completed full-genome candidate benchmark in the inspected development results directory.

## Evidence

- Benchmark default: `run_epykit.py:32` sets `empirical_fdr=False`. The profile switch at line 38 enables permutations only for explicit permutation profiles.
- Commit `1f67ddf` (2026-10-06 22:52 +0200) changes permutation aggregation, population handling and documentation in `src/epykit/dmr.py`; its tests target that optional path. Between preserved baseline `0566519` and committed HEAD `ce7ef4c`, `dmr.py` is the only changed package source file.
- The saved permutation pilot is documented in `reports/README.md:7`: 99 shuffled scans on chr22. Its improvements compare inference modes and cannot establish improvement in the unchanged default path.
- The working tree adds `_region_screen.py`, calls its geometry-window screen from `src/epykit/dmr.py:1010`, updates combined q-values at line 1148, and adds `region_scoring` metadata in `tl.py`. These changes are uncommitted; the helper and regression test are untracked.
- `scripts/check_position_screen.py:14` reads frozen full-genome CpG test results, but line 30 reads candidates from a separate chr22-only recomputation pilot. It evaluates a standalone filter and does not invoke the modified package caller. Its 32 signal calls / 24 matches / 75% precision must not be attributed to the current package implementation.

## Fresh checks

Imported package confirmed as `/scratch/epykit-calibration/src/epykit/dmr.py`.

Using identical frozen chr22 DMC parquet inputs and the original chain-merge parameters, calling the default region caller and applying q <= .05 produced:

| Source | Signal calls | Null calls |
|---|---:|---:|
| Preserved baseline 0566519 | 47 | 18 |
| Committed HEAD ce7ef4c | 47 | 18 |
| Current uncommitted working tree | 19 | 0 |

Baseline and committed-HEAD result tables were exactly equal on both inputs. Working-tree signal intervals differed from the saved prototype: 19 versus 32 calls, with 9 package-only and 22 prototype-only intervals. Null outputs were empty in both. These are chr22 extraction checks, not genome-wide benchmarks or held-out validation; the regional BH testing family is restricted to chr22 in these checks.

Fresh targeted regression run: `tests/test_region_position_screen.py`, `tests/test_dmr_chain_merge.py`, `tests/test_chain_merge_empirical_fdr.py`, and `tests/test_region_permutation_population.py`: **65 passed, 5 warnings**, 9.74 seconds.

An older broader test log (`results/default_path_tests.log`) records 51 passes and one failure caused by missing `bioframe`; the targeted passing checks do not establish a fully passing suite.

## Implication

The user's recollection is supported: the committed fix and main saved development report focused on optional permutations, which the original benchmark does not use. The default-path working-tree change is a separate, incompletely evaluated candidate. CpG calibration remains unresolved by this region-only patch. Evaluate the actual package under unchanged default benchmark settings before using prototype performance claims or launching long comparative runs.
