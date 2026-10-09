# GSE263850 real-data benchmark

Run: `python3 real_benchmark.py run --name GSE263850_paper_20261009 --profile paper_aligned`.
Status: `python3 real_benchmark.py status --name GSE263850_paper_20261009`.
Rebuild comparison: `python3 real_benchmark.py compare --name GSE263850_paper_20261009`.

GEO SOFT verifies hg38/GRCh38, three heterozygous AKAP11-KO clones versus
three WT replicates, combined-strand counts and a minimum coverage of five.
The adapter rebases obsolete sample-sheet paths, checks BED-to-COV offsets,
validates every supplied COV row and loads 813 published DMRs from table 5.
Raw inputs and earlier benchmark runners/results are preserved.

DSS uses native DMLfit.multiFactor(~group), smoothing with span 500,
DMLtest.multiFactor(coef=2), and callDMR(delta=0,p.threshold=1e-5,
minlen=50,minCG=3,dis.merge=100,pct.sig=.5). DSS 2.58.0 internally clamps
the merge distance to minlen (50); requested and effective values are recorded.

The other tools use the closest mapped profile, documented in each run's
settings. Epykit retains the accepted region correction and uses raw seeds
at 1e-5, smoothing and zero effect minimum. MethylKit uses native CpG tests
and explicit raw-p aggregation instead of its former fixed windows.
BSmooth uses a normal-score-equivalent cutoff, which is not a calibrated
CpG p-value. dmrseq rejects zero candidate cutoffs; this profile uses 0.01,
500 bp smoothing, 100 bp gap and native region q-values. DMRcate uses the
raw-CpG candidate count with its native smoothed rank threshold; its 100 bp
lambda controls both region joining and the kernel. These are approximate
parameter mappings, not claims that the packages implement identical tests.

Full callers run sequentially without imposed memory/time caps. Failures
are recorded while later tools continue. Paper comparisons and a workbook
update after each completed caller. Final gene concordance uses supplied
hg38 RefGene, not the GSE64177 hg19 annotation.

Verification: 12 adapter/geometry/source-compatibility checks passed, and
all six callers plus comparison/annotation finished on the first 10,000
coverage rows per sample. This smoke run tests integration, not performance.
