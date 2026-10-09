"""Actual block/store/API/CLI execution of the experimental count score."""

import importlib
import importlib.util
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from scipy.stats import false_discovery_control

import epykit as ep
from epykit.methyldata import MethylData


def make_md(root, sites=80, samples=8, sparse=False):
    root = Path(root)
    rng = np.random.default_rng(713)
    pos = np.arange(1, sites + 1, dtype=np.int32) * 10
    ids = [f"s{i}" for i in range(samples)]
    n = rng.integers(8, 61, size=(sites, samples))
    m = rng.binomial(n, rng.beta(3, 7, size=n.shape))
    if sparse:
        n[:4, : samples // 2], m[:4, : samples // 2] = 0, 0
    for i, sample in enumerate(ids):
        p = root / f"sample={sample}" / "chrom=chr1"
        p.mkdir(parents=True)
        pl.DataFrame(
            {"pos": pos, "strand": ["*"] * sites, "N_meth": m[:, i], "coverage": n[:, i]}
        ).write_parquet(p / "part-0.parquet", row_group_size=4096)
    obs = pl.DataFrame(
        {
            "sample_id": ids,
            "treatment": [1] * (samples // 2) + [0] * (samples // 2),
            "group": ["case"] * (samples // 2) + ["control"] * (samples // 2),
            "donor": [f"d{i}" for i in range(samples // 2)] * 2,
            "batch": ([0, 1] * (samples // 4)) * 2,
        }
    )
    md = MethylData(obs=obs, store=str(root))
    md.uns["unite"] = {"type": "union" if sparse else "intersect"}
    return md, pos, m, n


def run(md, **kwargs):
    # Keep actual prior training and fitting, using a bounded test fixture.
    return ep.tl.dmc(md, test="bb_score", score_prior_max_groups=80, tsv=False, **kwargs)


@pytest.mark.parametrize(
    "formula,contrast",
    [
        (None, None),
        ("~ donor + group", "group"),
        ("~ batch + group", "group"),
        ("~ treatment + batch", np.array([0.0, 1.0, 0.0])),
    ],
)
def test_public_api_uses_count_kernel_for_binary_and_nuisance_designs(
    tmp_path, monkeypatch, formula, contrast
):
    assert importlib.util.find_spec("epykit._count_score_store") is not None, (
        "Count store runner missing"
    )
    cs = importlib.import_module("epykit._count_score")
    observed = []
    original = cs.test_counts

    def capture(*args, **kwargs):
        observed.append(np.asarray(args[2]).shape[1])
        return original(*args, **kwargs)

    monkeypatch.setattr(cs, "test_counts", capture)
    md, _, _, _ = make_md(tmp_path / "input")
    run(md, formula=formula, contrast=contrast, score_ci=False, materialize=False)
    assert observed and observed[0] >= 2
    assert md.uns["dmc"]["test_used"] == "bb_score"
    assert md.uns["dmc"]["last_key"] == "dmc_bb_score"
    assert not md.uns["dmc"]["materialized"]
    assert "dmc_glm_contrast" not in md.varm
    assert md.uns["dmc"]["reference"] is None
    assert md.uns["dmc"]["score_reference"] == "normal"
    assert md.uns["dmc"]["inference_status"] == "experimental"
    assert md.dmc_store.test == "bb_score"
    result = md.dmc_store.read_chrom("chr1")
    assert {
        "score_stat",
        "score_variance",
        "rho",
        "rho_uncertainty",
        "score_valid_samples",
        "score_converged",
        "meth_diff_raw",
        "interval_available",
    }.issubset(result.columns)
    assert result["meth_diff_ci_lo"].is_nan().all()
    assert md.dmc_store._manifest["chroms"][0]["prior_sha256"]


def test_block_loader_aligns_sparse_union_without_zero_methylation_substitution(tmp_path):
    assert importlib.util.find_spec("epykit._count_blocks") is not None, "Count blocks missing"
    blocks = importlib.import_module("epykit._count_blocks")
    md, pos, m, n = make_md(tmp_path / "input", sites=37, sparse=True)
    # Entire absent rows, not just zero depth records, must align properly.
    path = Path(md.store) / "sample=s0/chrom=chr1/part-0.parquet"
    frame = pl.read_parquet(path).filter(pl.col("pos") % 30 != 0)
    frame.write_parquet(path)
    expected_m, expected_n = m.copy(), n.copy()
    expected_m[pos % 30 == 0, 0], expected_n[pos % 30 == 0, 0] = 0, 0
    outputs = list(
        blocks.iter_count_blocks(md.store, "chr1", md.obs["sample_id"].to_list(), pos, block_size=7)
    )
    assert max(len(p) for p, _, _ in outputs) <= 7
    np.testing.assert_array_equal(np.concatenate([p for p, _, _ in outputs]), pos)
    np.testing.assert_array_equal(np.vstack([a for _, a, _ in outputs]), expected_m)
    np.testing.assert_array_equal(np.vstack([b for _, _, b in outputs]), expected_n)


def test_duplicate_input_coordinates_are_refused(tmp_path):
    assert importlib.util.find_spec("epykit._count_blocks") is not None, "Count blocks missing"
    blocks = importlib.import_module("epykit._count_blocks")
    md, pos, _, _ = make_md(tmp_path / "input", sites=20)
    path = Path(md.store) / "sample=s0/chrom=chr1/part-0.parquet"
    frame = pl.read_parquet(path)
    pl.concat([frame.head(1), frame]).write_parquet(path)
    with pytest.raises(ValueError, match="duplicate|strictly increasing"):
        list(blocks.iter_count_blocks(md.store, "chr1", ["s0"], pos, block_size=7))


def test_large_chromosome_fitting_is_bounded_and_chunk_invariant(tmp_path, monkeypatch):
    assert importlib.util.find_spec("epykit._count_score_store") is not None, "Count store missing"
    cs = importlib.import_module("epykit._count_score")
    observed = []
    original = cs.test_counts

    def capture(m, n, *args, **kwargs):
        observed.append(m.shape)
        return original(m, n, *args, **kwargs)

    monkeypatch.setattr(cs, "test_counts", capture)
    md, _, _, _ = make_md(tmp_path / "input", sites=33001, samples=8)
    run(md, score_ci=False, materialize=False, score_block_size=32768)
    first = md.dmc_store
    a = first.read_chrom("chr1")
    assert sum(s[0] for s in observed) == 33001 and max(s[0] for s in observed) <= 32768
    observed.clear()
    run(md, score_ci=False, materialize=False, score_block_size=4096)
    second = md.dmc_store
    b = second.read_chrom("chr1")
    assert max(s[0] for s in observed) <= 4096
    for column in ("pvalue", "qvalue", "rho", "meth_diff"):
        np.testing.assert_allclose(a[column], b[column], atol=1e-9, rtol=1e-7, equal_nan=True)
    assert first.manifest["prior_sha256"] == second.manifest["prior_sha256"]


def test_correction_retains_unestimable_sites_in_complete_family(tmp_path):
    md, _, _, _ = make_md(tmp_path / "input", sites=80, sparse=True)
    run(md, score_ci=False, materialize=False)
    frame = md.dmc_store.read_chrom("chr1")
    p = frame["pvalue"].to_numpy()
    assert np.isnan(p[:4]).all()
    expected = false_discovery_control(np.nan_to_num(p, nan=1.0), method="bh")
    np.testing.assert_allclose(frame["qvalue"], expected, rtol=1e-10, atol=1e-12)
    assert md.dmc_store.manifest["hypothesis_count"] == 80
    assert md.dmc_store.manifest["untestable_sites"] >= 4


def test_model_prior_eligibility_and_inference_parameters_invalidate_store_cache(
    tmp_path, monkeypatch
):
    assert importlib.util.find_spec("epykit._count_score_store") is not None, "Count store missing"
    prior_module = importlib.import_module("epykit._dispersion_prior")
    md, _, _, _ = make_md(tmp_path / "input")
    run(md, score_ci=False, materialize=False)
    first = md.dmc_store.path
    run(md, score_ci=False, materialize=False)
    assert md.dmc_store.path == first and md.uns["dmc"]["resumed"]
    run(md, score_ci=True, materialize=False)
    assert md.dmc_store.path != first
    run(md, score_ci=False, materialize=False, formula="~ treatment + batch", contrast="treatment")
    assert md.dmc_store.path != first
    run(md, score_ci=False, materialize=False, min_samples_treatment=3)
    assert md.dmc_store.path != first
    monkeypatch.setattr(prior_module, "ALGORITHM_REVISION", "changed-prior-revision")
    run(md, score_ci=False, materialize=False)
    assert md.dmc_store.path != first


@pytest.mark.parametrize(
    "kwargs",
    [
        {"score_block_size": 0},
        {"score_block_size": 32769},
        {"score_dispersion": "invalid"},
        {"score_reference": "chi2"},
        {"score_prior_max_groups": 0},
        {"smoothing": True},
        {"neighbour_combine": True},
        {"sep_fallback": True},
        {"empirical_fdr": True},
        {"glm_backend": "gpu"},
    ],
)
def test_invalid_score_requests_are_refused_before_output(tmp_path, kwargs):
    from epykit._dmc_config import DMCConfig

    with pytest.raises(ValueError, match="score|bb_score"):
        DMCConfig(test="bb_score", **kwargs).validate()


def test_cli_score_options_parse_and_execute_actual_engine(tmp_path):
    from epykit.cli import build_parser

    md, _, _, _ = make_md(tmp_path / "input")
    sheet = tmp_path / "samples.csv"
    md.obs.write_csv(sheet)
    output = tmp_path / "out.parquet"
    args = build_parser().parse_args(
        [
            "dmc",
            "--methylstore",
            md.store,
            "--samplesheet",
            str(sheet),
            "--treatment-group",
            "case",
            "--control-group",
            "control",
            "--output",
            str(output),
            "--test",
            "bb_score",
            "--formula",
            "~ donor + group",
            "--contrast",
            "group",
            "--score-dispersion",
            "group",
            "--no-score-ci",
            "--score-reference",
            "normal",
            "--score-block-size",
            "32",
            "--no-tsv",
        ]
    )
    args.func(args)
    assert output.is_file()
    frame = pl.read_parquet(output)
    assert "score_stat" in frame.columns and frame["meth_diff_ci_lo"].is_nan().all()


def test_nonzero_contrast_null_is_refused_instead_of_ignored(tmp_path):
    md, _, _, _ = make_md(tmp_path / "input")
    with pytest.raises(ValueError, match="homogeneous"):
        run(md, formula="~ treatment", contrast="treatment = 1")
    assert "dmc" not in md.uns


def test_registry_distinguishes_blockwise_store_dispatch():
    from epykit._dmc_engines import engine_spec
    from epykit.dmc import _ENGINE_RUNNERS

    assert engine_spec("bb_score").blockwise_store
    assert "bb_score" not in _ENGINE_RUNNERS


def test_direct_process_api_accepts_score_and_bounded_controls(tmp_path):
    from epykit.dmc import process_chromosomes_dmc

    md, _, _, _ = make_md(tmp_path / "input")
    store = process_chromosomes_dmc(
        md.store,
        md.treatment_ids,
        md.control_ids,
        test="bb_score",
        score_ci=False,
        score_prior_max_groups=80,
        score_block_size=32,
        return_store=True,
        out_dir=tmp_path / "direct-output",
    )
    assert store.test == "bb_score" and store.path == tmp_path / "direct-output"
    assert "score_stat" in store.read_chrom("chr1").columns
