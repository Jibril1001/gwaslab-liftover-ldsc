"""Liftover tutorial: https://cloufield.github.io/gwaslab/harmonization_liftover/"""

import pandas as pd
import pytest

import liftover_walkthrough as lw
from conftest import T2T_CHAIN


@pytest.fixture(scope="module")
def lifted(manifest, tmp_path_factory):
    """Tutorial run: 5,000 rows, hg19 -> hg38 (-> CHM13 when the chain is present)."""
    chain = T2T_CHAIN if T2T_CHAIN.exists() else None
    return lw.run(manifest["sumstats"]["t2d"], chain, nrows=5000, outdir=tmp_path_factory.mktemp("lift"))


@pytest.fixture(scope="module")
def panel(manifest):
    return pd.read_csv(manifest["panel"], sep="\t")


def test_loads_tutorial_row_count(lifted):
    assert lifted["steps"][0]["variants"] == 5000


def test_hg19_to_hg38_sets_build(lifted):
    step = lifted["steps"][1]
    assert step["step"] == "hg19 -> hg38"
    assert step["build"] == "38"


def test_tutorial_unmapped_variants_are_removed(lifted, truth):
    coords = lifted["coords"]
    for vid in truth["unmapped_hg19_to_hg38"]:
        chrom_pos = vid.split("_")[0]
        row = coords[coords["SNPID"] == chrom_pos]
        assert len(row) == 1, vid
        assert row["POS_hg38"].isna().all(), f"{vid} should be unmapped"
    assert lifted["steps"][1]["unmapped"] >= len(truth["unmapped_hg19_to_hg38"])


def test_hg38_positions_match_dbsnp_hapmap3(lifted, panel):
    """Lifted positions agree with dbSNP151 hg38 coordinates of the same SNPs."""
    coords = lifted["coords"].dropna(subset=["POS_hg38"])
    panel = panel.assign(SNPID=panel["CHR"].astype(str) + ":" + panel["BP"].astype(str))
    merged = coords.merge(panel, on="SNPID")
    assert len(merged) > 4900
    agree = (merged["POS_hg38"].astype(int) == merged["POS38"]).mean()
    assert agree > 0.995, f"only {agree:.4f} of positions match dbSNP hg38"


def test_status_code_build_digits_change(lifted):
    coords = lifted["coords"].dropna(subset=["STATUS_hg38"])
    assert coords["STATUS_hg38"].astype(int).astype(str).str[:2].eq("38").all()


def test_statistics_unchanged_by_liftover(manifest, lifted):
    raw = pd.read_csv(manifest["sumstats"]["t2d"], sep="\t", nrows=5000)
    raw["SNPID"] = raw["SNP"]
    out = (
        lifted["sumstats"]
        .data[["SNPID", "BETA", "SE", "P"]]
        .merge(raw[["SNPID", "BETA", "SE", "P"]], on="SNPID", suffixes=("", "_raw"))
    )
    assert len(out) == len(lifted["sumstats"].data)
    for col in ["BETA", "SE", "P"]:
        assert (out[col] - out[f"{col}_raw"]).abs().max() < 1e-12


def test_round_trip_hg38_to_hg19_restores_positions(manifest):
    sumstats = lw.load_t2d(manifest["sumstats"]["t2d"], nrows=2000)
    original = sumstats.data.set_index("SNPID")["POS"].copy()
    sumstats.liftover(from_build="19", to_build="38", verbose=False)
    sumstats.liftover(from_build="38", to_build="19", verbose=False)
    assert sumstats.build == "19"
    back = sumstats.data.set_index("SNPID")["POS"]
    same = (back == original.loc[back.index]).mean()
    assert same > 0.999


@pytest.mark.skipif(not T2T_CHAIN.exists(), reason="run scripts/download_data.sh for the T2T chain")
def test_hg38_to_chm13_with_downloaded_chain(lifted):
    step = lifted["steps"][2]
    assert step["step"] == "hg38 -> T2T-CHM13"
    assert step["build"] == "13"
    assert step["variants"] > 0.98 * lifted["steps"][1]["variants"]
    # First tutorial-comparable SNP on chr1 shifts to the CHM13 coordinate system.
    first = lifted["coords"].dropna(subset=["POS_chm13"]).iloc[0]
    assert first["POS_chm13"] != first["POS_hg38"]


def test_outputs_written(lifted):
    outdir = lifted["outdir"]
    steps = pd.read_csv(outdir / "liftover_steps.tsv", sep="\t")
    assert list(steps["variants"]) == [s["variants"] for s in lifted["steps"]]
    assert (outdir / f"t2d_lifted_build{lifted['sumstats'].build}.tsv.gz").exists()
