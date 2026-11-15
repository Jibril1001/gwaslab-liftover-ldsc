"""The synthetic inputs have the layout the tutorials expect."""

from pathlib import Path

import pandas as pd

BBJ_COLUMNS = ["SNP", "CHR", "POS", "REF", "ALT", "Frq", "BETA", "SE", "P", "Dir", "N"]


def test_sumstats_have_bbj_columns(manifest):
    for path in manifest["sumstats"].values():
        df = pd.read_csv(path, sep="\t", nrows=5)
        assert list(df.columns) == BBJ_COLUMNS


def test_t2d_snp_is_chr_pos_and_bmi_snp_is_rsid(manifest):
    t2d = pd.read_csv(manifest["sumstats"]["t2d"], sep="\t", nrows=5)
    bmi = pd.read_csv(manifest["sumstats"]["bmi_female"], sep="\t", nrows=5)
    assert (t2d["SNP"] == t2d["CHR"].astype(str) + ":" + t2d["POS"].astype(str)).all()
    assert bmi["SNP"].str.startswith("rs").all()


def test_t2d_contains_tutorial_unmapped_variants(manifest, truth):
    t2d = pd.read_csv(manifest["sumstats"]["t2d"], sep="\t")
    ids = set(t2d["CHR"].astype(str) + ":" + t2d["POS"].astype(str) + "_" + t2d["REF"] + "_" + t2d["ALT"])
    assert set(truth["unmapped_hg19_to_hg38"]) <= ids


def test_ldscore_files_cover_all_autosomes(manifest):
    ld = manifest["ld"]
    for prefix in [ld["ldscores"], ld["baseline"], ld["weights"]]:
        for c in range(1, 23):
            assert Path(f"{prefix}{c}.l2.ldscore.gz").exists()
            assert Path(f"{prefix}{c}.l2.M_5_50").exists()
    for c in range(1, 23):
        assert Path(f"{ld['baseline']}{c}.annot.gz").exists()


def test_baseline_m_matches_annotation_counts(manifest):
    prefix = manifest["ld"]["baseline"]
    annot = pd.read_csv(f"{prefix}1.annot.gz", sep="\t")
    m = [int(x) for x in Path(f"{prefix}1.l2.M_5_50").read_text().split()]
    assert m == [int(annot["base"].sum()), int(annot["Enhancer"].sum())]


def test_ldcts_lists_each_tissue_with_control(manifest):
    lines = Path(manifest["ld"]["ldcts"]).read_text().splitlines()
    names = [line.split()[0] for line in lines]
    assert names == ["Pancreas", "Liver", "Brain"]
    for line in lines:
        prefixes = line.split()[1].split(",")
        assert len(prefixes) == 2 and prefixes[1].endswith("control.")
