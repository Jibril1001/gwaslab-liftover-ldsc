"""LDSC tutorial: https://cloufield.github.io/gwaslab/ldsc_in_gwaslab/

The synthetic z-scores are drawn from the LDSC model with known parameters, so
each GWASLab estimate is checked against the simulated truth (within 3 SE).
"""

import pytest

import ldsc_walkthrough as lw

Z = 3.0  # tolerance in standard errors


@pytest.fixture(scope="module")
def ldsc(manifest, tmp_path_factory):
    s = manifest["sumstats"]
    return lw.run(
        s["t2d"],
        s["bmi_female"],
        s["bmi_male"],
        manifest["ld"],
        outdir=tmp_path_factory.mktemp("ldsc"),
    )


def test_hapmap3_filter_keeps_panel_snps(ldsc, truth):
    # Every panel SNP is HapMap3; only the 5 off-panel tutorial variants drop out.
    assert ldsc["n_variants"] == truth["n_panel_snps"] + 5
    assert ldsc["n_hapmap3"] == truth["n_panel_snps"]


def test_h2_recovers_truth(ldsc, truth):
    row = ldsc["h2"].iloc[0]
    h2, se = float(row["h2_obs"]), float(row["h2_se"])
    assert abs(h2 - truth["h2"]["T2D"]) < Z * se
    assert se < 0.01


def test_h2_intercept_recovers_inflation(ldsc, truth):
    row = ldsc["h2"].iloc[0]
    intercept, se = float(row["Intercept"]), float(row["Intercept_se"])
    assert abs(intercept - truth["intercept"]["T2D"]) < Z * se
    assert float(row["Mean_chi2"]) > float(row["Intercept"])


def test_rg_table_has_tutorial_columns(ldsc):
    rg = ldsc["rg"]
    for col in ["p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "gcov_int"]:
        assert col in rg.columns
    assert list(rg["p2"]) == ["BMI_female", "BMI_male"]
    assert set(rg["p1"]) == {"T2D"}


@pytest.mark.parametrize("trait", ["BMI_female", "BMI_male"])
def test_rg_recovers_truth(ldsc, truth, trait):
    row = ldsc["rg"].set_index("p2").loc[trait]
    assert abs(float(row["rg"]) - truth["rg"][trait]) < Z * float(row["se"])
    assert float(row["p"]) < 1e-6
    assert abs(float(row["h2_obs"]) - truth["h2"][trait]) < Z * float(row["h2_obs_se"])


def test_rg_no_sample_overlap_intercept(ldsc):
    # Independent synthetic cohorts: genetic-covariance intercept ~ 0.
    for _, row in ldsc["rg"].iterrows():
        assert abs(float(row["gcov_int"])) < Z * float(row["gcov_int_se"])


def test_rg_plot_written(ldsc):
    assert ldsc["rg_plot"].exists() and ldsc["rg_plot"].stat().st_size > 0


def test_partitioned_h2_enrichment(ldsc, truth):
    res = ldsc["partitioned"].set_index("Category")
    enh = res.loc["EnhancerL2_0"]
    enrichment, se = float(enh["Enrichment"]), float(enh["Enrichment_std_error"])
    assert enrichment > 1 + Z * se
    assert abs(enrichment - truth["enhancer_enrichment"]) < Z * se
    assert float(res.loc["baseL2_0", "Prop._SNPs"]) == pytest.approx(1.0)


def test_partitioned_summary_matches_total_h2(ldsc, truth):
    row = ldsc["partitioned_summary"].iloc[0]
    assert abs(float(row["h2_obs"]) - truth["h2"]["T2D"]) < Z * float(row["h2_se"])


def test_cell_type_analysis_ranks_causal_tissue_first(ldsc, truth):
    cts = ldsc["cts"].reset_index(drop=True)
    assert set(cts["Name"]) == {"Pancreas", "Liver", "Brain"}
    assert cts.loc[0, "Name"] == truth["causal_tissue"]
    assert cts.loc[0, "Coefficient_P_value"] < 1e-6


def test_result_tables_written(ldsc):
    outdir = ldsc["rg_plot"].parent
    for name in [
        "ldsc_h2.tsv",
        "ldsc_rg.tsv",
        "ldsc_partitioned_h2_summary.tsv",
        "ldsc_partitioned_h2.tsv",
        "ldsc_h2_cts.tsv",
    ]:
        assert (outdir / name).exists(), name


def test_h2_is_invariant_to_liftover(manifest, ldsc, tmp_path, monkeypatch):
    """Liftover + LDSC together: h2 on hg38 coordinates matches the hg19 estimate."""
    monkeypatch.chdir(tmp_path)  # LDSC writes scratch files to the cwd
    t2d = lw.load_t2d(manifest["sumstats"]["t2d"])
    t2d.liftover(from_build="19", to_build="38", verbose=False)
    assert t2d.build == "38"
    ref = manifest["ld"]["ldscores"]
    t2d.estimate_h2_by_ldsc(ref_ld_chr=ref, w_ld_chr=ref, verbose=False)
    h2_38 = float(t2d.ldsc_h2.iloc[0]["h2_obs"])
    h2_19 = float(ldsc["h2"].iloc[0]["h2_obs"])
    assert h2_38 == pytest.approx(h2_19, rel=0.02)
