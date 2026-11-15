"""Synthetic stand-ins for the data used by the GWASLab liftover and LDSC tutorials.

The tutorials use BioBank Japan summary statistics (``t2d_bbj.txt.gz``,
``bbj_bmi_female.txt.gz``, ``bbj_bmi_male.txt.gz``) and pre-computed EAS LD
Score references. Those are large downloads, so this module generates small
files with the same layout:

* Summary statistics with the BBJ columns ``SNP CHR POS REF ALT Frq BETA SE P
  Dir N``, placed on real HapMap3 SNPs (hg19) taken from GWASLab's bundled
  HapMap3 list, so ``liftover`` and ``filter_hapmap3`` behave as on real data.
* LD Score reference files in the LDSC on-disk format (per-chromosome
  ``.l2.ldscore.gz``, ``.l2.M_5_50`` and ``.annot.gz`` files plus an
  ``.ldcts`` list) for the h2/rg, partitioned h2 and cell-type analyses.

Z-scores are drawn from the LDSC model, ``E[chi2_j] = N * sum_c tau_c l(j,c) +
intercept``, so the true heritability, genetic correlation, enrichment and
causal tissue are known and recorded in ``truth.json``.
"""

from __future__ import annotations

import argparse
import json
from importlib import resources
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SEED = 2026
CHROMS = list(range(1, 23))

# Variants the tutorial reports as unmapped by hg19 -> hg38 liftover.
TUTORIAL_UNMAPPED = [
    ("1", 1581586, "T", "C"),
    ("1", 1584117, "G", "C"),
    ("1", 1584218, "T", "C"),
    ("1", 1585257, "A", "G"),
    ("1", 1585313, "CACAA", "C"),
]

TISSUES = ["Pancreas", "Liver", "Brain"]
CAUSAL_TISSUE = "Pancreas"

# Trait design (all chosen so mean chi2 is close to the tutorial's ~1.5).
N_T2D = 25_000
N_BMI_FEMALE = 70_000  # matches n=70000 in the tutorial
N_BMI_MALE = 80_000  # matches n=80000 in the tutorial
H2_BASE_T2D = 0.06
H2_ENHANCER_T2D = 0.02
H2_TISSUE_T2D = 0.02
H2_BMI = 0.015
INTERCEPT_T2D = 1.03
RG_BMI_FEMALE = 0.4
RG_BMI_MALE = 0.25


def hapmap3_path(build: str = "19") -> Path:
    name = "hapmap3_db150_hg19.snplist.gz" if build == "19" else "hapmap3_db151_hg38.snplist.gz"
    return Path(str(resources.files("gwaslab") / "data" / "hapmap3_SNPs" / name))


def load_panel(step: int = 12) -> pd.DataFrame:
    """Every ``step``-th autosomal biallelic HapMap3 SNP present in both hg19 and hg38 lists."""
    hg19 = pd.read_csv(hapmap3_path("19"), sep="\t", dtype={"#CHROM": str})
    hg38 = pd.read_csv(hapmap3_path("38"), sep="\t", dtype={"#CHROM": str})
    hg19 = hg19[hg19["#CHROM"].isin([str(c) for c in CHROMS])]
    hg19 = hg19[hg19["A1"].str.len().eq(1) & hg19["A2"].str.len().eq(1)]
    hg19 = hg19.drop_duplicates("rsid", keep=False).drop_duplicates(["#CHROM", "POS"], keep=False)
    hg38 = hg38.drop_duplicates("rsid", keep=False)[["rsid", "#CHROM", "POS"]]
    hg38 = hg38.rename(columns={"#CHROM": "CHR38", "POS": "POS38"})
    panel = hg19.merge(hg38, on="rsid", how="inner")
    panel = panel[panel["#CHROM"] == panel["CHR38"]]
    panel = panel.rename(columns={"rsid": "SNP", "#CHROM": "CHR", "POS": "BP"})
    panel["CHR"] = panel["CHR"].astype(int)
    panel = panel.sort_values(["CHR", "BP"]).iloc[::step].reset_index(drop=True)
    return panel[["CHR", "SNP", "BP", "A1", "A2", "POS38"]]


def _ld_component(rng, base: np.ndarray, annot: np.ndarray, share: float) -> np.ndarray:
    """LD Score to an annotation: proportional to the SNP's total LD, larger for annotated SNPs."""
    mix = share * (1 + 4 * annot) / (1 + 4 * annot.mean())
    return base * mix * rng.uniform(0.5, 1.5, size=base.size)


def build_ldscores(panel: pd.DataFrame, rng) -> dict:
    n = len(panel)
    base = rng.gamma(shape=2.0, scale=10.0, size=n)  # mean LD Score ~20
    enh = rng.binomial(1, 0.10, size=n).astype(float)
    tissue_annot = {t: rng.binomial(1, 0.10, size=n).astype(float) for t in TISSUES}
    control = rng.binomial(1, 0.30, size=n).astype(float)
    return {
        "base": base,
        "enhancer_annot": enh,
        "enhancer": _ld_component(rng, base, enh, 0.10),
        "tissue_annot": tissue_annot,
        "tissue": {t: _ld_component(rng, base, a, 0.10) for t, a in tissue_annot.items()},
        "control_annot": control,
        "control": _ld_component(rng, base, control, 0.30),
    }


def _bbj_table(panel, z, n, rng, swap_alleles=False, snp_as_rsid=False) -> pd.DataFrame:
    frq = rng.uniform(0.05, 0.95, size=len(panel))
    se = 1.0 / np.sqrt(2 * frq * (1 - frq) * n)
    beta = z * se
    p = 2 * stats.norm.sf(np.abs(z))
    signs = np.where(rng.random((len(panel), 3)) < 0.5, "+", "-")
    signs[np.arange(len(panel)), 0] = np.where(z > 0, "+", "-")
    ref, alt = panel["A2"].to_numpy(), panel["A1"].to_numpy()
    if swap_alleles:  # BBJ BMI files label the effect allele as REF (tutorial uses ea="REF")
        ref, alt = alt, ref
    snp = (
        panel["SNP"].to_numpy()
        if snp_as_rsid
        else (panel["CHR"].astype(str) + ":" + panel["BP"].astype(str)).to_numpy()
    )
    return pd.DataFrame(
        {
            "SNP": snp,
            "CHR": panel["CHR"].to_numpy(),
            "POS": panel["BP"].to_numpy(),
            "REF": ref,
            "ALT": alt,
            "Frq": frq.round(4),
            "BETA": beta,
            "SE": se,
            "P": p,
            "Dir": ["".join(r) for r in signs],
            "N": n,
        }
    )


def simulate(panel: pd.DataFrame, ld: dict, rng) -> tuple[dict, dict]:
    m = len(panel)
    m_enh = ld["enhancer_annot"].sum()
    m_tis = ld["tissue_annot"][CAUSAL_TISSUE].sum()
    tau_base = H2_BASE_T2D / m
    tau_enh = H2_ENHANCER_T2D / m_enh
    tau_tis = H2_TISSUE_T2D / m_tis
    # Per-SNP expected genetic chi2 contribution / N for T2D and BMI
    g_t2d = tau_base * ld["base"] + tau_enh * ld["enhancer"] + tau_tis * ld["tissue"][CAUSAL_TISSUE]
    g_bmi = (H2_BMI / m) * ld["base"]

    tables = {}
    u_t2d = rng.standard_normal(m)
    z_t2d = np.sqrt(N_T2D * g_t2d) * u_t2d + np.sqrt(INTERCEPT_T2D) * rng.standard_normal(m)
    t2d = _bbj_table(panel, z_t2d, N_T2D, rng)
    extra = pd.DataFrame(
        {
            "SNP": [f"{c}:{p}" for c, p, _, _ in TUTORIAL_UNMAPPED],
            "CHR": [int(c) for c, *_ in TUTORIAL_UNMAPPED],
            "POS": [p for _, p, _, _ in TUTORIAL_UNMAPPED],
            "REF": [r for *_, r, _ in TUTORIAL_UNMAPPED],
            "ALT": [a for *_, a in TUTORIAL_UNMAPPED],
            "Frq": 0.2,
            "BETA": 0.01,
            "SE": 0.02,
            "P": 0.62,
            "Dir": "+-+",
            "N": N_T2D,
        }
    )
    tables["t2d"] = pd.concat([t2d, extra]).sort_values(["CHR", "POS"]).reset_index(drop=True)

    for key, n, rg in [("bmi_female", N_BMI_FEMALE, RG_BMI_FEMALE), ("bmi_male", N_BMI_MALE, RG_BMI_MALE)]:
        u = rg * u_t2d + np.sqrt(1 - rg**2) * rng.standard_normal(m)
        z = np.sqrt(n * g_bmi) * u + rng.standard_normal(m)
        df = _bbj_table(panel, z, n, rng, swap_alleles=True, snp_as_rsid=True)
        # Each study passes QC on a slightly different SNP set, as real files do.
        # (This also avoids a GWASLab 4.2.3 bug in LDSC's smart_merge fast path,
        # which is only taken when two traits have identical SNP lists; see docs.)
        tables[key] = df[rng.random(m) > 0.03].reset_index(drop=True)

    h2_t2d = H2_BASE_T2D + H2_ENHANCER_T2D + H2_TISSUE_T2D
    # Enhancer SNPs also carry base and tissue heritability.
    h2_in_enh = (
        tau_base * m_enh
        + H2_ENHANCER_T2D
        + tau_tis * (ld["enhancer_annot"] * ld["tissue_annot"][CAUSAL_TISSUE]).sum()
    )
    truth = {
        "n_panel_snps": int(m),
        "h2": {"T2D": h2_t2d, "BMI_female": H2_BMI, "BMI_male": H2_BMI},
        "intercept": {"T2D": INTERCEPT_T2D, "BMI_female": 1.0, "BMI_male": 1.0},
        "rg": {"BMI_female": RG_BMI_FEMALE, "BMI_male": RG_BMI_MALE},
        "enhancer_enrichment": float((h2_in_enh / h2_t2d) / (m_enh / m)),
        "causal_tissue": CAUSAL_TISSUE,
        "n": {"T2D": N_T2D, "BMI_female": N_BMI_FEMALE, "BMI_male": N_BMI_MALE},
        "unmapped_hg19_to_hg38": [f"{c}:{p}_{r}_{a}" for c, p, r, a in TUTORIAL_UNMAPPED],
    }
    return tables, truth


def _write_chr_files(prefix: str, panel: pd.DataFrame, ldcols: dict, annot: dict | None = None):
    """Write ``{prefix}{chr}.l2.ldscore.gz``, ``.l2.M_5_50`` and optionally ``.annot.gz``.

    ``prefix`` is an LDSC ``--ref-ld-chr`` style prefix: a directory ending in
    ``/`` or a directory plus a file stem such as ``.../baseline.``.
    """
    Path(f"{prefix}1").parent.mkdir(parents=True, exist_ok=True)
    for c in CHROMS:
        idx = (panel["CHR"] == c).to_numpy()
        sub = panel.loc[idx, ["CHR", "SNP", "BP"]].reset_index(drop=True)
        ldf = sub.copy()
        for name, values in ldcols.items():
            ldf[name] = np.round(values[idx], 3)
        ldf.to_csv(f"{prefix}{c}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        if annot is not None:
            m = [str(int(a[idx].sum())) for a in annot.values()]
            af = sub[["CHR", "BP", "SNP"]].copy()
            af["CM"] = 0
            for name, values in annot.items():
                af[name] = values[idx].astype(int)
            af.to_csv(f"{prefix}{c}.annot.gz", sep="\t", index=False, compression="gzip")
        else:
            m = [str(int(idx.sum()))] * len(ldcols)
        Path(f"{prefix}{c}.l2.M_5_50").write_text("\t".join(m) + "\n")
        Path(f"{prefix}{c}.l2.M").write_text("\t".join(m) + "\n")


def write_ldscores(outdir: Path, panel: pd.DataFrame, ld: dict) -> dict:
    """Write reference LD Score files mirroring the tutorial's EAS directories."""
    ones = np.ones(len(panel))
    cts_dir = outdir / "Multi_tissue_gene_expr"
    paths = {
        "ldscores": f"{outdir / 'eas_ldscores'}/",
        "baseline": f"{outdir / 'eas_baseline' / 'baseline'}.",
        "weights": f"{outdir / 'eas_weights' / 'weights.EAS.hm3_noMHC'}.",
        "ldcts": str(cts_dir / "Multi_tissue_gene_expr.EAS.ldcts"),
    }
    _write_chr_files(paths["ldscores"], panel, {"L2": ld["base"]})
    _write_chr_files(
        paths["baseline"],
        panel,
        {"baseL2": ld["base"], "EnhancerL2": ld["enhancer"]},
        annot={"base": ones, "Enhancer": ld["enhancer_annot"]},
    )
    _write_chr_files(paths["weights"], panel, {"L2": ld["base"]})
    _write_chr_files(f"{cts_dir / 'control'}.", panel, {"L2": ld["control"]})
    lines = []
    for t in TISSUES:
        _write_chr_files(f"{cts_dir / t}.", panel, {"L2": ld["tissue"][t]})
        # LDSC splits .ldcts lines on whitespace, so paths must not contain spaces.
        lines.append(f"{t}\t{cts_dir / t}.,{cts_dir / 'control'}.")
    Path(paths["ldcts"]).write_text("\n".join(lines) + "\n")
    return paths


def generate(outdir: str | Path, step: int = 12, seed: int = SEED) -> dict:
    """Generate all synthetic inputs under ``outdir`` and return a manifest of paths."""
    outdir = Path(outdir).resolve()  # absolute: LDSC runs chdir and .ldcts embeds paths
    rng = np.random.default_rng(seed)
    panel = load_panel(step)
    ld = build_ldscores(panel, rng)
    tables, truth = simulate(panel, ld, rng)

    sample = outdir / "sample_data"
    sample.mkdir(parents=True, exist_ok=True)
    files = {
        "t2d": sample / "t2d_bbj.synthetic.txt.gz",
        "bmi_female": sample / "bbj_bmi_female.synthetic.txt.gz",
        "bmi_male": sample / "bbj_bmi_male.synthetic.txt.gz",
    }
    for key, path in files.items():
        tables[key].to_csv(path, sep="\t", index=False, compression="gzip")
    panel[["SNP", "CHR", "BP", "POS38"]].to_csv(
        sample / "panel_hg19_hg38.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    manifest = {
        "sumstats": {k: str(v) for k, v in files.items()},
        "panel": str(sample / "panel_hg19_hg38.tsv.gz"),
        "ld": write_ldscores(outdir / "ldscores", panel, ld),
        "truth": truth,
        "seed": seed,
    }
    with open(outdir / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    return manifest


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default="data/synthetic")
    ap.add_argument("--step", type=int, default=12, help="keep every N-th HapMap3 SNP")
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()
    manifest = generate(args.outdir, args.step, args.seed)
    print(json.dumps(manifest["truth"], indent=2))
    print(f"Wrote synthetic data to {args.outdir}/ (manifest.json)")


if __name__ == "__main__":
    main()
