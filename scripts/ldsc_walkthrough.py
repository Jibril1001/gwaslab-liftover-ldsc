"""Runnable walkthrough of the GWASLab "LD score regression" tutorial.

Docs: https://cloufield.github.io/gwaslab/ldsc_in_gwaslab/

Steps (same calls and arguments as the tutorial):

1. Load T2D sumstats (hg19) and run ``basic_check()``.
2. ``filter_hapmap3()`` (shown for inspection; the LDSC methods also extract
   HapMap3 SNPs internally).
3. ``estimate_h2_by_ldsc`` -> ``ldsc_h2``.
4. ``estimate_rg_by_ldsc`` against BMI (female, male) -> ``ldsc_rg`` and
   ``gl.plot_rg``.
5. ``estimate_partitioned_h2_by_ldsc`` with a baseline model ->
   ``ldsc_partitioned_h2_summary`` / ``ldsc_partitioned_h2_results``.
6. ``estimate_h2_cts_by_ldsc`` with a multi-tissue ``.ldcts`` file ->
   ``ldsc_h2_cts``.

By default it runs on synthetic data (``scripts/synthetic.py``). To use the
real BBJ files and EAS LD Score references, pass the paths explicitly (see
``--help`` and ``docs/data.md``).

Every LDSC run writes a scratch ``ldsc.delete``/``ldsc.part_delete`` file in
the working directory, so the script ``chdir``s into ``--outdir``.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from contextlib import contextmanager
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import gwaslab as gl  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import synthetic  # noqa: E402


def load_t2d(path: str | Path) -> gl.Sumstats:
    t2d = gl.Sumstats(
        str(path),
        snpid="SNP",
        chrom="CHR",
        pos="POS",
        ea="ALT",
        nea="REF",
        beta="BETA",
        se="SE",
        p="P",
        direction="Dir",
        build="19",
        n="N",
        verbose=False,
    )
    t2d.basic_check(verbose=False)
    return t2d


def load_bmi(path: str | Path, n: int) -> gl.Sumstats:
    # BBJ BMI files report the effect allele in the REF column (hence ea="REF").
    return gl.Sumstats(
        str(path),
        fmt="auto",
        ea="REF",
        nea="ALT",
        rsid="SNP",
        n=n,
        sep="\t",
        build="19",
        verbose=False,
    )


@contextmanager
def _cwd(path: Path):
    old = Path.cwd()
    path.mkdir(parents=True, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def run(
    t2d_path,
    bmi_female_path,
    bmi_male_path,
    ref: dict,
    n_bmi_female: int = 70000,
    n_bmi_male: int = 80000,
    outdir: str | Path = "results/ldsc",
    steps=("h2", "rg", "partitioned", "cts"),
    verbose: bool = False,
) -> dict:
    """Run the tutorial's LDSC analyses.

    ``ref`` holds LDSC-style prefixes: ``ldscores`` (h2/rg ref and weights),
    ``baseline``, ``weights`` (partitioned/cts), ``frq`` (optional) and
    ``ldcts``.
    """
    outdir = Path(outdir).resolve()
    t2d = load_t2d(Path(t2d_path).resolve())
    out: dict = {"t2d": t2d}

    hm3 = t2d.filter_hapmap3()
    out["n_variants"] = len(t2d.data)
    out["n_hapmap3"] = len(hm3.data)

    with _cwd(outdir):
        if "h2" in steps:
            t2d.estimate_h2_by_ldsc(ref_ld_chr=ref["ldscores"], w_ld_chr=ref["ldscores"], verbose=verbose)
            out["h2"] = t2d.ldsc_h2
            t2d.ldsc_h2.to_csv("ldsc_h2.tsv", sep="\t", index=False)

        if "rg" in steps:
            bmi_female = load_bmi(Path(bmi_female_path).resolve(), n_bmi_female)
            bmi_male = load_bmi(Path(bmi_male_path).resolve(), n_bmi_male)
            t2d.estimate_rg_by_ldsc(
                other_traits=[bmi_female, bmi_male],
                rg="T2D,BMI_female,BMI_male",
                ref_ld_chr=ref["ldscores"],
                w_ld_chr=ref["ldscores"],
                verbose=verbose,
            )
            out["rg"] = t2d.ldsc_rg
            t2d.ldsc_rg.to_csv("ldsc_rg.tsv", sep="\t", index=False)
            fig, *_ = gl.plot_rg(t2d.ldsc_rg, verbose=verbose)
            fig.savefig("ldsc_rg.png", dpi=150, bbox_inches="tight")
            out["rg_plot"] = outdir / "ldsc_rg.png"

        if "partitioned" in steps:
            kwargs = {}
            if ref.get("frq"):
                kwargs["frqfile_chr"] = ref["frq"]
            t2d.estimate_partitioned_h2_by_ldsc(
                ref_ld_chr=ref["baseline"],
                w_ld_chr=ref["weights"],
                overlap_annot=True,
                print_coefficients=True,
                print_delete_vals=True,
                verbose=verbose,
                **kwargs,
            )
            out["partitioned_summary"] = t2d.ldsc_partitioned_h2_summary
            out["partitioned"] = t2d.ldsc_partitioned_h2_results
            t2d.ldsc_partitioned_h2_summary.to_csv("ldsc_partitioned_h2_summary.tsv", sep="\t", index=False)
            t2d.ldsc_partitioned_h2_results.to_csv("ldsc_partitioned_h2.tsv", sep="\t", index=False)

        if "cts" in steps:
            t2d.estimate_h2_cts_by_ldsc(
                ref_ld_chr=ref.get("baseline_cts", ref["baseline"]),
                ref_ld_chr_cts=ref["ldcts"],
                w_ld_chr=ref["weights"],
                verbose=verbose,
            )
            out["cts"] = t2d.ldsc_h2_cts
            t2d.ldsc_h2_cts.to_csv("ldsc_h2_cts.tsv", sep="\t", index=False)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="GWASLab LDSC tutorial walkthrough")
    ap.add_argument("--synthetic-dir", default="data/synthetic")
    ap.add_argument("--t2d", help="T2D sumstats (default: synthetic)")
    ap.add_argument("--bmi-female", help="BMI female sumstats")
    ap.add_argument("--bmi-male", help="BMI male sumstats")
    ap.add_argument("--ldscores", help="--ref-ld-chr/--w-ld-chr prefix for h2 and rg, e.g. eas_ldscores/")
    ap.add_argument("--baseline", help="baseline-LD prefix, e.g. .../baselineLD2_2/baselineLD.")
    ap.add_argument("--baseline-cts", help="baseline prefix for cts (tutorial: baseline1_2/baseline.)")
    ap.add_argument("--weights", help="weights prefix, e.g. .../weights.EAS.hm3_noMHC.")
    ap.add_argument("--frq", help="frequency prefix, e.g. .../1000G.EAS.QC.")
    ap.add_argument("--ldcts", help=".ldcts file for cell-type-specific analysis")
    ap.add_argument("--steps", default="h2,rg,partitioned,cts")
    ap.add_argument("--outdir", default="results/ldsc")
    ap.add_argument("--verbose", action="store_true", help="print GWASLab logs")
    args = ap.parse_args()

    if args.t2d:
        files = {"t2d": args.t2d, "bmi_female": args.bmi_female, "bmi_male": args.bmi_male}
        ref = {
            "ldscores": args.ldscores,
            "baseline": args.baseline,
            "baseline_cts": args.baseline_cts or args.baseline,
            "weights": args.weights,
            "frq": args.frq,
            "ldcts": args.ldcts,
        }
        truth = None
    else:
        manifest_path = Path(args.synthetic_dir) / "manifest.json"
        if not manifest_path.exists():
            synthetic.generate(args.synthetic_dir)
        manifest = json.loads(manifest_path.read_text())
        files, ref, truth = manifest["sumstats"], manifest["ld"], manifest["truth"]

    out = run(
        files["t2d"],
        files["bmi_female"],
        files["bmi_male"],
        ref,
        outdir=args.outdir,
        steps=tuple(args.steps.split(",")),
        verbose=args.verbose,
    )

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print(f"Variants: {out['n_variants']:,}  HapMap3: {out['n_hapmap3']:,}")
    for key in ["h2", "rg", "partitioned_summary", "partitioned", "cts"]:
        if key in out:
            print(f"\n== {key} ==")
            print(out[key].to_string(index=False))
    if truth:
        print("\n== simulated truth ==")
        keys = ["h2", "rg", "enhancer_enrichment", "causal_tissue"]
        print(json.dumps({k: truth[k] for k in keys}, indent=2))
    print(f"\nResults written to {args.outdir}/")


if __name__ == "__main__":
    main()
