"""Runnable walkthrough of the GWASLab "Liftover" tutorial.

Docs: https://cloufield.github.io/gwaslab/harmonization_liftover/

Steps (same calls as the tutorial):

1. Load the first ``nrows`` variants of a BBJ-style T2D file (hg19).
2. ``basic_check()``.
3. ``liftover(from_build="19", to_build="38")`` with GWASLab's bundled
   UCSC ``hg19ToHg38.over.chain.gz``.
4. ``liftover(from_build="38", to_build="13", chain_path=...)`` with the
   T2T-CHM13 chain ``grch38-chm13v2.chain`` (optional; see
   ``scripts/download_data.sh``).

Usage::

    python scripts/liftover_walkthrough.py                       # synthetic data
    python scripts/liftover_walkthrough.py --sumstats t2d_bbj.txt.gz
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import gwaslab as gl
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import synthetic  # noqa: E402

DEFAULT_CHAIN = Path("data/reference/grch38-chm13v2.chain")


def load_t2d(path: str | Path, nrows: int | None = 5000) -> gl.Sumstats:
    """Load T2D sumstats exactly as the liftover tutorial does."""
    sumstats = gl.Sumstats(
        str(path),
        snpid="SNP",
        chrom="CHR",
        pos="POS",
        ea="ALT",
        nea="REF",
        neaf="Frq",
        beta="BETA",
        se="SE",
        p="P",
        nrows=nrows,
        verbose=False,
    )
    sumstats.basic_check(verbose=False)
    return sumstats


def _snapshot(sumstats: gl.Sumstats, label: str) -> pd.DataFrame:
    df = sumstats.data[["SNPID", "CHR", "POS", "STATUS"]].copy()
    df.columns = ["SNPID", f"CHR_{label}", f"POS_{label}", f"STATUS_{label}"]
    return df


def run(
    sumstats_path: str | Path,
    chain_t2t: str | Path | None = None,
    nrows: int | None = 5000,
    outdir: str | Path | None = None,
    verbose: bool = False,
) -> dict:
    """Run the tutorial's two liftovers and return per-step counts and coordinates."""
    sumstats = load_t2d(sumstats_path, nrows)
    steps = [{"step": "loaded", "build": sumstats.build, "variants": len(sumstats.data), "unmapped": 0}]
    coords = _snapshot(sumstats, "hg19")

    sumstats.liftover(from_build="19", to_build="38", verbose=verbose)
    steps.append({"step": "hg19 -> hg38", "build": sumstats.build, "variants": len(sumstats.data)})
    coords = coords.merge(_snapshot(sumstats, "hg38"), on="SNPID", how="left")

    if chain_t2t is not None and Path(chain_t2t).exists():
        sumstats.liftover(from_build="38", to_build="13", chain_path=str(chain_t2t), verbose=verbose)
        steps.append({"step": "hg38 -> T2T-CHM13", "build": sumstats.build, "variants": len(sumstats.data)})
        coords = coords.merge(_snapshot(sumstats, "chm13"), on="SNPID", how="left")

    for i in range(1, len(steps)):
        steps[i]["unmapped"] = steps[i - 1]["variants"] - steps[i]["variants"]

    result = {"steps": steps, "coords": coords, "sumstats": sumstats, "outdir": None}
    if outdir is not None:
        outdir = Path(outdir)
        result["outdir"] = outdir
        outdir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(steps).to_csv(outdir / "liftover_steps.tsv", sep="\t", index=False)
        coords.to_csv(outdir / "liftover_coordinates.tsv.gz", sep="\t", index=False)
        sumstats.data.to_csv(outdir / f"t2d_lifted_build{sumstats.build}.tsv.gz", sep="\t", index=False)
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description="GWASLab liftover tutorial walkthrough")
    ap.add_argument("--sumstats", help="BBJ-style T2D file (default: generate synthetic data)")
    ap.add_argument("--synthetic-dir", default="data/synthetic")
    ap.add_argument("--chain-t2t", default=str(DEFAULT_CHAIN), help="hg38 -> CHM13 chain file")
    ap.add_argument("--nrows", type=int, default=5000, help="rows to load (tutorial: 5000)")
    ap.add_argument("--outdir", default="results/liftover")
    ap.add_argument("--verbose", action="store_true", help="print GWASLab logs")
    args = ap.parse_args()

    path = args.sumstats
    if path is None:
        manifest = Path(args.synthetic_dir) / "manifest.json"
        if not manifest.exists():
            synthetic.generate(args.synthetic_dir)
        path = json.loads(manifest.read_text())["sumstats"]["t2d"]
    if not Path(args.chain_t2t).exists():
        print(
            f"T2T chain not found at {args.chain_t2t}; skipping hg38 -> CHM13 (run scripts/download_data.sh)."
        )

    result = run(path, args.chain_t2t, args.nrows, args.outdir, args.verbose)
    pd.set_option("display.width", 160)
    print(pd.DataFrame(result["steps"]).to_string(index=False))
    print()
    print(result["coords"].head(5).to_string(index=False))
    print(f"\nResults written to {args.outdir}/")


if __name__ == "__main__":
    main()
