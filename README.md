# gwaslab-liftover-ldsc

Runnable walkthroughs and tests for two GWASLab docs pages:

- **Liftover** — <https://cloufield.github.io/gwaslab/harmonization_liftover/>
- **LD score regression in GWASLab** — <https://cloufield.github.io/gwaslab/ldsc_in_gwaslab/>

Each script uses the same calls and arguments as the tutorial. By default it
runs on **synthetic data with a known answer**, so every estimate can be checked
against the truth. Pass real file paths to run it on the tutorial's BBJ data.

## What is covered

| Docs section | Script / test |
| --- | --- |
| Load 5,000 T2D variants, `basic_check()` | `liftover_walkthrough.py`, `test_liftover.py` |
| `liftover(from_build="19", to_build="38")` (bundled UCSC chain) | `test_liftover.py` |
| `liftover(from_build="38", to_build="13", chain_path=...)` (T2T-CHM13) | `test_liftover.py` |
| `filter_hapmap3()` | `ldsc_walkthrough.py`, `test_ldsc.py` |
| `estimate_h2_by_ldsc` → `ldsc_h2` | `test_ldsc.py` |
| `estimate_rg_by_ldsc` (T2D vs BMI female / male) → `ldsc_rg`, `gl.plot_rg` | `test_ldsc.py` |
| `estimate_partitioned_h2_by_ldsc` → `ldsc_partitioned_h2_summary/results` | `test_ldsc.py` |
| `estimate_h2_cts_by_ldsc` → `ldsc_h2_cts` | `test_ldsc.py` |
| Liftover + LDSC together: h2 is unchanged after hg19 → hg38 | `test_ldsc.py` |

## Synthetic data

The tutorials use BioBank Japan T2D/BMI sumstats and EAS LD Score references.
These are large downloads, and JENGER, which hosts the BBJ files, could not be
reached when this repository was written. `scripts/synthetic.py` therefore builds
stand-ins:

- **Sumstats.** The files have the BBJ columns (`SNP CHR POS REF ALT Frq BETA SE P Dir N`)
  and sit on ~100k real HapMap3 SNPs from GWASLab's bundled hg19 list.
  - The T2D file uses `CHR:POS` IDs.
  - The BMI files use rsIDs, with the effect allele in `REF`, matching the tutorial's `ea="REF"`.
  - The T2D file also includes the five variants that the tutorial reports as
    unmapped by hg19 → hg38 (for example `1:1581586_T_C`).
- **LDSC references.** These use the LDSC on-disk format:
  - per-chromosome `.l2.ldscore.gz`, `.l2.M_5_50` and `.annot.gz` files;
  - a `baseline.` model with an `Enhancer` annotation;
  - `weights.EAS.hm3_noMHC.` weights;
  - a `Multi_tissue_gene_expr` `.ldcts` file listing Pancreas, Liver and Brain.
- **Simulation.** Z-scores are drawn from the LDSC model, `E[χ²] = N·Σ τ_c ℓ(j,c) + a`.
  The true values are T2D h2 = 0.10 with intercept 1.03, rg = 0.40 and 0.25, an
  enhancer enrichment of about 2.8, and Pancreas as the causal tissue. They are
  written to `manifest.json`.

## Setup

GWASLab 4.2.x pins `pysam==0.22.1` and `matplotlib<3.9`. These have prebuilt
wheels only up to **Python 3.12**, so use a 3.12 environment:

```bash
uv venv --python 3.12 .venv        # or: python3.12 -m venv .venv
source .venv/bin/activate
uv pip install -r requirements.txt # or: pip install -r requirements.txt
```

### Older CPUs (no AVX2)

If `import gwaslab` dies with `Illegal instruction`, the default `polars` build
needs AVX/AVX2. Install the compatibility runtime instead:

```bash
pip install "polars[rtcompat]"
```

## Run

```bash
scripts/download_data.sh                 # T2T-CHM13 chain (6 MB); enables the CHM13 test
pytest                                   # 28 tests, ~2 min

python scripts/liftover_walkthrough.py   # -> results/liftover/
python scripts/ldsc_walkthrough.py       # -> results/ldsc/
```

Without the chain file, the CHM13 test is skipped and the walkthrough stops at hg38.

### On the real tutorial data

```bash
scripts/download_data.sh --bbj --ldsc --baseline   # see docs/data.md for sizes/caveats

python scripts/liftover_walkthrough.py --sumstats data/sample_data/t2d_bbj.txt.gz

W=data/reference/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC.
python scripts/ldsc_walkthrough.py \
  --t2d data/sample_data/t2d_bbj.txt.gz \
  --bmi-female data/sample_data/bbj_bmi_female.txt.gz \
  --bmi-male data/sample_data/bbj_bmi_male.txt.gz \
  --ldscores "$W" --weights "$W" \
  --baseline <baselineLD_v2.2 prefix> --baseline-cts <baseline_v1.2 prefix> \
  --ldcts <Multi_tissue_gene_expr.EAS.ldcts> --steps h2,rg,partitioned,cts
```

## Results

See [`docs/results.md`](docs/results.md). In short, on synthetic data:

| Analysis | Estimate | Truth |
| --- | --- | --- |
| hg19 → hg38 → CHM13 (5,000 rows) | 5000 → 4991 → 4990 | the tutorial's 5 unmapped variants are dropped |
| T2D h2 (intercept) | 0.097 ± 0.004 (1.040 ± 0.016) | 0.10 (1.03) |
| rg T2D–BMI female / male | 0.40 ± 0.04 / 0.26 ± 0.04 | 0.40 / 0.25 |
| Enhancer enrichment | 3.1 ± 0.24 | 2.8 |
| Top cell type | Pancreas (p = 4e-24) | Pancreas |

## Layout

```
scripts/  synthetic.py              synthetic sumstats + LDSC reference files
          liftover_walkthrough.py   liftover tutorial (CLI + run())
          ldsc_walkthrough.py       LDSC tutorial (CLI + run())
          download_data.sh          T2T chain, EAS LDSC references, BBJ sumstats
tests/    test_synthetic.py test_liftover.py test_ldsc.py
docs/     data.md results.md img/
```
