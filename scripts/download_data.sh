#!/usr/bin/env bash
# Download the external files used by the GWASLab liftover and LDSC tutorials.
#
#   scripts/download_data.sh            # T2T-CHM13 chain only (6 MB, used by the tests)
#   scripts/download_data.sh --ldsc     # + EAS LDSC weights (12 MB)
#   scripts/download_data.sh --baseline # + EAS baseline / baselineLD v2.2 (~850 MB)
#   scripts/download_data.sh --bbj      # + BBJ T2D / BMI sumstats from JENGER
#
# Files go under data/reference/ and data/sample_data/ (both git-ignored).
set -euo pipefail

cd "$(dirname "$0")/.."
REF=data/reference
SAMPLE=data/sample_data
ZENODO=https://zenodo.org/records/10515792/files   # S-LDSC reference files
mkdir -p "$REF" "$SAMPLE"

fetch() {  # fetch URL DEST
  if [[ -s "$2" ]]; then echo "exists: $2"; return; fi
  echo "downloading $1"
  curl -fsSL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"
}

# hg38 -> T2T-CHM13 v2 chain (liftover tutorial, second liftover)
fetch https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain \
  "$REF/grch38-chm13v2.chain"

for arg in "$@"; do
  case "$arg" in
    --ldsc)
      fetch "$ZENODO/1000G_Phase3_EAS_weights_hm3_no_MHC.tgz?download=1" "$REF/eas_weights.tgz"
      tar -xzf "$REF/eas_weights.tgz" -C "$REF"
      ;;
    --baseline)
      fetch "$ZENODO/1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz?download=1" "$REF/eas_baseline_v1.2.tgz"
      fetch "$ZENODO/1000G_Phase3_EAS_baselineLD_v2.2_ldscores.tgz?download=1" "$REF/eas_baselineLD_v2.2.tgz"
      tar -xzf "$REF/eas_baseline_v1.2.tgz" -C "$REF"
      tar -xzf "$REF/eas_baselineLD_v2.2.tgz" -C "$REF"
      ;;
    --bbj)
      # Same sources as GWASLab's examples/0_sample_data/download_sample_data.sh
      fetch http://jenger.riken.jp/14/ "$SAMPLE/t2d_bbj.txt.gz"
      fetch http://jenger.riken.jp/2analysisresult_qtl_download/ "$SAMPLE/bbj_bmi_male.txt.gz"
      fetch http://jenger.riken.jp/4analysisresult_qtl_download/ "$SAMPLE/bbj_bmi_female.txt.gz"
      ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done
ls -la "$REF"
