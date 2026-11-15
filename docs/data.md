# Data sources

## Used by the tutorials

| File | Source | Size | Notes |
| --- | --- | --- | --- |
| `t2d_bbj.txt.gz` | `http://jenger.riken.jp/14/` | large | BBJ T2D. **JENGER timed out when this repo was written (2026-11-15).** |
| `bbj_bmi_male.txt.gz` | `http://jenger.riken.jp/2analysisresult_qtl_download/` | large | BBJ BMI male (tutorial: `n=80000`) |
| `bbj_bmi_female.txt.gz` | `http://jenger.riken.jp/4analysisresult_qtl_download/` | large | BBJ BMI female (tutorial: `n=70000`) |
| `hg19ToHg38.over.chain.gz` | bundled with GWASLab (`gwaslab/data/chains/`) | — | used by `liftover(from_build="19", to_build="38")` |
| `grch38-chm13v2.chain` | [human-pangenomics S3](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain) | 6.0 MB | hg38 → T2T-CHM13 v2 |
| EAS weights `weights.EAS.hm3_noMHC.` | [Zenodo 10515792](https://zenodo.org/records/10515792) `1000G_Phase3_EAS_weights_hm3_no_MHC.tgz` | 11.5 MB | verified to extract to `1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC.*` |
| EAS baseline v1.2 | Zenodo 10515792 `1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz` | 280 MB | tutorial's `baseline1_2/baseline.` (cell-type analysis) |
| EAS baselineLD v2.2 | Zenodo 10515792 `1000G_Phase3_EAS_baselineLD_v2.2_ldscores.tgz` | 569 MB | tutorial's `baselineLD2_2/baselineLD.` (partitioned h2) |

The old Broad `alkesgroup/LDSCORE` Google Cloud URLs return 404. The S-LDSC
reference files now live on Zenodo.

Two of the tutorial's references are not in that Zenodo record:

- The plain `eas_ldscores/` files used for h2 and rg. The EAS HapMap3 weights
  above are a standard substitute, used as both `ref_ld_chr` and `w_ld_chr`.
- The `Multi_tissue_gene_expr.EAS.ldcts` cell-type files. These come with the
  LDSC cell-type-specific analysis (Finucane et al. 2018) and must be obtained
  separately.

## Synthetic stand-ins

`scripts/synthetic.py` (seed 2026) writes the files below to `data/synthetic/`.
The tests write them to a temporary directory instead.

```
manifest.json                          paths + simulated truth
sample_data/t2d_bbj.synthetic.txt.gz   100,559 rows (100,554 HapMap3 + 5 tutorial-unmapped)
sample_data/bbj_bmi_{female,male}.synthetic.txt.gz   ~97% of the panel each
sample_data/panel_hg19_hg38.tsv.gz     rsID, hg19 and dbSNP151 hg38 positions
ldscores/eas_ldscores/{1..22}.l2.*     single L2 (h2, rg)
ldscores/eas_baseline/baseline.{1..22}.{l2.ldscore.gz,l2.M_5_50,annot.gz}
ldscores/eas_weights/weights.EAS.hm3_noMHC.{1..22}.l2.*
ldscores/Multi_tissue_gene_expr/{Pancreas,Liver,Brain,control}.{1..22}.l2.*
ldscores/Multi_tissue_gene_expr/Multi_tissue_gene_expr.EAS.ldcts
```

The panel is every 12th autosomal, biallelic SNP from GWASLab's HapMap3 list.
Only SNPs that are unique and present in both the hg19 (dbSNP150) and hg38
(dbSNP151) lists are kept, so lifted positions can be checked against dbSNP.
