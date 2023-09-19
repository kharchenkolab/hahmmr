# HaHMMR

<!-- badges: start -->
[![<kharchenkolab>](https://circleci.com/gh/kharchenkolab/hahmmr.svg?style=svg)](https://app.circleci.com/pipelines/github/kharchenkolab/hahmmr)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8342630.svg)](https://doi.org/10.5281/zenodo.8342630)
<!-- badges: end -->

<img src="hahmmr_logo.png" align="right" width="150">

**H**aplotype-**a**ware **H**idden **M**arkov **M**odel for **R**NA (HaHMMR) is a method for detecting CNVs from bulk RNA-seq data. Extending the haplotype-aware HMM in [Numbat](https://github.com/kharchenkolab/numbat) for single-cell RNA-seq, HaHMMR offers enhanced capabilities for detecting low-clonality CNVs from bulk data.

# Usage

## Preparing data
First, obtain expression counts and phased allele counts from the RNA-seq sample. The expression counts can be prepared using a transcript quantification tool such as [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html). The phased allele counts can be prepared using the [pileup_and_phase.R](https://kharchenkolab.github.io/numbat/articles/numbat.html#preparing-data) pipeline from Numbat. A [Docker](https://kharchenkolab.github.io/numbat/articles/numbat.html#docker) container is available for running this pipeline.

For example, within the Numbat Docker you can run `pileup_and_phase` in bulk RNA-seq mode like this:
```
Rscript /numbat/inst/bin/pileup_and_phase.R \
    --bulk \
    --label {sample} \
    --samples {sample} \
    --bams /mnt/mydata/{sample}.bam \
    --outdir /mnt/mydata/{sample} \
    --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
    --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --paneldir /data/1000G_hg38 \
    --ncores ncores
```

The **integer** expression counts (`count_mat`) should be a one-column matrix where rownames are genes and colname is the sample name. The phased allele counts (`df_allele`) should be a dataframe containing columns `snp_id`, `CHROM`, `POS`, `cM` (genetic distance in centimorgan), `REF`, `ALT`, `AD` (ALT allele count), `DP` (total allele count), `GT` (phased genotype), `gene`. 

## Running HaHMMR
Here is an example using the RNA-seq samples from a Meningioma [study](https://pubmed.ncbi.nlm.nih.gov/27548314/).

```
library(dplyr)
library(hahmmr)
allele_counts = data.table::fread('http://pklab.med.harvard.edu/teng/data/hmm_example/MN-5_TUMOR_allele_counts.tsv.gz')
gene_counts = readRDS(url('http://pklab.med.harvard.edu/teng/data/hmm_example/MN_gene_counts.rds'))
```

Sample MN-1037 has a diploid genome so we can use it to create a reference expression profile.

```
ref_internal = gene_counts[,'MN-1037_TUMOR',drop=F] %>% {./sum(.)}
head(ref_internal)
##          MN-1037_TUMOR
## 7SK       0.000000e+00
## A1BG      1.107976e-06
## A1BG-AS1  5.003764e-07
## A1CF      3.574117e-08
## A2ML1     3.931529e-07
## A4GALT    9.314150e-05
```

We can now analyze it using HaHMMR.

```
sample = 'MN-5_TUMOR'

bulk = get_bulk(
        count_mat = gene_counts[,sample,drop=F],
        df_allele = allele_counts,
        lambdas_ref = ref_internal,
        gtf = gtf_hg38
    ) %>% 
    analyze_joint()

bulk %>% plot_psbulk(min_depth = 15)
```
