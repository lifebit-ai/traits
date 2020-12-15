# lifebit-ai/traits: Usage

## Introduction

Runs heritability in your GWAS summary statistics as well as computing the genetic correlation between your trait of interest and a second trait with gwas summary statistics. 

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run lifebit-ai/traits -profile binary_h2
```

## 1 - Information about the method & how it is used

Vignette: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

**IMPORTANT NOTE**: it uses `BETA` column as `BETA`, outputs of SAIGE provide this column by default. 

## 2 - Parameters:
### 2.1 - Required parameters
- **--post_analysis** : String with `genetic_correlation_h2` or `heritability` for running genetic correlation analysis or heritability after GWAS.
- **--gwas_summary** : Path/URL to external gwas summary statistics to run genetic correlation analysis between cohort of interest and external GWAS summary statistics. The following column names and format (can also be comma-separated instead of whitespace-separated) are required to ensure it works:
```
snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
rs3131972	1	742584	A	G	1.092	0.0817	0.2819	0.694	0	0.16055
rs3131969	1	744045	A	G	1.087	0.0781	0.2855	0.939	0	0.133028
rs3131967	1	744197	T	C	1.093	0.0835	0.2859	0.869	0	.
rs1048488	1	750775	T	C	0.9158	0.0817	0.2817	0.694	0	0.836449
rs12562034	1	758311	A	G	0.9391	0.0807	0.4362	0.977	0	0.0925926
rs4040617	1	769185	A	G	0.9205	0.0777	0.2864	0.98	0	0.87156
rs28576697	1	860508	T	C	1.079	0.2305	0.7423	0.123	0	0.74537
rs1110052	1	863421	T	G	1.088	0.2209	0.702	0.137	0	0.752294
rs7523549	1	869180	T	C	1.823	0.8756	0.4929	0.13	0	0.0137615
``` 

### 2.2 - Optional parameters

- **--post_analysis** : String with `genetic_correlation_h2` or `heritability` for running genetic correlation analysis or heritability after GWAS.
- **--gwas_summary** : Path/URL to external gwas summary statistics to run genetic correlation analysis between cohort of interest and external GWAS summary statistics. The following column names and format (can also be comma-separated instead of whitespace-separated) are required to ensure it works:
- **--output_tag** : String with tag for files
- **--gwas_cat_study_id** : String with ID from GWAS catalogue study to be used as input for genetic correlation
- **--gwas_cat_study_size**  : Integer with size of study being used for genetic correlation
- **--gwas_catalogue_ftp** : Path to ftp locations of harmonized GWAS catalogue studies. Defaults to "https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/gel/prs/ftp_locations_harmonized.csv"
- **--hapmap3_snplist** : Path to snp list from hapmap3 to be used for analysis. Defaults to "https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/gel/gel-gwas/assets/w_hm3.snplist"
- **--ld_scores_tar_bz2** : Path to precomputed LD scores from Defaults to the European 1000 Genomes cohort at "https://lifebit-featured-datasets.s3-eu-west-1.amazonaws.com/projects/gel/gel-gwas/assets/eur_w_ld_chr.tar.bz2"
- **--outdir** : Path to output directory. Defaults to './results'