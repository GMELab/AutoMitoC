# <p align="center"> ![image](https://github.com/GMELab/AutoMitoC/blob/main/AutoMitoC_LOGO_2021_11_26_SMALL.png) </p>
# <p align="center"> Automatic Mitochondrial DNA Copy Number (AutoMitoC) Pipeline to infer mtDNA-CN from genotyping arrays </p>

## Table of Contents

 - [Background](#Background)
 - [Key Features](#key-features)
 - [Method Overview](#method-overview)
 - [Pre-Processing](#pre-processing)
 - [Necessary Inputs](#necessary-inputs)
 - [Quick Start](#quick-start)
 - [Areas of Future Development](#areas-of-future-development)
 - [Contact information](#contact-information)
 - [Citation](#citation)
 - [Data Availability](#data-availability)


## Background

Mitochondria have their own unique DNA due to their origins as ancient bacteria. Due to this property, the ratio of mitochondrial to nuclear DNA signals, dubbed as the mitochondrial DNA copy number (mtDNA-CN) has emerged as a popular disease biomarker of mitochondrial etiology. One way to quantify mtDNA-CN is through genotyping arrays which enables not only the estimation of mtDNA-CN but also the capability to simultaneously interrogate its common genetic determinants. Previous array-based methods involve steps that are subjective, that require manual curation or extra lab testing, and that have not been extensively validated in different ethnicities. Accordingly, this served as the motivation to develop the AutoMitoC pipeline in order to enable fast, automated, and robust estimation of mtDNA-CN from genotyping arrays in all populations. 

## Key Features

AutoMitoC builds off of the [MitoPipeline (Lane _et al_, 2014)](http://genvisis.org/MitoPipeline/) with three key distinguishing features: 

* Nuclear signal is approximated by globally rare variants in place of common variants 
* Cross-hybridizing probes are identified by empirically assessing evidence for cross-hybridization
* The primary estimate of mitochondrial (MT) signal is ascertained using principal component analysis (PCA)
 
For further details on how AutoMitoC was developed or the rationale behind these features, please see the [main manuscript](https://www.medrxiv.org/content/10.1101/2021.04.08.21255031v1).

## Method Overview

![AutoMitoC_CHART_2021_11_26_FINAL](https://user-images.githubusercontent.com/30928727/143615288-fd36ea55-23b5-4ead-85dc-3e9c5dc2c1cd.jpeg)

This figure was created using [**lucid charts**](https://www.lucidchart.com/pages/landing?utm_source=google&utm_medium=cpc&utm_campaign=_chart_en_tier1_mixed_search_brand_exact_&km_CPC_CampaignId=1490375427&km_CPC_AdGroupID=55688909257&km_CPC_Keyword=lucid%20charts&km_CPC_MatchType=e&km_CPC_ExtensionID=&km_CPC_Network=g&km_CPC_AdPosition=&km_CPC_Creative=442433236001&km_CPC_TargetID=aud-921551091782:kwd-64262996435&km_CPC_Country=9000910&km_CPC_Device=c&km_CPC_placement=&km_CPC_target=&mkwid=stwZuDiXw_pcrid_442433236001_pkw_lucid%20charts_pmt_e_pdv_c_slid__pgrid_55688909257_ptaid_aud-921551091782:kwd-64262996435_&gclid=CjwKCAiAqIKNBhAIEiwAu_ZLDgo5PDn3MCLTVBtuG_3E53Epo4dl48Wgh9LAfXd2qfiKrLjN31UBQhoCHEkQAvD_BwE)

## Pre-processing

1. We recommend applying some basic PLINK QC for autosomal variants:
```sh
plink --bfile your_study \
--max-maf 0.01 \
--geno 0.05 \
--out your_study_qc
```

*Note that we have seen comparable performance at even lower max-maf thresholds (e.g. --max-maf 0.001 with ~ 10000 autosomal probes) but this is also dependent on your sample size. 

2. Follow step 1 of the PennCNV pipeline to generate [signal intensity files containing L2R and BAF values](http://penncnv.openbioinformatics.org/en/latest/user-guide/affy/#step-1-generate-the-signal-intensity-data-based-on-raw-cel-files).

3. Correct signal intensities (MT and Auto L2R) for GC-waves using the PennCNV implementation of the [Diskin _et al_ (2008) adjustment](https://github.com/WGLab/PennCNV/blob/master/genomic_wave.pl). 

4. Remove "wavy" samples with a post GC-corrected autosomal L2R standard deviation of > 0.35. 

## Necessary Inputs

The core AutoMitoC Rscript assumes all [pre-processing](#pre-processing) steps have been completed and requires the following inputs to be specified:

1. Path to the **Autosomal L2R matrix** ("EXAMPLE.AUTO_MATRIX.csv.gz"; comma separated; rows as samples; columns as autosomal probe L2R; no header or sample IDs)
2. Path to the **MT L2R matrix** ("EXAMPLE.MITO_MATRIX.csv.gz" comma separated; rows as samples; columns as mitochondrial probe L2R; no header or sample IDs)
3. Path to the **Sample descriptor file** ("EXAMPLE.SAMPLE_FILE.csv.gz"; comma separated; rows as samples; header line) containing the following columns:
* Column 1: "SAMPLE" --> unique sample IDs **in the order that they appear in the MT and Autosomal L2R matrices**
* Column 2: "SEX" --> self-reported gender (male=1;female=2) 
* Column 3: "AGE" --> age in years 
* Column 4 (Optional): "MITO_BENCHMARK" --> secondary mtDNA-CN measurement for comparison (e.g. qPCR, WGS, digital PCR). Missing values are perimssible and should be coded as 'NA'. 
4. **Number of cores** to parallelize the PCA steps (10)
5. Output Path / **Output File Name Prefix** ("TOY_EXAMPLE")

For file formatting examples, please download the toy dataset in the subsequent "Quick Start" section. 

## Quick Start

1. Download the R package and decompress:
```sh
git clone https://github.com/GMELab/AutoMitoC
cd AutoMitoC/scripts
```

2. Install requisite R libraries:
```R
R
install.packages(c("data.table","ggplot2","parallel"))
q()
n
```
3. Run the R script as follows:
```sh
Rscript AUTOMITOC_V1_RSCRIPTS_CLEAN.r EXAMPLE.AUTO_MATRIX.csv.gz EXAMPLE.MITO_MATRIX.csv.gz EXAMPLE.SAMPLE_FILE.csv.gz 10 EXAMPLE_TOY
```
4. If AutoMitoC is run successfully on the toy dataset, you should see something similar to the following output:

```R
[1] "AutoMitoC Script v1.0 Initialized @ 2021-11-26 02:26:47"
Loading required package: data.table
Loading required package: ggplot2
Loading required package: parallel
[1] "Reading in Input Arguments..."
[1] "Step 1. Preliminary PCA Correction of Autosomal Probe Intensities ... "
[1] "The Top 74 Autosomal Principal Components Explain ~70% variance in Probe Intensities and will be corrected for "
null device
          1
[1] "... Running Probe Correction with 10 core(s)"
[1] "Step 2. Probe Cross-hybridization Check ... "
[1] "Removing 127 Autosomal probes with potential cross-hybridization to the sex (R > 0.05) or mitochondrial (R > 0.05) genomes "
[1] "Removing 3 Mitochondrial probes with potential cross-hybridization to the sex (R > 0.20) or autosomal (R > 0.05) genomes "
[1] "*** Note that due to the robust epidemiological association between mtDNA-CN and sex, there is an expectation that mitochondrial probes will be associated with sex and therefore the correlation coefficient threshold for removing mitochondrial probes with evidence of cross-hybridization to sex chromosomes is more stringent"
[1] "Step 3. Final PCA of Clean Autosomal Probe Intensities ..."
[1] "The Top 72 Autosomal Principal Components Explain ~70% variance in Probe Intensities and will be corrected for "
null device
          1
[1] "... AutoMitoC mtDNA-CN estimates are finished!"
[1] "... Benchmarking Complete!"
[1] "AutoMitoC vs. user-supplied mtDNA-CN measurement performance: R = 0.647; Association P-value = 1.67e-164"
[1] "AutoMitoC Script Finished @ 2021-11-26 02:27:02"
[1] "Script Duration: 0 hrs, 0 mins, 16 secs ..."

```
## Areas of Future Development

AutoMitoC remains under active development and we plan to make improvements in the following areas in the near future:

1. The adjustment for autosomal principal components selects _k_ top PCs accounting for 70% of the variance in autosomal signal intensities, which approximated the inflection point on the scree plot and  worked well for us in practice on two independent datasets and different arrays; however, in the future, we plan to make this threshold less arbitrary and empirically determine the _k_ PCs corresponding to the inflection point. 
2. AutoMitoC has been benchmarked using multiple Affymetrix arrays and we plan to also test using Illumina arrays. 
3. We will allow for greater flexibility to enable benchmarking of AutoMitoC estimates to other phenotypic correlates of mtDNA-CN (e.g. blood cell composition)

## Contact information
Any queries pertaining to the AutoMitoC scripts can be addressed to:
**Michael Chong (michael.chong@phri.ca)** and **Guillaume Paré (pareg@mcmaster.ca)**

## Citation 
### Logo 

Created with [BioRender.com](https://biorender.com/)

### Manuscript

[Chong M, Mohammadi-Shemirani P, Perrot N, Nelson W, Morton R, Narula S, Lali R, Khan I, Khan M, Judge C, Machipisa T, Cawte N, O'Donnell M, Pigeyre M, Akhabir L, Paré G. GWAS and ExWAS of blood Mitochondrial DNA copy number identifies 71 loci and highlights a potential causal role in dementia. _eLife_ (2022) doi:10.7554/eLife.70382](https://elifesciences.org/articles/70382)

## Data Availability 
We applied AutoMitoC to the UK Biobank and performed GWAS for buffy coat mtDNA-CN including up to 395,781 participants. Full genome-wide summary statistics (even the sub-bonferroni associations!) are available for download from the [GWAS catalogue](https://www.ebi.ac.uk/gwas/publications/35023831)
