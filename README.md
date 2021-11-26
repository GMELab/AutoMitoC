# AutoMitoC
Automatic Mitochondrial DNA Copy Number (AutoMitoC) Pipeline to infer mtDNA-CN from SNP array intensities 

## Background

Mitochondria have their own unique DNA due to their origins as ancient bacteria. Due to this property, the ratio of mitochondrial to nuclear DNA signals, dubbed as the mitochondrial DNA copy number (mtDNA-CN) has emerged as a popular disease biomarker of mitochondrial etiology. One way to quantify mtDNA-CN is through genotyping arrays which enables not only the estimation of mtDNA-CN but also the capability to simultaneously interrogate its common genetic determinants. Previous array-based methods involve steps that are subjective, that require manual curation or extra lab testing, and that have not been extensively validated in different ethnicities. Accordingly, this served as the motivation to develop the AutoMitoC pipeline in order to enable fast, automated, and robust estimation of mtDNA-CN from genotyping arrays in all populations. 

## Key Features

AutoMitoC builds off of the MitoPipeline (Lane _et al_, 2014) with three key distinguishing features: 

* Nuclear signal is approximated by globally rare variants in place of common variants 
* Cross-hybridizing probes are identified by empirically assessing evidence for cross-hybridization
* The primary estimate of mitochondrial (MT) signal is ascertained using principal component analysis (PCA)
 
For further details on how AutoMitoC was developed or the rationale behind these features, please see the main manuscript (link).

## Method Overview

![image](https://user-images.githubusercontent.com/30928727/143525953-4f39541d-53e0-4f3e-a5bf-4850ad2f1b10.png)

## Necessary Input Files

The AutoMitoC pipeline requires the following inputs.

1. MT L2R matrix (comma separated; rows as samples; columns as mitochondrial probes; no header or sample IDs)
2. Autosomal L2R matrix (comma separated; rows as samples; columns as mitochondrial probes; no header or sample IDs)
3. A sample descriptor file (comma separated; rows as samples; header line) containing:
i) sample IDs #in the order that they appear in the MT and Autosomal L2R matrices
ii) self-reported gender
iii) age in years
iv) (optional) secondary mtDNA-CN measurement (e.g. qPCR, WGS, digital PCR)

For file formatting examples, please download the toy dataset in the subsequent "Quick Start" section. 





Currently, the AutoMitoC pipeline starts at the second major step after pre-processing. Pre-processing can be done using previously developed tools.

For Affymetrix data, L2R values can be derived by XXX

Filtering for rare autosomal variants can be done as follows. 

Filtering for common autosomal variants can be done as follows. 

## Quick Start

1. Download the toy dataset  
2. Install requisite R libraries

install.packages("data.table")
    install.packages("ggplot2")
    install.packages("parallel")

3. Run the R script as follows

## Caveats and Areas of Future Development

AutoMitoC remains under active development and we encourage users to provide feedback to michael.chong@phri.ca for suggestions. We plan to make improvements in the following areas in the near future:

1. The adjustment for autosomal principal components selects _k_ top PCs accounting for 70% of the variance in autosomal signal intensities, which approximated the inflection point on the scree plot and  worked well for us in practice on two independent datasets and different arrays; however, in the future, we plan to make this threshold less arbitrary and empirically determine the _k_ PCs corresponding to the inflection point. 
2. AutoMitoC has been benchmarked using multiple Affymetrix arrays and we plan to soon test using Illumina arrays. 
3. In the future, we will allow for greater flexibility to enable benchmarking of AutoMitoC estimates to other phenotypic correlates of mtDNA-CN (e.g. blood cell composition)
