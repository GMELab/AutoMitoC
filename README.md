# AutoMitoC
Automatic Mitochondrial DNA Copy Number (AutoMitoC) Pipeline to infer mtDNA-CN from SNP array intensities 

## Background

Mitochondria have their own unique DNA due to their origins as ancient bacteria. Due to this property, the ratio of mitochondrial to nuclear DNA signals, dubbed as the mitochondrial DNA copy number (mtDNA-CN) has emerged as a popular disease biomarker of mitochondrial etiology. One way to quantify mtDNA-CN is through genotyping arrays which enables not only the estimation of mtDNA-CN but also the capability to simultaneously interrogate its common genetic determinants. Previous array-based methods involve steps that are subjective, that require manual curation or extra lab testing, and that have not been extensively validated in different ethnicities. Accordingly, this served as the motivation to develop the AutoMitoC pipeline in order to enable fast, automated, and robust estimation of mtDNA-CN from genotyping arrays in all populations. 

## Key Features

AutoMitoC builds off of the MitoPipeline (Lane _et al_, 2014) with three key distinguishing features: 

1. Autosomal signal normalization utilises globally rare variants in place of common variants which confers advantages in terms of both speed and portability to ethnically diverse studies
2. Cross-hybridizing probes are identified by empirically assessing evidence for cross-hybridization via association of signal intensities
3. (iii) The primary estimate of mitochondrial (MT) signal is ascertained using principal component analysis (PCA) as opposed to using the median signal intensity of MT probes. For further details on how AutoMitoC was developed or the rationale behind these features, please see the main manuscript (link).

## Method Overview

![image](https://user-images.githubusercontent.com/30928727/143525953-4f39541d-53e0-4f3e-a5bf-4850ad2f1b10.png)

## Required Inputs ##

## ## 


## Dependencies ## 

##  ##

## Example with toy dataset ## 

## Caveats ##

AutoMitoC remains under active development and we encourage users to provide feedback to continually improve the method. That being said, we plan to make improvements in the following areas in the near future:

1. The adjustment for autosomal principal components selects _k_ top PCs accounting for 70% of the variance in autosomal signal intensities, which approximated the inflection point on the scree plot and also worked well for us in practice on two independent datasets and two different arrays; however, in the future, we plan to make this threshold less arbitrary and empirically determine the _k_ PCs corresponding to the inflection point. 
2. AutoMitoC has been benchmarked using multiple Affymetrix arrays and we plan to validate using Illumina arrays next. 
