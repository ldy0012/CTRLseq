# CTRLseq

CTRLseq is a VAE-based method for latent confounder correction in RNA-seq differential expression analysis.

## Overview

RNA-seq data are often affected by hidden confounding factors such as batch effects, tissue heterogeneity, and technical variation.  
CTRLseq estimates these latent confounders directly from residual structure using a variational autoencoder (VAE) and incorporates them into differential expression models.

## Method

1. Fit GLM to RNA-seq count data  
2. Extract deviance residuals  
3. Estimate latent factors using VAE  
4. Incorporate latent factors into adjusted GLM  
5. Perform differential expression analysis  

## Installation

```r
devtools::install_github("ldy0012/CTRLseq")
