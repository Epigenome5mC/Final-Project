# Final-Project
Final Project for TRGN510

## Description of Project

My project consists of looking into DNA methylation patterns of various cell lines. I will be using data from my lab that was investigated by using an EPIC (850K) methylation array. This array interrogates about 850,000 CpG sites that focus on primarily promoters and enhancers. Our data compares hypoxia of various oxygen concentrations vs. normoxic samples (21% Oxygen). I am curious to identify what genes fall under the differentially methylated probes (DMP) within this data. I will be using the R package called minfi, with the preprocessing steps of ssnoob (for EPIC data). My major goal for this project is to create a heat map of the top 1,000 DMPs, and if time allows, to identify the genes that fall under these probes and identify if they fall primarily at promoters or enhancers.

## R programming steps

This will load the required packages for this project.
```
source('http://bioconductor.org/biocLite.R')
biocLite('minfi')
biocLite("IlluminaHumanMethylationEPICmanifest")
biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
require(minfi)
```

Creating your directory

```
baseDir <- file.path("~/Desktop/HypNorEPIC")
```

Creating RGset data set.

```
targets2 <- read.csv(file.path(baseDir, "SampleSheet.csv"), stringsAsFactors = FALSE)
targets2$Basename <- file.path(baseDir, targets2$Sentrix_ID, paste0(targets2$Sentrix_ID, targets2$Sentrix_Position))
targets <- targets2
RGset <- read.metharray.exp(targets = targets)
```
