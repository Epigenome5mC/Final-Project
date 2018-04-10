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

You should create a SampleSheet, which explains the way the data is organized. This is extremely important. For my data, it was organized the following way:

| Sample_Name | Well_Position | Sentrix_Position | Sentrix_ID | Complete_Barcode |
| --- | --- | --- | --- | --- |
| MCF10A N 2d | A01 | _R01C01 | 201172520042 | 201172520042_R01C01|
| MCF10A H 0.5% 2d | B01 | _R02C01 | 201172520042 | 201172520042_R02C01 |
| 68 N 2d |	C01 |	_R03C01 |	201172520042 |	201172520042_R03C01 |
| 68 H 4% 2d |	D01 |	_R04C01 |	201172520042 |	201172520042_R04C01 |
| 68 H 0.5% 2d |	E01 |	_R05C01 |	201172520042 |	201172520042_R05C01 |
| 68 N 8 wk |	F01 |	_R06C01 |	201172520042 |	201172520042_R06C01 |
| 68 N 3d |	G01 |	_R07C01 |	201172520042 |	201172520042_R07C01 |
| 68 H 4% 8wk |	H01 |	_R08C01 |	201172520042 |	201172520042_R08C01 |


Creating the Basement, which will eventually create the RGset data set.

```
targets2 <- read.csv(file.path(baseDir, "SampleSheet.csv"), stringsAsFactors = FALSE)
targets2$Basename <- file.path(baseDir, targets2$Sentrix_ID, paste0(targets2$Sentrix_ID, targets2$Sentrix_Position))
targets <- targets2
RGset <- read.metharray.exp(targets = targets)
```
