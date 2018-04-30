# Final-Project
Final Project for TRGN510

## Description of Project

My project consists of looking into the DNA methylation patterns of various circulating tumor cell (CTC) lines. I will be using data from my lab that was investigated by using an EPIC (850K) methylation array. This array interrogates about 850,000 CpG sites that focus on primarily promoters and enhancers. I will be comparing how the methylation patterns of all our CTC lines compare to each other. Our lines display various degrees of heterogeneity, but I am interested in identifying this same phenomenon in our lines. I will be using the R package called minfi, with the preprocessing steps of ssnoob (for EPIC data). My major goal for this project is to create a heat map of the top DMPs (differentially methylated probes).

### 1. The necessary packages and tools.

This will load the required packages for this project.
```
source('http://bioconductor.org/biocLite.R')
biocLite('minfi')
biocLite('minfidata')
biocLite("IlluminaHumanMethylationEPICmanifest")
biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
require(minfi)
```
### 2. Due to our data not being published, I have provided practice data which you can use in order to practice the following steps. I have also provided the data for veiwing purposes within this repository, so feel free to download it.

```
biocLite('minfidata')

```

### 3. Creating your directory.

#### By using your own data. With this command, you can add a path directly from your computer.
```
baseDir <- file.path("~/Desktop/HypNorEPIC")
```

#### This will set up a directory for the minfi data for practice.
```
baseDir <- system.file("extdata", package = "minfiData")
```

You should create a SampleSheet, which explains the way the data is organized. This is extremely important. I made mine like below and saved it as a .csv file. For my data, it was organized the following way:

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
The Basename column is important because the column you make will contain the paths for each .idat file. You cannot generate this on your own through excel. You must run the code below to generate it. Make sure that where the .idat files are, you put them in a folder called by /organized by their 
Sentrix_ID (for example: 201172520042).

```
targets2 <- read.csv(file.path(baseDir, "SampleSheet.csv"), stringsAsFactors = FALSE)
targets2$Basename <- file.path(baseDir, targets2$Sentrix_ID, paste0(targets2$Sentrix_ID, targets2$Sentrix_Position))
targets <- targets2
RGset <- read.metharray.exp(targets = targets)
```

Checking things look ok.

```
pd <- pData(RGset)
pd[,1:4]
```

Set up the annotation of the data.

```
RGset@annotation=c(array='IlluminaHumanMethylationEPIC', annotation='ilm10b2.hg19')
```

This will give you a PDF containing Quality Control data:

```
qcReport(RGset, sampNames = pd$Sample_Name, sampGroups = pd$Well_Position, pdf = "qcReport.pdf")
```

Processing the data (normalizing):

```
MSet.raw <- preprocessRaw(RGset)
MSet.norm <- preprocessNoob(RGset, offset = 15, dyeCorr = TRUE, verbose = TRUE, dyeMethod=c("single", "reference"))
```

Creating a mdsPlot:

```
mdsPlot(MSet.norm, numPositions = 1000, sampGroups = pd$Well_Position, sampNames = pd$Sample_Name)
```




