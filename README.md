# Final-Project
Final Project for TRGN510

## Description of Project

My project consists of looking into the DNA methylation patterns of various circulating tumor cell (CTC) lines. I will be using data from my lab that was investigated by using an EPIC (850K) methylation array. This array interrogates about 850,000 CpG sites that focus on primarily promoters and enhancers. I will be comparing how the methylation patterns of all our CTC lines compare to each other. Our lines display various degrees of heterogeneity, but I am interested in identifying this same phenomenon in our CTC's methylomes. I will be using the R package called minfi, with the preprocessing steps of ssnoob (for EPIC data). My major goal for this project is to look at the DMPs (differentially methylated probes).

### 1. The necessary packages and tools.

This will load the required packages for this project.
```
source('http://bioconductor.org/biocLite.R')
biocLite('minfi')
biocLite("IlluminaHumanMethylationEPICmanifest")
biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
require(minfi)
```
### 2. Due to our data not being published, I have provided practice data which you can use in order to practice the following steps. I have also provided the data for veiwing purposes within this repository, so feel free to download it.

```
biocLite('minfidata')
```

### 3. Creating your directory.

#### This is for your own data. With this command, you can add a path directly from your computer.
```
baseDir <- file.path("~/Desktop/CTCs_EPIC/")
```

#### This will set up a directory for the minfi data for practice.
```
baseDir <- system.file("extdata", package = "minfiData")
```

### 4. Setting up the SampleSheet. 

I highly suggest you create a SampleSheet, which explains the way the data is organized. This is extremely important. If a sample sheet is not provided, you will not be able to proceed forward. I made mine like below and saved it as a .csv file. The following is an example from minfi's data, which I used to organize mine:

| Sample_Name | Sample_Well | Sample_Group | Sentrix_ID | Sentrix_Position | 
| --- | --- | --- | --- | --- |
| GroupA_3 | H5 | GroupA | 5723646052 | R02C02 |
| GroupA_2 | D5 | GroupA | 5723646052 | R04C01 |
| GroupB_3 | C6 |	GroupB | 5723646052 |	R05C02 |
| GroupB_1 | F7 |	GroupB | 5723646052 |	R04C02 |
| GroupA_1 | G7 |	GroupA | 5723646052 |	R05C02 |
| GroupB_2 | H7 |	GroupB | 5723646052 |	R06C02 |


### 5. Setting up your data and the Basement.
```
targets <- read.metharray.sheet(baseDir)
```
#### This will help you see if the Basement has been created.
```
sub(baseDir, "", targets$Basename)
```

#### This will create a Basement if you are unable to later. 
This will help in createing the RGset data set. The Basename column is vital because the column you make will contain the paths for each .idat file. You cannot generate this on your own through excel. You must run the code below to generate it. Make sure that where the .idat files are, you put them in a folder called by /organized by their 
Sentrix_ID (for example: 201172520042).

```
targets2 <- read.csv(file.path(baseDir, "SampleSheet.csv"), stringsAsFactors = FALSE)
targets2$Basename <- file.path(baseDir, targets2$Sentrix_ID, paste0(targets2$Sentrix_ID, targets2$Sentrix_Position))
targets <- targets2
```

### 6. Reading the data.
```
RGset <- read.metharray.exp(targets = targets)
```

#### Checking if annotation created.
```
RGset
```
#### If annotation not created, load:
```
RGset@annotation=c(array='IlluminaHumanMethylationEPIC', annotation='ilm10b2.hg19')
```

#### Preparing for PDF.
```
pd <- pData(RGset)
pd[,1:4]
```

### 7. PDF containing Quality Control data:

```
qcReport(RGset, sampNames = pd$Sample_Name, sampGroups = pd$Well_Position, pdf = "qcReport.pdf")
```

#### In the event that you want this data separately, I have provided the code for each part:

##### Beta density plots
```
densityPlot(RGset, sampGroups = pd$Sample_Group, main = "Beta", xlab = "Beta")
```

##### Beta beanplots
```
densityBeanPlot(RGset, sampGroups = pd$Sample_Group, sampNames = pd$Sample_Name)
```

##### Beta stripplot - Control: BISULFITE CONVERSION
```
controlStripPlot(RGset, controls="BISULFITE CONVERSION II", sampNames = pd$Sample_Name)
```

### 8. Normalizing the data.
```
MSet.raw <- preprocessRaw(RGset)
```

##### Normalizing with ssNoob (for EPIC data).
```
MSet.norm <- preprocessNoob(RGset, offset = 15, dyeCorr = TRUE, verbose = TRUE, dyeMethod=c("single", "reference"))
```

##### Normalizing with SWAN for 450K and the minfidata.
```
MSet.norm <- preprocessSWAN(RGset)
```


### 9. Creating a MDS plot with top 1,000 positions.
```
mdsPlot(MSet.norm, numPositions = 1000, sampGroups = pd$Sample_Group, sampNames = pd$Sample_Name)
```

### 10. Finding differentially methylated positions (DMPs).
This will create a 20,000 CpG subset of our dataset to speed up the demo:
```
mset <- MSet.norm[1:20000,]
```

##### This will load how many groups you decided to look at.
```
table(pd$Sample_Group)
```

##### Settig up data for finding DMPs.
```
M <- getM(mset, type = "beta", betaThreshold = 0.001)
dmp <- dmpFinder(M, pheno=pd$Sample_Group, type="categorical")
cpgs <- rownames(dmp)[1:4]
par(mfrow=c(2,2))
plotCpg(mset, cpg=cpgs, pheno=pd$Sample_Group)
```

### 11. Uploading data with ShinyMethyl.
Please make sure to run the following first:
```
source("https://bioconductor.org/biocLite.R")
biocLite("shinyMethyl")
biocLite("minfi")
```
Then run this on a R markdown using R studio:

```
library(minfi)
library(shinyMethyl)
baseDir <- file.path("~/Desktop/CTCs_EPIC/")
targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets)
summary <- shinySummarize(RGset)
runShinyMethyl(summary)
```





