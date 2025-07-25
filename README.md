# Chapter 1_Script
Chapter 1 analysis with Dart data


# 1 Data wrangling and filtering


# Description

The dataset analysed in the following script is from this publication:
Erdmann et al., 2025. Population structure of Acropora kenti at inshore reefs of the Great Barrier Reef. Molecular Ecology.

This script analyses a SNP's dataset from DArT using the package darTRverse (version 2024), which contains 7 subpackages:
https://github.com/green-striped-gecko/dartRverse

# Abbreviations

Reef names:
BRR = Bramble Reef
EP = East Pelorus
GB = Geoffrey Bay
HaR = Havannah Reef
HB = Horseshoe Bay
HFB = Huntingfield Bay
JB = John Brewer Reef
KR = Keeper Reef
LPB = Little Pioneer Bay
MB = Maud Bay
MR = Middle Reef
PB = Picnic Bay
SO = South Orpheus
WB = Wilson Bay
WP = West Pelorus

Data sets:
ak = Acropora kenti SNP dataset
ak.filtered = used ak, which has been filtered for all downstream analysis
ak.filtered.nr = used ak.filtered and removed replicates (nr = no replicates)
ak.filtered.nr.nm = = used ak.filtered.nr and removed migrants (nm = no migrants)
ak.pop = used ak.filtered.nr to run population structure analysis
ak.pop.nm = used ak.pop and removed migrants
ak.pop.mag = population structure dataset for Magnetic Island only
ak.rel = used ak.filtered.nr.nm to run relatedness analysis
ak.gen = used ak.filtered.nr.nm to calculate genetic diversity estimates (Fst, Fis, He, Ho)
ak.gen.mag = used ak.gen to calculate genetic diversity metrics for Magnetic Island only
ak.gen.adj = used ak.gen to calculate genetic diversity metrics for adjacent reefs only

# Installation

https://green-striped-gecko.github.io/kioloa/install.html

```{r}
# Install the necessary bioconductor packages
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("SNPRelate")

# Install dartRverse (dartRverse) & core (dartR.base, dartR.data)
install.packages("dartRverse")
```

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("snpStats")
```

```{r}
# Install dartR packages
install.packages("dartR.sim")
install.packages("dartR.popgen")
install.packages("dartR.spatial")
install.packages("dartR.captive")

# Update dartR.base
dartRverse_install(package = "dartR.base", rep = "github", branch = "dev")
# Check which packages and versions are installed
dartRverse_install()
```

## Load library

```{r}
library(dartRverse)
```

# Load and inspect raw dataset

The SNP input dart file is named Report_DAc23-8597_5_moreOrders_SNP_2.csv.
The metadata file is called "ind_metrics_SNP2.csv" and contains information for all individuals, such as, ID, name of populations, locations and coordinates of sample locations.

```{r}
# Merge the DartR SNP dataset with the metadata file
ak <- gl.read.dart(filename="Report_DAc23-8597_5_moreOrders_SNP_2.csv", ind.metafile = "ind_metrics_SNP2.csv")
```

```{r}
# Save merged SNP and metadata file
save(ak, file="ak.rdata")    
```

```{r}
# Get names of individuals, location and population
nInd(ak)
nLoc(ak)
nPop(ak)
```

```{r}
# Create a table showing the number of individuals per population
table(pop(ak))
```

```{r}
# Create a barplot from that table
barplot(table(pop(ak)), las=2)
```

```{r}
# Run a smearplot to assess quality of dataset
gl.smearplot(ak)
```

# Remove duplicates

```{r}
# Get names of individual samples
ak$ind.names
```

```{r}
# Inspect the individual call rates and remove duplicates with the lower call rate
gl.report.callrate(ak, method = "ind", ind.to.list = 528)
```

Remove the following duplicates:
At_162_J20.1', 'At_163_J20.1', 'B1', 'A2', 'A3', 'A7', 'B8', 'A9', 'A10', 'B11', 'A15', 'B17', 'B18', 'A19', 'A21', 'WB_03/04/2022_36'

```{r}
# Drop duplicates from current dataset
ak1 <- gl.drop.ind(ak, ind.list=c('At_162_J20.1', 'At_163_J20.1', 'B1', 'A2', 'A3', 'A7', 'B8', 'A9', 'A10', 'B11', 'A15', 'B17', 'B18', 'A19', 'A21', 'WB_03/04/2022_36'))
```

```{r}
# Filter for monomorphs after removing duplicates
ak1 <- gl.filter.monomorphs(ak1)
```

```{r}
# Recalculate locus metrics
ak1 <- gl.recalc.metrics(ak1)
```

# Filtering

## ak.clean filtering

Filter dataset for (1) reproducibility 0.99, (2) read depth >10 and <100, (3) maf (0.005), (4) call rate for loci (0.8), (5) call rate for individuals (0.75), (6) monomorphs, (7) secondaries and (8) remove missing loci.

```{r}
# Run a smearplot to assess quality of dataset
gl.smearplot(ak1)
```

```{r}
nInd(ak1)
nLoc(ak1)
```

### Reproducibility 0.99

A metrics that investigates if calling a SNP twice results in the same outcome.

```{r}
# Inspect the current reproducibility to remove unreliable loci
gl.report.reproducibility(ak1)
```

```{r}
# Filter for reproducibility due to high missing rate
ak2 <- gl.filter.reproducibility(ak1)
```

```{r}
# Run a smearplot to assess quality of dataset
gl.smearplot(ak2)
```

```{r}
nInd(ak2)
nLoc(ak2)
```

### Read depth >10 <100

Loci with high (= unreliable) and low read depth (= non-confident) should be removed.

```{r}
# Inspect the current read depth
gl.report.rdepth(ak2)
```

```{r}
# Filter for a read depth greater than 10 and lower than 100.
ak3 <- gl.filter.rdepth(ak2, lower = 10, upper = 100)
```

```{r}
nInd(ak3)
nLoc(ak3)
```

```{r}
# Run a smearplot to assess quality of dataset
gl.smearplot(ak3)
```

### Minor allele frequency (maf) 0.005

This filter removes singletons, especially if read depth is low. Despite that our read depth is high, we apply this filter with a low value to keep more loci, but improve the signal in our populations structure analysis.

```{r}
# Inspect the current minor allele frequency
gl.report.maf(ak3)
```

```{r}
# Filter for a minor allele frequency of 0.005
ak4 <- gl.filter.maf(ak3,  threshold = 0.005)
```

```{r}
nInd(ak4)
nLoc(ak4)
```

```{r}
# Run a smearplot to assess quality of dataset
gl.smearplot(ak4)
```

### Call rate for loci 0.8

A call rate for loci of 0.8 removes very poorly sequenced loci.

```{r}
# Inspect the current call rate for loci
gl.report.callrate(ak4)
```

```{r}
# Filter fora call rate of 0.8
ak5 <- gl.filter.callrate(ak4, method = "loc", threshold = 0.8)
```

```{r}
nInd(ak5)
nLoc(ak5)
```

```{r}
# Run a smearplot to assess quality of dataset
gl.smearplot(ak5)
```

### Call rate for individual 0.75

A call rate for individuals of 0.75 removes very poorly sequenced loci.

```{r}
# Inspect the current callrate for individuals
gl.report.callrate(ak5, method="ind")
```

```{r}
# Filter for a call rate for individuals of 0.75
ak6 <-gl.filter.callrate(ak5, method ="ind", threshold=0.75)
```

```{r}
# Recalculate locus metrics
ak6 <- gl.recalc.metrics(ak6)
```

```{r}
nInd(ak6)
nLoc(ak6)
```

```{r}
# Run a smearplot to assess quality of dataset
gl.smearplot(ak6)
```

### Monomorphs

This step deletes monomorphic loci from a genlight {adegenet} object. A DArT dataset will not have monomorphic loci, but they can arise, along with loci when populations or individuals are deleted.

```{r}
# Inspect the current number of monomorphs
gl.report.monomorphs(ak6)
```

```{r}
# Filter for monomorphs
ak7 <- gl.filter.monomorphs(ak6 ,verbose = 5)
```

```{r}
nInd(ak7)
nLoc(ak7)
```

```{r}
# Run a smearplot to assess quality of dataset
gl.smearplot(ak7)
```

### Secondaries

Secondaries are fragments with more than one SNP. These multiple SNP loci within a fragment are likely to be linked, and so you may wish to remove secondaries. The filter removes loci that may not be inherited independently. Filtering secondaries at the end of the filtering process ensures that highly quality loci remain in the dataset.

```{r}
# Inspect the current number of secondaries
gl.report.secondaries(ak7)
```

```{r}
# Filter for secondaries
ak8 <- gl.filter.secondaries(ak7, verbose = 5)
```

```{r}
nInd(ak8)
nLoc(ak8)
```

```{r}
# Run a smearplot to assess quality of dataset
gl.smearplot(ak8)
```

### Remove missing loci

A DArT dataset will not have individuals for which all calls are scored as missing (NA) across all loci, but such individuals may sneak in to the dataset when loci are deleted. Retaining individual or loci with all NAs can cause issues for several functions.

```{r}
# Remove missing loci
ak9 <-gl.filter.allna(ak8)
```

```{r}
nInd(ak9)
nLoc(ak9)
```

```{r}
# Run a smearplot to assess quality of dataset
gl.smearplot(ak9)
```

```{r}
# Recalculate locus metrics
ak9 <- gl.recalc.metrics(ak9)
```

# Inspect filtered dataset

```{r}
# A compliance check ensures the dataset has no monomorphic loci, no missing data, properly assigned populations, individual names are unique, coordinates are correct, and recalculates locus metrics.
gl.compliance.check(ak9)
```

```{r}
nInd(ak9)
nLoc(ak9)
nPop(ak9)
```

```{r}
# Create a table showing the number of individuals per population
table(pop(ak9))
```

```{r}
# Create a barplot from that table
barplot(table(pop(ak9)), las=2)
```

```{r}
# rename the dataset after the filtering process is completed
ak.filtered <- ak9
```

```{r}
# Get names of individuals
ak.filtered$ind.names
```

```{r}
# Rename sample "20" to WB_36 for consistency
ak.filtered$ind.names[ak.filtered$ind.names == "20"] <- "WB_36"
```

```{r}
# Confirm name change
ak.filtered$ind.names
```

```{r}
# Save filtered and clean dataset
save(ak.filtered, file="ak.filtered.rdata")
```

## Assign colours to populations

Assign colours based on these site locations:
Magnetic Island vs. Palm Island vs. Mid-shelf

```{r}
# Get order of population names
#levels(ak.gen.mi@pop)
#  "BRR" "EP"  "MI"  "HaR" "JB"  "KR"  "LPB" "SO"  "WP"

# Assign colours to populations in that order:
# Mid-shelf = yellow "#FDD835"; Palms = green "#74c476"; Magnetic Island = orange "orange";
cols.mi <- c("#FDD835", "#74c476", "orange", "#74c476", "#FDD835", "#FDD835", "#74c476", "#74c476", "#74c476")
```

```{r}
# Get order of population names
#levels(ak.gen.mi.pi@pop)
#  "BRR" "PI"  "MI"  "JB"  "KR"

# Assign colours to populations in that order: 
# Mid-shelf = yellow "#FDD835"; Palms = green "#74c476"; Magnetic Island = orange "orange";
cols.mi.pi <- c("#FDD835", "#74c476", "orange", "#FDD835", "#FDD835")
```

```{r}
# Get order of population names
#levels(ak.gen.mi.pi.ms@pop)
#  "MS" "PI" "MI"

# Assign colours to populations in that order: 
# Mid-shelf = yellow "#FDD835"; Palms = green "#74c476"; Magnetic Island = orange "orange";
cols.mi.pi.ms <- c("#FDD835", "#74c476", "orange")
```
