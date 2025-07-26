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


# 2 Population structure

## Load Libraries

```{r}
library(ggplot2)      # to plot results
library(ggpubr)       # to use ggarrange and other plotting functions
library(gridExtra)    # to use grid.arrange
library(cowplot)      # to use plot_grid and ggdraw
library(tidyr)        # general functions
```

## Identify close relatives/replicates

Identify close relatives with a percentage of 99% similarity in their genotypes. 

Ideally, in a large dataset with related and unrelated individuals and several replicated individuals, the first histogram should have four "peaks". The first peak should represent unrelated individuals, the second peak should correspond to second-degree relationships (such as cousins), the third peak should represent first-degree relationships (like parent/offspring and full siblings), and the fourth peak should represent replicated individuals.

In order to ensure that replicated individuals are properly identified, it's important to have a clear separation between the third and fourth peaks in the second histogram. This means that there should be bins with zero counts between these two peaks.

```{r}
# Inspect the number of replicates
dartR.base::gl.report.replicates(ak.filtered, perc_geno=0.99)
```

Drop these 24 individuals: "At_11_J20", "A154", "A24", "At_9_J20", "B9", "At_29_J20", "At_8_J20", "At_7_J20", "A99", "At_34_N", "At_4_F", "At_41_J20", "At_58_F20", "At_6_N", "A11", "A25", "A8", "GB_03/03/2022_17" "HB_4", "MR_14/03/2022_12" "MR_14/03/2022_13", "MR_14/03/2022_8", "PB_05/04/22_45", "WB_03/04/2022_14"

## Remove replicates

Close relatives (with a percentage of 99% similarity in their genotypes) were identified using gl.report.replicates and will be removed from this dataset.

Drop individuals from the replicates function:

```{r}
ak.filtered.nr <- gl.drop.ind(ak.filtered, ind.list=c("At_11_J20", "A154", "A24", "At_9_J20", "B9", "At_29_J20", "At_8_J20", "At_7_J20", "A99", "At_34_N", "At_4_F", "At_41_J20", "At_58_F20", "At_6_N", "A11", "A25", "A8", "GB_03/03/2022_17", "HB_4", "MR_14/03/2022_12", "MR_14/03/2022_13", "MR_14/03/2022_8", "PB_05/04/22_45", "WB_03/04/2022_14"))
```

```{r}
# Filter for monomorphs
ak.filtered.nr <- gl.filter.monomorphs(ak.filtered.nr)
```

```{r}
# Recalculate locus metrics
ak.filtered.nr <- gl.recalc.metrics(ak.filtered.nr)
```

```{r}
nInd(ak.filtered.nr)
nLoc(ak.filtered.nr)
```

```{r}
gl.smearplot(ak.filtered.nr)
```

## Remove missing loci

A DArT dataset will not have individuals for which the calls are scored as missing (NA) across all loci, but such individuals may sneak in to the dataset when loci are deleted. Retaining individual or loci with all NAs can cause issues for several functions

```{r}
# Filter for missing loci
ak.pop <-gl.filter.allna(ak.filtered.nr)
```

## Inspect filtered file

```{r}
nInd(ak.pop)
nLoc(ak.pop)
```

```{r}
table(pop(ak.pop))
```

```{r}
indNames(ak.pop)
```

```{r}
save(ak.pop, file="ak.pop.rdata")
```

## Convert to structure

Convert ak.pop to structure format to run externally with structure program

```{r}
gl2structure(ak.pop, ind.names = NULL, add.columns = NULL, ploidy = 2, export.marker.names = TRUE, outfile = "ak.pop.str", outpath = NULL, verbose = NULL)
```

```{r}
nInd(ak.pop)
nLoc(ak.pop)
```

## Admixture plot for ak.pop

### Load csv

Added the location to csv, so it is possible to group by location for plots.

```{r}
# Load csv with data from run # 12 in the structure analysis
run12 <-  read.csv("ak.pop_run12.csv")
```

```{r}
head(run12)
```

### Data wrangling

For the plots to show both clusters the matrix needs to be amended.

```{r}
# Create 2 rows per ID for each cluster, assign the fst value and arrange by ID
run12<-run12 %>% gather(key="cluster", value = "fst", c(-ID, -loc)) %>% arrange(ID)
```

```{r}
head(run12)
```

For a first glimpse, create a barplot that shows the bays on the x-axis, fst on y-axis and all 4 clusters stacked.

```{r}
strrun12 <- ggplot(run12, aes(x = ID, y = fst, fill=cluster)) +
   geom_bar(position="stack", stat = "identity", width=0.2)
strrun12
```

For a more detailed investigation into single bays, select the location from the dataset and make subsets for all locations ("HFB", "WB", "MB", "HB", "GB", "PB", "MR", "HaR", "SO", "LPB", "WP", "EP", "BRR", "JB", "KR")

```{r}
HFBrun12 <-run12[grep("^HFB", run12$loc), ] 
HFBrun12
```

```{r}
WBrun12<-run12[grep("^WB", run12$ID), ]
WBrun12
```

```{r}
MBrun12<-run12[grep("^MB", run12$ID), ]   
MBrun12
```

```{r}
HBrun12<-run12[grep("^HB", run12$ID), ]
HBrun12
```

```{r}
GBrun12<-run12[grep("^GB", run12$ID), ]
GBrun12
```

```{r}
PBrun12<-run12[grep("^PB", run12$ID), ]
PBrun12
```

```{r}
MRrun12<-run12[grep("^MR", run12$ID), ]
MRrun12
```

```{r}
HaRrun12<-run12[grep("^HaR", run12$ID), ]   
HaRrun12
```

```{r}
SOrun12 <-run12[grep("^SO", run12$loc), ] 
SOrun12
```

```{r}
LPBrun12 <-run12[grep("^LPB", run12$loc), ]  
LPBrun12
```

```{r}
WPrun12 <-run12[grep("^WP", run12$loc), ] 
WPrun12
```

```{r}
EPrun12 <-run12[grep("^EP", run12$loc), ] 
EPrun12
```

```{r}
BRRrun12<-run12[grep("^BRR", run12$ID), ]     #  samples
BRRrun12
```

```{r}
JBrun12<-run12[grep("^JB", run12$ID), ]  
JBrun12
```

```{r}
KRrun12<-run12[grep("^KR", run12$ID), ] 
KRrun12
```

### Plot

```{r}
# Assign a certain colour to the clusters, so that cluster1 is consistently blue and cluster2 orange
cluster_colors <- c("orange", "blue")
```

Make a plot for all locations ("HFB", "WB", "MB", "HB", "GB", "PB", "MR", "HaR", "SO", "LPB", "WP", "EP", "BRR", "JB", "KR").

```{r}
HFB12 <- ggplot(HFBrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
   theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
  ggtitle("Huntingfield Bay")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/HFB12.png",
width=800, height=400)

HFB12
```

```{r}
WB12 <- ggplot(WBrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
  ggtitle("Wilson Bay")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/WB12.png",
width=800, height=400)

WB12
```

```{r}
MB12 <- ggplot(MBrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
   theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Maud Bay")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/MB12.png",
width=800, height=400)

MB12
```

```{r}
HB12 <- ggplot(HBrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
    theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
  ggtitle("Horseshoe Bay")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/HB12.png",
width=800, height=400)

HB12
```

```{r}
GB12 <- ggplot(GBrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Geoffrey Bay")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/GB12.png",
width=800, height=400)

GB12
```

```{r}
PB12 <- ggplot(PBrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Picnic Bay")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/PB12.png",
width=800, height=400)

PB12
```

```{r}
MR12 <- ggplot(MRrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Middle Reef")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/MR12.png",
width=800, height=400)

MR12
```

```{r}
HaR12 <- ggplot(HaRrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors)  +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Havannah Reef")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/HaR12.png",
width=800, height=400)

HaR12
```

```{r}
SO12 <- ggplot(SOrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("South Orpheus")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/SO12.png",
width=800, height=400)

SO12
```

```{r}
LPB12 <- ggplot(LPBrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Little Pioneer Bay")


png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/LPB12.png",
width=800, height=400)

LPB12
```

```{r}
WP12 <- ggplot(WPrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
      ggtitle("West Pelorus")


png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/WP12.png",
width=800, height=400)

WP12
```

```{r}
EP12 <- ggplot(EPrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
     ggtitle("East Pelorus")


png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/EP12.png",
width=800, height=400)

EP12
```

```{r}
BRR12 <- ggplot(BRRrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Bramble Reef")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/BRR12.png",
width=800, height=400)

BRR12
```

```{r}
JB12 <- ggplot(JBrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("John Brewer Reef")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/JB12.png",
width=800, height=400)

JB12
```

```{r}
KR12 <- ggplot(KRrun12, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Keeper Reef")

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/KR12.png",
width=800, height=400)

KR12
```

Combine all structure plots into one plot.

```{r}
# Combine all Magnetic Island subpopulations
# Remove both x and y axis elements (including labels, ticks, and lines) from all but the first column and combine all plots
strmag12_combined <- ggarrange(
  HFB12 + theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()),  # Remove both axes
  WB12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  MB12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  HB12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  GB12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  PB12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  MR12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  ncol = 2, nrow = 4,
  widths = c(1, 1),
  heights = c(1, 1),
  legend = "none"
)

# Add a common y-axis label using cowplot::plot_grid()
strmag12_labeled <- plot_grid(
  ggdraw() + 
    draw_label("Admixture proportions", angle = 90, size = 15, vjust = 0.5),
  strmag12_combined,
  ncol = 2,
  rel_widths = c(0.05, 1)  # 5% for label, 95% for plot
)

# Save the plot
ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures//strmag12.pdf", plot = strmag12_labeled, width = 8, height = 10)
```

```{r}
# Combine all adjacent reef populations
# Remove both x and y axis elements (including labels, ticks, and lines) from all but the first column and combine all plots
stradj12_combined <- ggarrange(
  HaR12 + theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()),  # Remove both axes
  SO12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  LPB12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  WP12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  EP12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  BRR12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  JB12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  KR12 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  ncol = 2, nrow = 4,
  widths = c(1, 1),
  heights = c(1, 1),
  legend = "none"
)

# Add a common y-axis label using cowplot::plot_grid()
stradj12_labeled <- plot_grid(
  ggdraw() + 
    draw_label("Admixture proportions", angle = 90, size = 15, vjust = 0.5),
  stradj12_combined,
  ncol = 2,
  rel_widths = c(0.05, 1)  # 5% for label, 95% for plot
)

# Save the plot
ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures//stradj12.pdf", plot = stradj12_labeled, width = 8, height = 10)
```

# Fine-scale population structure for Magnetic Island

## Remove migrants

For fine scale structure at Magnetic Island migrants from adjacent reefs will be removed to create a subset only for the Magnetic Island population called ak.pop.nm.mag

```{r}
# Get names for all individuals
indNames(ak.pop)
```

The following migrants, which have allele frequency differences >0.5, have been identified with structure plots:
HB: 1, 2, 14, 20, 21, 2021_11
WB: 8, 25, 37
HFB: 4, 6, 9, 10, 22
HaR: 7, 10, 12, 24, 43

```{r}
# Remove 19 migrants with a proportion >5%
ak.pop.nm <- gl.drop.ind(ak.pop, ind.list=c("HB_1", "HB_2", "HB_27/11/21_11", "HB_27/11/21_14", "HB_27/11/21_20", "HB_27/11/21_21", "WB_03/04/2022_8", "WB_03/04/2022_25", "WB_03/04/2022_37", "HFB_4", "HFB_6", "HFB_9", "HFB_10", "HFB_22", "HaR_06/05/2022_7", "HaR_06/05/2022_10", "HaR_06/05/2022_12", "HaR_06/05/2022_24", "HaR_06/05/2022_43"))
```

```{r}
# Filter for monomorphs
ak.pop.nm <- gl.filter.monomorphs(ak.pop.nm)
```

```{r}
# Recalculate locus metrics
ak.pop.nm <- gl.recalc.metrics(ak.pop.nm)
```

## Remove missing loci

A DArT dataset will not have individuals for which the calls are scored all as missing (NA) across all loci, but such individuals may sneak in to the dataset when loci are deleted. Retaining individual or loci with all NAs can cause issues for several functions

```{r}
# Filter for missing loci
ak.pop.nm <- gl.filter.allna(ak.pop.nm)
```

```{r}
nInd(ak.pop.nm)
nLoc(ak.pop.nm)
```

```{r}
gl.smearplot(ak.pop.nm)
```

```{r}
save(ak.pop.nm, file="ak.pop.nm.rdata")
```

## Keep only Magnetic Island subpopulations

```{r}
indNames(ak.pop.nm)
```

```{r}
# Keep only subpopulations from Magnetic Island for fine structure analysis
ak.pop.mag <- gl.keep.pop(ak.pop.nm, pop.list = list("GB", "HB", "HFB", "MB", "MR", "PB", "WB"))
```

```{r}
save(ak.pop.mag, file="ak.pop.mag.rdata")
```

```{r}
nInd(ak.pop.mag)
nLoc(ak.pop.mag)
```

## Convert to structure

Convert ak.pop.mag to structure format to run externally with structure program

```{r}
gl2structure(ak.pop.mag, ind.names = NULL, add.columns = NULL, ploidy = 2, export.marker.names = TRUE, outfile = "ak.pop.mag.str", outpath = NULL, verbose = NULL)
```

```{r}
nInd(ak.pop.mag)
nLoc(ak.pop.mag)
```

## Admixture plot for ak.pop.nm.mag

### Load csv

Added the location to csv, so it is possible to group by location for plots.

```{r}
# Load csv with data from run # 12 in the structure analysis
run10 <-  read.csv("ak.pop.nm.mag_run10.csv")
```

```{r}
head(run10)
```

### Data wrangling

For the plots to show both clusters the matrix needs to be amended.

```{r}
# Create 2 rows per ID for each cluster, assign the fst value and arrange by ID
run10<-run10 %>% gather(key="cluster", value = "fst", c(-ID, -loc)) %>% arrange(ID)
```

```{r}
head(run10)
```

For a first glimpse, create a barplot that shows the bays on the x-axis, fst on y-axis and all 4 clusters stacked.

```{r}
strrun10 <- ggplot(run10, aes(x = ID, y = fst, fill=cluster)) +
   geom_bar(position="stack", stat = "identity", width=0.2)
strrun10
```

For a more detailed investigation into single bays, select the location from the dataset and make subsets for all locations ("HFB", "WB", "MB", "HB", "GB", "PB", "MR").

```{r}
HFBrun10 <-run10[grep("^HFB", run10$loc), ] 
HFBrun10
```

```{r}
WBrun10<-run10[grep("^WB", run10$ID), ]
WBrun10
```

```{r}
MBrun10<-run10[grep("^MB", run10$ID), ]   
MBrun10
```

```{r}
HBrun10<-run10[grep("^HB", run10$ID), ]
HBrun10
```

```{r}
GBrun10<-run10[grep("^GB", run10$ID), ]
GBrun10
```

```{r}
PBrun10<-run10[grep("^PB", run10$ID), ]
PBrun10
```

```{r}
MRrun10<-run10[grep("^MR", run10$ID), ]
MRrun10
```

### Plot

```{r}
# Assign a certain colour to the clusters, so that cluster1 is consistently blue and cluster2 orange
cluster_colors <- c("blue", "orange")
```

Make a plot for all locations ("HFB", "WB", "MB", "HB", "GB", "PB", "MR", "HaR", "SO", "LPB", "WP", "EP", "BRR", "JB", "KR").

```{r}
HFB10 <- ggplot(HFBrun10, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
  ggtitle("Huntingfield Bay") +
   theme_classic()

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/HFB10.png",
width=800, height=400)

HFB10
```

```{r}
WB10 <- ggplot(WBrun10, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
  ggtitle("Wilson Bay") +
    theme_classic()

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/WB10.png",
width=800, height=400)

WB10
```

```{r}
MB10 <- ggplot(MBrun10, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Maud Bay") +
   theme_classic()

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/MB10.png",
width=800, height=400)

MB10
```

```{r}
HB10 <- ggplot(HBrun10, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
  ggtitle("Horseshoe Bay") +
    theme_classic()

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/HB10.png",
width=800, height=400)
HB10
```

```{r}
GB10 <- ggplot(GBrun10, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position="stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Geoffrey Bay") +
  theme_classic()

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/GB10.png",
width=800, height=400)
GB10
```

```{r}
PB10 <- ggplot(PBrun10, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Picnic Bay") +
  theme_classic()

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/PB10.png",
width=800, height=400)

PB10
```

```{r}
MR10 <- ggplot(MRrun10, aes(x = ID, y = fst, fill = cluster)) +
   geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90)) +
  font("x.text", size = 8) +
    labs(y = "admixture proportions") +
    ggtitle("Middle Reef") +
  theme_classic()

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/MR10.png",
width=800, height=400)
MR10
```

Combine all structure plots into one plot.

```{r}
# Combine all Magnetic Island subpopulations
# Remove both x and y axis elements (including labels, ticks, and lines) from all but the first column and combine all plots
strmaggie10_combined <- ggarrange(
  HFB10 + theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()),  # Remove both axes
  WB10 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  MB10 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  HB10 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  GB10 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  PB10 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  MR10 + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank()),  # Remove both axes
  ncol = 2, nrow = 4,
  widths = c(1, 1),
  heights = c(1, 1),
  legend = "none"
)

# Add a common y-axis label using cowplot::plot_grid()
strmaggie10_labeled <- plot_grid(
  ggdraw() + 
    draw_label("Admixture proportions", angle = 90, size = 15, vjust = 0.5),
  strmaggie10_combined,
  ncol = 2,
  rel_widths = c(0.05, 1)  # 5% for label, 95% for plot
)

# Save the plot
ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/strmaggie10.pdf", plot = strmaggie10_labeled, width = 8, height = 10)
```

