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

# 3 Relatedness

## Load Libraries

```{r}
library(ggpubr)       # to use ggarrange
library(tidyverse)
```

## Replicates removed

Close relatives (with a percentage of 99% similarity in their genotypes) were identified using gl.report.replicates in script "Population Structure" and were removed from this dataset.
This dataset was called ak.filtered.nr (nr = no replicates).

## Remove migrants

After running population structure analyses, migrants were identified and removed for the subsequent analysis.

```{r}
# Remove 19 migrants with a proportion >5%
ak.filtered.nr.nm <- gl.drop.ind(ak.filtered.nr, ind.list=c("HB_1", "HB_2", "HB_27/11/21_11", "HB_27/11/21_14", "HB_27/11/21_20", "HB_27/11/21_21", "WB_03/04/2022_8", "WB_03/04/2022_25", "WB_03/04/2022_37", "HFB_4", "HFB_6", "HFB_9", "HFB_10", "HFB_22", "HaR_06/05/2022_7", "HaR_06/05/2022_10", "HaR_06/05/2022_12", "HaR_06/05/2022_24", "HaR_06/05/2022_43"))
```

```{r}
# Filter for monomorphs
ak.filtered.nr.nm <- gl.filter.monomorphs(ak.filtered.nr.nm)
```

```{r}
# Recalculate locus metrics
ak.filtered.nr.nm <- gl.recalc.metrics(ak.filtered.nr.nm)
```

```{r}
nInd(ak.filtered.nr.nm)
nLoc(ak.filtered.nr.nm)
```

\[1\] 422 \[1\] 5265

```{r}
gl.smearplot(ak.filtered.nr.nm)
```

```{r}
save(ak.filtered.nr.nm, file="ak.filtered.nr.nm.rdata")
```

## Remove missing loci

A DArT dataset will not have individuals for which the calls are scored all as missing (NA) across all loci, but such individuals may sneak in to the dataset when loci are deleted. Retaining individual or loci with all NAs can cause issues for several functions.

```{r}
# Filter for missing loci
ak.rel <-gl.filter.allna(ak.filtered.nr.nm)
```

Starting gl.filter.allna Processing genlight object with SNP data Identifying and removing loci and individuals scored all missing (NA) Deleting loci that are scored as all missing (NA) *Zero loci that are missing (NA) across all individuals* Deleting individuals that are scored as all missing (NA) *Zero individuals that are missing (NA) across all loci* Completed: gl.filter.allna

## Inspect filtered file

```{r}
nInd(ak.rel)
nLoc(ak.rel)
```

```{r}
gl.smearplot(ak.rel)
```

```{r}
table(pop(ak.rel))
```

```{r}
save(ak.rel, file="ak.rel.rdata")
```

## Relatedness in vcftools

```{r}
# Convert dart file into vcf file
gl2vcf(ak.rel, plink.bin.path = "C:/Users/sandr", outfile = "ak.rel", outpath = "C:/Users/sandr")
```

Run relatedness analysis with vcftools using this code:

vcftools --vcf ak.rel.nm.vcf --relatedness2 --out relatedness_ak.rel.nm

### Plot Magnetic Island vs. Adjacent

```{r}
# Load csv file that contains vcftools output
relatedness_ak.rel <- read.csv('C:/Users/sandr/Documents/Chapter 1 - Github Script/relatedness_ak.rel.nm.csv')
```

```{r}
head(relatedness_ak.rel)
```

#### Boxplot Magnetic Island

```{r}
# Select only entries for the Magnetic Island populations
entries_mag = c("HB", "GB", "PB", "MR", "WB", "HFB", "MB")

mag.relate.box <- relatedness_ak.rel %>%
  filter(sapply(INDV1, function(x) strsplit(x, "_")[[1]][1]) %in% entries_mag) %>%
  filter(sapply(INDV2, function(x) strsplit(x, "_")[[1]][1]) %in% entries_mag)
```

```{r}
# Create a boxplot
boxplot_mag = ggplot(data = mag.relate.box, aes(y = RELATEDNESS_PHI)) + 
            geom_boxplot(outlier.size = 1) +
            theme_minimal()
boxplot_mag
```

#### Boxplot Adjacent

```{r}
# Select only entries for the adjacent reef populations
entries_adj = c("SO", "LPB", "EP", "WP", "KR", "JB", "HaR", "BRR")

adj.relate.box <- relatedness_ak.rel %>%
  filter(sapply(INDV1, function(x) strsplit(x, "_")[[1]][1]) %in% entries_adj) %>%
  filter(sapply(INDV2, function(x) strsplit(x, "_")[[1]][1]) %in% entries_adj)
```

```{r}
# Create a boxplot
boxplot_adj = ggplot(data = adj.relate.box, aes(y = RELATEDNESS_PHI)) + 
                geom_boxplot(outlier.size = 1) +
                theme_minimal()
boxplot_adj
```

```{r}
# Combine both boxplots into one plot
rel.box.adj.mag <- ggarrange(boxplot_mag, boxplot_adj, 
          ncol = 2, nrow = 1,
          widths = c(1,1),                # first number stand for the first column and the second for the second column
          heights = c(1,1))

rel.box.adj.mag

# Save the plot
ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/rel.box.adj.mag.pdf", plot = rel.box.adj.mag, width = 8, height = 10)
```

#### ANOVA

```{r}
# Calculate the means
mean(mag.relate.box$RELATEDNESS_PHI)
mean(adj.relate.box$RELATEDNESS_PHI)
```

```{r}
# Combine both dataframes into one and add the location
phi.adj.vs.mag <- rbind(
  mutate(adj.relate.box, value = RELATEDNESS_PHI, group = "Adjacent"),
  mutate(mag.relate.box, value = RELATEDNESS_PHI, group = "Magnetic"))
```

```{r}
# Run the ANOVA analysis
anova_phi <- aov(value ~ group, data = phi.adj.vs.mag)
```

```{r}
summary(anova_phi)
```

# 4 Genetic Diversity

Genetic diversity analyses are sensitive to low read depth and losing rare alleles.

Genetic diversity analyses are not very sensitive to - secondaries, because we don't worry about loci that are independently inherited or recombination - missing data and singletons so, no maf filter or imputing

## Load Libraries

```{r}
library(dartRverse)
library(dartR)   # to run gl.basic.stats
library(dartR.base)
library(hierfstat)        # allelic richness
library(ggplot2)
library(gplots)
library(directlabels)     # to display PCA's
library(vegan)            # to run a PERMANOVA
#library(boot)
library(plotly)         # interactive PCoA
library(ggthemes)       # to use theme_tufte()
library(ggpubr)       # to use ggarrange
library(cowplot)      # to use ggdraw
```

## Replicates removed

Close relatives (with a percentage of 99% similarity in their genotypes) were identified using gl.report.replicates in script "Population Structure" and were removed from this dataset.
This dataset was called ak.filtered.nr (nr = no replicates).

## Migrants removed 

Migrants were identified via population structure analysis in script "Population Structure" and were removed from this dataset in the script "Relatedness".
This dataset was called ak.filtered.nr.nm (nr = no replicates, nm = no migrants).

## Remove missing loci

A DArT dataset will not have individuals for which the calls are scored all as missing (NA) across all loci, but such individuals may sneak in to the dataset when loci are deleted. Retaining individual or loci with all NAs can cause issues for several functions

```{r}
# Filter for missing loci
ak.gen <-gl.filter.allna(ak.filtered.nr.nm)
```
## Inspect filtered file

```{r}
nInd(ak.gen)
nLoc(ak.gen)
```

```{r}
gl.smearplot(ak.gen)
```

```{r}
table(pop(ak.gen))
```

```{r}
save(ak.gen, file="ak.gen.rdata")
```

## Subset dataset

### Magnetic Island

Subset the dataset ak.gen for downstream analysis.

```{r}
# Select only subpopulations from Magnetic Island
ak.gen.mag <- gl.keep.pop(ak.gen, pop.list = list("GB", "HB", "HFB", "MB", "MR", "PB", "WB"))
```

```{r}
# Filter for monomorphs
ak.gen.mag <- gl.filter.monomorphs(ak.gen.mag)
```

```{r}
# Recalculate locus metrics
ak.gen.mag <- gl.recalc.metrics(ak.gen.mag)
```

```{r}
# Save the dataset
save(ak.gen.mag, file="C:/Users/sandr/Documents/Chapter 1 - Github Script/ak.gen.mag.rdata")
```

```{r}
nInd(ak.gen.mag)
```

#### Magnetic Island - North and South

```{r}
# Select only subpopulations from north Magnetic Island 
ak.gen.no <- gl.keep.pop(ak.gen.mag, pop.list = list("HFB", "WB", "MB", "HB"))
# Select only subpopulations from south Magnetic Island
ak.gen.so <- gl.keep.pop(ak.gen.mag, pop.list = list("MR", "PB", "GB"))

# Merge subpopulations from north Magnetic Island into one population
ak.gen.no.so <- gl.merge.pop(ak.gen.mag, old=c("HB", "HFB", "MB", "WB"), new='no')
# Using this merged dataset, merge also subpopulations from south Magnetic Island into one population. Now it is possible to compare north and south with each other
ak.gen.no.so <- gl.merge.pop(ak.gen.no.so, old=c("GB", "MR", "PB"), new='so')
```

```{r}
# Filter all for monomorphs
ak.gen.no <- gl.filter.monomorphs(ak.gen.no)
ak.gen.so <- gl.filter.monomorphs(ak.gen.so)
ak.gen.no.so <- gl.filter.monomorphs(ak.gen.no.so)
```

```{r}
# Recalculate locus metrics
ak.gen.no <- gl.recalc.metrics(ak.gen.no)
ak.gen.so <- gl.recalc.metrics(ak.gen.so)
ak.gen.no.so <- gl.recalc.metrics(ak.gen.no.so)
```

```{r}
# Save the datasets
save(ak.gen.no, file="ak.gen.no")
save(ak.gen.so, file="ak.gen.so")
save(ak.gen.no.so, file="ak.gen.no.so")
```

```{r}
ak.gen.no$ind.names
ak.gen.so$ind.names
ak.gen.no.so$ind.names
```

#### ak.mi.pi.ms

This dataset merges all subpopulations from Magnetic Island into one big population called "mi" (Magnetic Island), all populations from the Palm Islands into "pi" (Palm Islands) and all mid-shelf populations into "ms" (mid-shelf).
This is used for the Fis calculations per group: Magnetic Island (mi), Palm Island (pi) and mid-shelf reefs (ms).

```{r}
ak.gen.mi <- gl.merge.pop(ak.gen, old=c("GB", "HB", "HFB", "MB", "MR", "PB", "WB"), new='MI')
ak.gen.mi.pi <- gl.merge.pop(ak.gen.mi, old=c("SO", "LPB", "EP", "WP", "HaR"), new='PI')
ak.gen.mi.pi.ms <- gl.merge.pop(ak.gen.mi.pi, old=c("BRR", "JB", "KR"), new='MS')
ak.gen.mi.adj <- gl.merge.pop(ak.gen.mi.pi.ms, old=c("PI", "MS"), new='ADJ')
```

```{r}
save(ak.gen.mi, file="ak.gen.mi.rdata")
save(ak.gen.mi.pi, file="ak.gen.mi.pi.rdata")
save(ak.gen.mi.pi.ms, file="ak.gen.mi.pi.ms.rdata")
save(ak.gen.mi.adj, file="ak.gen.mi.adj.rdata")
```

### Adjacent Reefs

```{r}
# Select only populations from adjacent reefs
ak.gen.adj <- gl.keep.pop(ak.gen, pop.list = list("BRR", "HaR", "JB", "KR", "SO", "LPB", "EP", "WP"))
```

```{r}
# Filter for monomorphs
ak.gen.adj <- gl.filter.monomorphs(ak.gen.adj)
```

```{r}
# Recalculate locus metrics
ak.gen.adj <- gl.recalc.metrics(ak.gen.adj)
```

```{r}
# Save the dataset
save(ak.gen.adj, file="C:/Users/sandr/Documents/Chapter 1 - Github Script/ak.gen.adj.rdata")
```

```{r}
nInd(ak.gen.adj)
```

## Basic diversity metrics

```{r}
# Calculate basic diversity metrics for each loci (He (= Hs), Ho, Fis etc.)
gl.basic.stats(ak.gen)    
```

## Pairwise Fst

```{r}
# Calculate pairwise Fst values between populations
nclusters <- min(4,parallel::detectCores())
fst.pop <- gl.fst.pop(ak.gen, nboots=100, percent=95, nclusters = nclusters)
```

```{r}
# Create a table from the output
knitr::kable((round(fst.pop$Fsts,3)))
```

In theory Fst should be between 0 and 1, but in practice Fst can become negative when the heterozygosity in the subpopulation is higher than in the total population; or when sample sizes are imbalanced.

### Heatmap of Fst values - Magnetic Island

```{r}
# Get names of populations
popNames(ak.gen.mag)
```

```{r}
# Reorder populations based on geographic location (north to south)
ak.gen.mag_reord <-gl.sort(ak.gen.mag, sort.by = "pop", order.by = c("HFB", "WB", "MB", "HB", "GB", "PB", "MR"))
```

```{r}
# Calculate Fst values for Magnetic Island
fst.mag_reord <- gl.fst.pop(ak.gen.mag_reord)
```

```{r}
# Display Fst values
fst.mag_reord
```

### Plot Heatmap

```{r}
# Convert dist object to data.frame
fst.mag_reord.matrix = as.matrix(fst.mag_reord)
```

```{r}
# Create a logical matrix of the same dimensions as fst.nm.mag_reord.matrix, with TRUE values in the upper triangle (excluding the diagonal) and FALSE elsewhere. The which(..., arr.ind = TRUE) function then returns the row and column indices of the TRUE values in this logical matrix. These indices are stored in ind.mag
ind.mag = which(lower.tri(fst.mag_reord.matrix), arr.ind = TRUE)
```

```{r}
# Construct a data frame with three columns:
# Site1: Contains the names of the columns of fst.nm.mag_reord.matrix corresponding to the second index (ind[,2]).
# Site2: Contains the names of the rows of fst.nm.mag_reord.matrix corresponding to the first index (ind[,1]).
# Fst: Contains the values from the fst.nm.mag_reord.matrix at the positions specified by ind, rounded to three decimal places.
fst.mag_reord.df = data.frame(Site1 = dimnames(fst.mag_reord.matrix)[[2]][ind.mag[,2]],
                    Site2 = dimnames(fst.mag_reord.matrix)[[1]][ind.mag[,1]],
                    Fst = fst.mag_reord.matrix[ ind.mag ] %>% round(digits = 4))
```

```{r}
# Convert minus values to zero
fst.mag_reord.df$Fst[fst.mag_reord.df$Fst < 0] = 0
```

```{r}
# Print data.frame summary
fst.mag_reord.df %>% str
```

```{r}
# Create an italic label for "Fst"
fst.label = expression(italic("F")[ST])
```

```{r}
# Reorder site 1 and 2
mag.geo.order1 <- c("HFB", "WB", "MB", "HB", "GB", "PB", "MR")
mag.geo.order2 <- c("MR", "PB", "GB", "HB", "MB", "WB" ,"HFB")
```

```{r}
# Convert sites to factors
fst.mag_reord.df$Site1 <- factor(fst.mag_reord.df$Site1, levels = mag.geo.order1)
fst.mag_reord.df$Site2 <- factor(fst.mag_reord.df$Site2, levels = mag.geo.order2)
```

```{r}
# Plot heatmap
fst.mag_reord_plot<- ggplot(data = fst.mag_reord.df, aes(x = Site1, y = Site2, fill = Fst)) +
  geom_raster(aes(fill = Fst)) +
  scale_fill_gradient(low = "grey", high = "blue", na.value = NA) +
  geom_text(aes(label = Fst), color="white", size = 3) +
  theme_tufte()

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/fst.mag_reord_plot.png", width=800, height=400)

ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/fst_mag_plot.pdf", plot = fst.mag_reord_plot, width = 8, height = 10)

fst.mag_reord_plot
```

### Heatmap of Fst values - Adjacent Reefs

```{r}
# Get names of populations
popNames(ak.gen.adj)
```

```{r}
# Reorder populations based on geographic location (inshore to mid-shelf)
ak.gen.adj_reord <-gl.sort(ak.gen.adj, sort.by = "pop", order.by = c("HaR", "SO", "LPB", "WP", "EP", "BRR", "JB",  "KR"))
```

```{r}
# Calculate Fst values for adjacent reefs
fst.adj_reord <- gl.fst.pop(ak.gen.adj_reord)
```

```{r}
# Display Fst values
fst.adj_reord
```

### Plot Heatmap

Convert dist object to data.frame

```{r}
# Convert dist object to data.frame
fst.adj_reord.matrix = as.matrix(fst.adj_reord)
```

```{r}
# Create a logical matrix of the same dimensions as fst.nm.mag_reord.matrix, with TRUE values in the upper triangle (excluding the diagonal) and FALSE elsewhere. The which(..., arr.ind = TRUE) function then returns the row and column indices of the TRUE values in this logical matrix. These indices are stored in ind
ind.adj = which(lower.tri(fst.adj_reord.matrix), arr.ind = TRUE)
```

```{r}
# Construct a data frame with three columns:
# Site1: Contains the names of the columns of fst.nm.mag_reord.matrix corresponding to the second index (ind[,2]).
# Site2: Contains the names of the rows of fst.nm.mag_reord.matrix corresponding to the first index (ind[,1]).
# Fst: Contains the values from the fst.nm.mag_reord.matrix at the positions specified by ind, rounded to three decimal places.
fst.adj_reord.df = data.frame(Site1 = dimnames(fst.adj_reord.matrix)[[2]][ind.adj[,2]],
                    Site2 = dimnames(fst.adj_reord.matrix)[[1]][ind.adj[,1]],
                    Fst = fst.adj_reord.matrix[ ind.adj ] %>% round(digits = 4))
```

```{r}
# Convert minus values to zero
fst.adj_reord.df$Fst[fst.adj_reord.df$Fst < 0] = 0
```

```{r}
# Print data.frame summary
fst.adj_reord.df %>% str
```

```{r}
# Reorder site 1 and 2
adj.geo.order1 <- c("HaR", "SO", "LPB", "WP", "EP", "BRR", "JB", "KR")
adj.geo.order2 <- c("KR", "JB", "BRR", "EP", "WP", "LPB", "SO", "HaR")
```

```{r}
# Convert sites to factors
fst.adj_reord.df$Site1 <- factor(fst.adj_reord.df$Site1, levels = adj.geo.order1)
fst.adj_reord.df$Site2 <- factor(fst.adj_reord.df$Site2, levels = adj.geo.order2)
```

```{r}
# Plot heatmap
fst.adj_reord_plot<- ggplot(data = fst.adj_reord.df, aes(x = Site1, y = Site2, fill = Fst)) +
  geom_raster(aes(fill = Fst)) +
  scale_fill_gradient(low = "grey", high = "blue", na.value = NA) +
  geom_text(aes(label = Fst), color="white", size = 3) +
  theme_tufte()

png(file="C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/fst.adj_reord_plot.png", width=800, height=400)    #this code must sit within the chunk to work

ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/fst_adj_plot.pdf", plot = fst.adj_reord_plot, width = 8, height = 10)

fst.adj_reord_plot
```

## Boxplot Fst comparisons

```{r}
# Boxplot Magnetic Island
ggplot(fst.mag_reord.df, aes(y = Fst)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 0.012), breaks = seq(0, 0.012, by = 0.005)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
```

```{r}
# Boxplot adjacent reefs
ggplot(fst.adj_reord.df, aes(y = Fst)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 0.012), breaks = seq(0, 0.012, by = 0.005)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
```

##### Combined Boxplot Magnetic Island and adjacent reefs

```{r}
# First, we need to combine the data frames and add a group identifier
fst.mag_reord.df$Group <- "Magnetic Island"
fst.adj_reord.df$Group <- "Adjacent Reefs"
combined_mag_adj_df <- rbind(fst.mag_reord.df, fst.adj_reord.df)

# Create a single plot with two boxplots
boxplot_mag_adj <- ggplot(combined_mag_adj_df, aes(x = Group, y = Fst)) +
  geom_boxplot() +
  theme_classic() +
  labs(y = expression(F[st])) +
  theme(
    axis.title.x = element_blank())

ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/boxplot_mag_adj.pdf", plot = boxplot_mag_adj, width = 8, height = 10)
ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/boxplot_mag_adj.jpeg", width = 15, height = 10, units = "cm")

boxplot_mag_adj
```

#### Boxplot Fst North

```{r}
# Get names of populations
popNames(ak.gen.no)
```

```{r}
# Reorder populations based on geographic location (west to east)
ak.gen.no_reord <-gl.sort(ak.gen.no, sort.by = "pop", order.by = c("HFB", "WB", "MB", "HB"))
```

```{r}
# Calculate Fst values for northern subpopulations
fst.gen.no_reord <- gl.fst.pop(ak.gen.no_reord)
```

```{r}
# Display Fst values
fst.gen.no_reord
```

```{r}
# Convert dist object to data.frame
fst.gen.no_reord.matrix = as.matrix(fst.gen.no_reord)
```

```{r}
# Create a logical matrix of the same dimensions as fst.nm.mag_reord.matrix, with TRUE values in the upper triangle (excluding the diagonal) and FALSE elsewhere. The which(..., arr.ind = TRUE) function then returns the row and column indices of the TRUE values in this logical matrix. These indices are stored in ind
ind.mag.no = which(lower.tri(fst.gen.no_reord.matrix), arr.ind = TRUE)
```

```{r}
# Construct a data frame with three columns:
# Site1: Contains the names of the columns of fst.nm.mag_reord.matrix corresponding to the second index (ind[,2]).
# Site2: Contains the names of the rows of fst.nm.mag_reord.matrix corresponding to the first index (ind[,1]).
# Fst: Contains the values from the fst.nm.mag_reord.matrix at the positions specified by ind, rounded to three decimal places.
fst.gen.no_reord.df = data.frame(Site1 = dimnames(fst.gen.no_reord.matrix)[[2]][ind.mag.no[,2]],
                    Site2 = dimnames(fst.gen.no_reord.matrix)[[1]][ind.mag.no[,1]],
                    Fst = fst.gen.no_reord.matrix[ ind.mag.no ] %>% round(digits = 4))
```

```{r}
# Plot northern bays
north_boxplot <- ggplot(fst.gen.no_reord.df, aes(y = Fst)) +
  geom_boxplot() +
#  scale_y_continuous(limits = c(0, 0.012), breaks = seq(0, 0.012, by = 0.005)) +
  theme_classic() +
  theme(
        #axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
north_boxplot
```

#### Boxplot Fst South

```{r}
# Get names of populations
popNames(ak.gen.so)
```

```{r}
# Reorder populations based on geographic location (north-south)
ak.gen.so_reord <-gl.sort(ak.gen.so, sort.by = "pop", order.by = c("GB", "PB", "MR"))
```

```{r}
# Calculate Fst values for southern subpopulations
fst.gen.so_reord <- gl.fst.pop(ak.gen.so_reord)
```

```{r}
# Display Fst values
fst.gen.so_reord
```

```{r}
# Convert dist object to data.frame
fst.gen.so_reord.matrix = as.matrix(fst.gen.so_reord)
```

```{r}
# Create a logical matrix of the same dimensions as fst.nm.mag_reord.matrix, with TRUE values in the upper triangle (excluding the diagonal) and FALSE elsewhere. The which(..., arr.ind = TRUE) function then returns the row and column indices of the TRUE values in this logical matrix. These indices are stored in ind
ind.mag.so = which(lower.tri(fst.gen.so_reord.matrix), arr.ind = TRUE)
```

```{r}
# Construct a data frame with three columns:
# Site1: Contains the names of the columns of fst.nm.mag_reord.matrix corresponding to the second index (ind[,2]).
# Site2: Contains the names of the rows of fst.nm.mag_reord.matrix corresponding to the first index (ind[,1]).
# Fst: Contains the values from the fst.nm.mag_reord.matrix at the positions specified by ind, rounded to three decimal places.
fst.gen.so_reord.df = data.frame(Site1 = dimnames(fst.gen.so_reord.matrix)[[2]][ind.mag.so[,2]],
                    Site2 = dimnames(fst.gen.so_reord.matrix)[[1]][ind.mag.so[,1]],
                    Fst = fst.gen.so_reord.matrix[ ind.mag.so ] %>% round(digits = 4))
```

```{r}
# Plot southern bays
south_boxplot <- ggplot(fst.gen.so_reord.df, aes(y = Fst)) +
  geom_boxplot() +
#  scale_y_continuous(limits = c(0, 0.012), breaks = seq(0, 0.012, by = 0.005)) +
  theme_classic() +
  theme(#axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
south_boxplot
```

### Combined Boxplot North and South

```{r}
# Combine the data frames and add a group identifier
fst.gen.no_reord.df$Group <- "North"
fst.gen.so_reord.df$Group <- "South"
combined_no_so_df <- rbind(fst.gen.no_reord.df, fst.gen.so_reord.df)

# Create a single plot with two boxplots
boxplot_no_so <- ggplot(combined_no_so_df, aes(x = Group, y = Fst)) +
  geom_boxplot() +
  theme_classic() +
  labs(y = expression(F[st])) +
  theme(
    axis.title.x = element_blank()
  )

# Save the plot
ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/boxplot_no_so.pdf", plot = boxplot_no_so, width = 4, height = 8)

boxplot_no_so
```

## Significance Testing Fst

```{r}
# Compare between Magnetic Island Fst values and adjacent reefs
Fst_signif <- t.test(fst.mag_reord.df$Fst, fst.adj_reord.df$Fst, var.equal = FALSE)

# Print the result
print(Fst_signif)
```
Conclusion: Magnetic Island has a mean fst of 0.0052142857, and adjacent reefs have a mean fst of 0.0005071429.

### North vs South

```{r}
# Compare between northern Magnetic Island Fst values and southern
Fst_signif_no_vs_so <- t.test(fst.gen.no_reord.df$Fst, fst.gen.so_reord.df$Fst, var.equal = FALSE)

# Print the result
print(Fst_signif_no_vs_so)
```
Conclusion: Northern bays have mean fst of -0.0008166667, southern bays have higher mean fst of 0.0079000000.

## Heatmap Euclidean Distance

### pairwise individual

```{r}
# Calculate an Euclidean Distance Matrix for individuals
ak.gen.EDM <- gl.dist.ind(ak.gen)
```

```{r}
# Create and save the heatmap showing Euclidean Distance, PCA and neighbour joining tree
pdf("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/heatmap_EDM.pdf", width = 10, height = 8)

# Plot the heatmap
gl.plot.heatmap(ak.gen.EDM)

dev.off()
```

### pairwise population

```{r}
# Calculate an Euclidean Distance Matrix for populations
ak.gen.EDM.pop <- gl.dist.pop(ak.gen)
```

```{r}
# Display Euclidean Distance Matrix per population
ak.gen.EDM.pop
```

```{r}
# Create and save the heatmap showing Euclidean Distance, PCA and neighbour joining tree
pdf("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/heatmap_EDM_pop.pdf", width = 10, height = 8)

# Plot the heatmap
gl.plot.heatmap(ak.gen.EDM.pop)

dev.off()
```

## Allelic richness

```{r}
# Convert the genlight object to genind format
ak.gen.gi <- gl2gi(ak.gen)
ak.gen.mag.gi <- gl2gi(ak.gen.mag)   # to compare Magnetic Island to adjacent reefs
ak.gen.adj.gi <- gl2gi(ak.gen.adj)   # to compare Magnetic Island to adjacent reefs
```

```{r}
# Save genind object
save(ak.gen.gi, file="ak.gen.gi.rdata")
save(ak.gen.mag.gi, file="ak.gen.mag.gi.rdata")
save(ak.gen.adj.gi, file="ak.gen.adj.gi.rdata")
```

```{r}
# Convert genind object to hierfstat format
ak.gen.hfstat <- genind2hierfstat(ak.gen.gi)
ak.gen.mag.hfstat <- genind2hierfstat(ak.gen.mag.gi)
ak.gen.adj.hfstat <- genind2hierfstat(ak.gen.adj.gi)
```

```{r}
# Calculate allelic richness
ar <- allelic.richness(ak.gen.hfstat)
ar.mag <- allelic.richness(ak.gen.mag.hfstat)
ar.adj <- allelic.richness(ak.gen.adj.hfstat)
```

```{r}
# Show summary of Ar statistics
summary(ar$Ar)     # gives mean AR for each population
summary(ar.mag$Ar) # gives mean AR for each population
summary(ar.adj$Ar) # gives mean AR for each population

```

### Significance test Ar

#### t-test

```{r}
# Compare Fis values between Magnetic Island and adjacent reefs
Ar_signif <- t.test(ar.mag$Ar, ar.adj$Ar, var.equal = FALSE)

# Print the result
print(Ar_signif)
```

Conclusion:

A very small p-value (typically less than 0.05) suggests strong evidence against the null hypothesis, i.e. there is no difference between the groups. Hence *significant difference*
95 percent confidence interval does not include 0, hence significant.

## Heterozygosity/Fis

Since we found very high FIS values for Magnetic Island, which are likely caused by the Wahlund effect, we analyse the Magnetic Island as one population against adjacent reefs.

```{r}
# Create colour scheme
# Get order of population names
#levels(ak.gen.mi.adj@pop)
#  "ADJ" "MI"
# Magnetic Island = orange "orange"; Adjacent = blue "blue";
cols.mi.adj <- c("blue", "orange")
```

```{r}
# Report expected and observed heterozygosity, and inbreeding coefficient
het.mi.adj <- gl.report.heterozygosity(ak.gen.mi.adj, method="pop", error.bar = "SE", plot.colors.pop = cols.mi.adj)
```

If you want to know how accurately your sample mean represents the population mean, you'd look at the SE, which is is small. 

## Significance Testing Fis

```{r}
# Create a data frame with FIS values from Maggie
het_Fis_mag <- het %>%
  filter(pop %in% c("GB", "HB", "HFB", "MB", "MR", "PB", "WB"))

# Print the results
print(het_Fis_mag)
```

```{r}
# Create a data frame with FIS values from adjacent reefs
het_Fis_adj <- het %>%
  filter(pop %in% c("SO", "LPB", "EP", "WP", "HaR", "BRR", "JB", "KR"))

# Print the results
print(het_Fis_adj)
```

```{r}
# Run a t-test to compare between Magnetic Island Fis values, and all other adjacent reefs
Fis_signif <- t.test(het_Fis_mag$FIS, het_Fis_adj$FIS, var.equal = FALSE)

# Print the result
print(Fis_signif)
```
Conclusion:
*not significantly different*

## Significance Testing Ho

```{r}
# Run a t-test to compare between Magnetic Island Ho values, and all other adjacent reefs
Ho_signif <- t.test(het_Fis_mag$Ho, het_Fis_adj$Ho, var.equal = FALSE)

# Print the result
print(Ho_signif)
```
Conclusion:
*not significantly different*

## Significance Testing He

```{r}
# # Run a t-test to compare between Magnetic Island He values, and all other adjacent reefs
He_signif <- t.test(het_Fis_mag$He, het_Fis_adj$He, var.equal = FALSE)

# Print the result
print(He_signif)
```
Conclusion:
*not significantly different*

## Fixed allelic differences/Private alleles

This time I have removed any hybrids and migrants. This analysis is population size dependent.

```{r}
pa_ak.gen <-dartR::gl.report.pa(ak.gen)
#gl.report.pa(ak.gen)       # if you run this without naming it, you get a graph of connections
```

Table of private alleles and fixed differences returned

*gl.fixed.diff* generates a matrix of fixed differences between populations taken pairwise. Fixed differences occur between Maggie subpopulations and others, unless amalgamated pops.

```{r}
fd <-gl.fixed.diff(ak.gen)
fd.mi.pi.ms <-gl.fixed.diff(ak.gen.mi.pi.ms)
fd.mi.pi.ms
```

### amalgamated pops

*Reran this for Maggie vs Adjacent, because some N are <10!*: no fixed differences detected, only private alleles!

```{r}
pa_ak.gen.mi <- gl.report.pa(ak.gen.mi)
```

When MI is amalgamated, there are fixed diff between Maggie population and adjacent populations!!! Plus many private alleles in all populations. But I have removed migrants and hybrids for this!

```{r}
pa_ak.gen.mi.pi.ms <- gl.report.pa(ak.gen.mi.pi.ms)
```

When all populations are amalgamated, there are fixed diff between Maggie and mid-shelf populations!!! Only many private alleles at Palm Islands! But I have removed migrants and hybrids for this!

### Heatmaps fixed diff

#### All pops non amalgamated

```{r}
fd.heatmap <- as.matrix(gl.fixed.diff(ak.gen)$fd)
```

```{r}
gl.plot.heatmap(fd.heatmap, margins = c(5, 10))
```

MB and WP show fixed diff!

#### Magnetic Island non amalgamated

```{r}
# You can create a heatmap of your fixed differences in one population (Magnetic Island only?) by first creating a matrix
fd.mag.heatmap <- as.matrix(gl.fixed.diff(ak.gen.mag)$fd)
```

Starting gl.fixed.diff Processing genlight object with SNP data Comparing populations for absolute fixed differences Monomorphic loci removed Populations, aggregations and sample sizes GB HB HFB MB MR PB WB 28 22 14 6 13 20 17 *Warning: Fixed differences can arise through sampling error if sample sizes are small* *Some sample sizes are small (N \< 10, minimum in dataset = 6 )* *Recommend manually amalgamating populations or setting test=TRUE to allow evaluation of statistical significance*

```{r}
gl.plot.heatmap(fd.mag.heatmap, margins = c(10, 10))
```

Error in (function (side, at = NULL, labels = TRUE, tick = TRUE, line = NA, : no locations are finite

#### Magnetic Island amalgamated

```{r}
# Create a heatmap of your fixed differences between Magnetic Island only and all others by first creating a matrix
fd.mi.heatmap <- as.matrix(gl.fixed.diff(ak.gen.mi)$fd)
```

Starting gl.fixed.diff Processing genlight object with SNP data Comparing populations for absolute fixed differences Monomorphic loci removed Populations, aggregations and sample sizes BRR EP MI HaR JB KR LPB SO WP 20 104 120 23 15 15 34 83 8 *Warning: Fixed differences can arise through sampling error if sample sizes are small* *Some sample sizes are small (N \< 10, minimum in dataset = 6 )* *Recommend manually amalgamating populations or setting test=TRUE to allow evaluation of statistical significance*

```{r}
gl.plot.heatmap(fd.mi.heatmap, margins = c(5, 10))
```

#### MI, PI and MS amalgamated

```{r}
# create a heatmap of fixed diff of amalgamated populations for MI, PI and MS
fd.mi.pi.ms.heatmap <- as.matrix(gl.fixed.diff(ak.gen.mi.pi.ms)$fd)
```

Starting gl.fixed.diff Processing genlight object with SNP data Comparing populations for absolute fixed differences Monomorphic loci removed Populations, aggregations and sample sizes MS PI MI 50 252 120 Comparing populations pairwise -- this may take time. Please be patient Completed: gl.fixed.diff

```{r}
gl.plot.heatmap(fd.mi.pi.ms.heatmap, margins = c(10, 10))
```

fixed diff between MI and MS!!!!

## Rare alleles

Rare alleles are typically defined as those with a minor allele frequency (MAF) below a certain threshold, often less than 1% or 5%.

```{r}
gl.report.maf(ak.gen)
```

```{r}
gl.report.maf(ak.gen.mi.pi.ms)
```

## Neighbour joining tree

```{r}
# Create and save the neighbour joining tree
pdf("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/nj.pdf", width = 8, height = 4)
nj <- gl.tree.nj(ak.gen, type="phylogram")
nj
dev.off()
```

```{r}



```

## PCoA

I think keeping migrants for this analysis is more plausible, since removing them shows better structure, but doesn't reflect the "true" situation.
PCA's are sensitive to missing data unlike PCoA's. So, I'm running a PCoA.

```{r}
# Run the PCoA calculations
pcoa <- gl.pcoa(ak.gen)
```

### Assign colours to populations

```{r}
# Get order of population names
# levels(ak@pop)
# "BRR" "EP"  "GB"  "HaR" "HB"  "HFB" "JB"  "KR"  "LPB" "MB"  "MR"  "PB"  "SO"  "WB"  "WP"

# Assign colours to populations in that order:
# Mid-shelf = yellow "#FDD835"; Palm Island = green "#74c476"; Magnetic Island = orange "orange";
cols <- c("#FDD835", "#74c476", "orange", "#74c476", "orange", "orange", "#FDD835", "#FDD835", "#74c476", "orange", "orange", "orange", "#74c476",  "orange", "#74c476")
```

```{r}
# Get order of population names
#levels(ak.gen@pop)
# "BRR" "EP"  "GB"  "HaR" "HB"  "HFB" "JB"  "KR"  "LPB" "MB"  "MR"  "PB"  "SO"  "WB"  "WP"

# Assign colours to populations in that order: 
# Mid-shelf and Palm Island = blue "blue"; Magnetic Island = orange "orange";
cols.gen <- c("blue", "blue", "orange", "blue", "orange", "orange", "blue", "blue", "blue", "orange", "orange", "orange", "blue",  "orange", "blue")
```

```{r}
# Get order of population names
#levels(ak.gen.mag@pop)
# "GB"  "HB"  "HFB"  "MB"  "MR"  "PB"   "WB"

# Assign colours to populations in that order: 
# south = blue "blue"; north = red "red";
cols.gen.mag <- c("darkgreen", "red", "red",  "red", "darkgreen", "darkgreen", "red")
```

### Plot PCA's

```{r}
# Plot the first two dimensions of the PCA
pcoa <- gl.pcoa.plot(pcoa, ak.gen, pt.colors = cols.gen, pop.labels="legend")
```

```{r}
# Save the plot
ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/pcoa.pdf", plot = pcoa, width = 10, height = 6)
```

```{r}
# Interactive PCoA
gl.pcoa.plot(pcoa, ak.gen, pop.labels="interactive", pt.colors = cols)
ggplotly()
```

```{r}
# 3D PCoA
gl.pcoa.plot(pcoa, ak.gen, pop.labels="pop", pt.colors = cols, xaxis=1, yaxis=2, zaxis=3)
```

## PCoA Magnetic Island

```{r}
# Calculate the PCA for Magnetic Island
pcoa.mag <- gl.pcoa(ak.gen.mag)
```

```{r}
# Compare northern vs southern subpopulations
pcoa_maggie <- gl.pcoa.plot(pcoa.mag, ak.gen.mag, pt.colors = cols.gen.mag, pop.labels="legend")
```

```{r}
# Save the plot
ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/pcoa_maggie.pdf", plot = pcoa_maggie, width = 10, height = 6)
```

```{r}
# Interactive PcoA
gl.pcoa.plot(pcoa.mag, ak.gen.mag, pop.labels="interactive", xaxis=1, yaxis=2)
ggplotly()
```

```{r}
# 3D PCoA
gl.pcoa.plot(pcoa.mag, ak.gen.mag, pop.labels="pop", xaxis=1, yaxis=2, zaxis=3)
```

### PERMANOVA Magnetic Island north/south

```{r}
# Get the distance matrix
dist_matrix_mag <- dist(pcoa.mag$scores)  # uses Euclidean distance
```

```{r}
# run a PERMANOVA test
permanova_result <- adonis2(dist_matrix_mag ~ ak.gen.mag@pop, permutations = 999)
print(permanova_result)
```
Significant difference between subpopulations, but how does it compare between north and south?

```{r}
# First, create a named vector of North/South classification
# Replace these with your actual population names
north_pops <- c("WB", "HFB", "MB", "HB")
south_pops <- c("GB", "PB", "MR")
```

```{r}
# Assign region to each individual based on their population
pop_names <- pop(ak.gen.mag)
region <- ifelse(pop_names %in% north_pops, "North", "South")
```

```{r}
# Add region info to your genlight object
ak.gen.mag@other$region <- factor(region)
```

```{r}
# Run PERMANOVA with the new region factor
adonis2(dist_matrix_mag ~ ak.gen.mag@other$region, permutations = 999)
```
statistically significant structure between North and South.

## Ne

### Ne all subpopulations

```{r}
gl.LDNe(ak.gen, outfile = "LD_ak.gen.txt")
```

Conclusion: Since there are some infinite numbers (HJB and MB), calculate Ne for Magnetic Island, Palm Island and mid-shelf reefs!

### Ne Magnetic Island, Palm Island and mid-shelf reefs

```{r}
# Calculate the population size at Magnetic Island, Palm Island and mid-shelf reefs
gl.LDNe(ak.gen.mi.pi.ms, outfile = "LD_ak.gen.mi.pi.ms.txt")
```

#### Plot Ne

```{r}
ne <- read.csv("Ne results.csv")
ne
```

```{r}
# Keep only rows where condition equals "No_S"
ne_noS <- ne %>% 
  filter(condition == "No_S")
```

```{r}
# plot shows correlation between all populations
ne.plot <- ggplot(ne_noS, aes(x=Population, y=value, fill=Population)) +
            geom_bar(position="dodge", stat="identity") +
            theme_classic() +
           scale_fill_manual(values=c('orange', 'blue','blue'))
     #       labs(x = "Geographic distance (km)") +
      #      labs(y = "Genetic distance (Fst)")

# Save the plot
ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/ne.plot.pdf", plot = ne.plot, width = 6, height = 4)

ne.plot
```

```{r}
# Test for significance
anova_ne <- aov(value ~ Population, data = ne)
summary(anova_ne)
```
            Df Sum Sq Mean Sq F value Pr(>F)  
Population   2 584566  292283   17.51 0.0222 *
Residuals    3  50085   16695                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```{r}
# Run a post-hoc test to see, which group differ significantly
TukeyHSD(anova_ne)
```
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = value ~ Population, data = ne)

$Population
                                diff        lwr       upr     p adj
Midshelf reefs-Magnetic Island 504.8  -35.13531 1044.7353 0.0594214
Palm Islands-Magnetic Island   749.7  209.76469 1289.6353 0.0206869
Palm Islands-Midshelf reefs    244.9 -295.03531  784.8353 0.2841806

### Ne Magnetic Island north/south
```{r}
# Calculate the population size at northern and southern bays at Magnetic Island
gl.LDNe(ak.gen.no.so, outfile = "LD_ak.gen.no.so.txt")
```

#### Plot Ne

```{r}
ne.no.so <- read.csv("Ne.no.so_results.csv")
ne.no.so
```

```{r}
# Keep only rows where condition equals "No_S"
ne.no.so_noS <- ne.no.so %>% 
  filter(condition == "No_S")
```

```{r}
# plot shows correlation between all populations
ne.no.so.plot <- ggplot(ne.no.so_noS, aes(x=Population, y=value, fill=Population)) +
            geom_bar(position="dodge", stat="identity") +
            theme_classic() +
           scale_fill_manual(values=c('darkgreen','red'))
     #       labs(x = "Geographic distance (km)") +
      #      labs(y = "Genetic distance (Fst)")

# Save the plot
ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/ne.no.so.plot.pdf", plot = ne.no.so.plot, width = 6, height = 4)

ne.no.so.plot
```

```{r}
t.test(value ~ Population, data = ne.no.so, var.equal = TRUE)
```

	Two Sample t-test

data:  value by Population
t = 6.7352, df = 2, p-value = 0.02134
alternative hypothesis: true difference in means between group North and group South is not equal to 0
95 percent confidence interval:
  256.9929 1166.1071
sample estimates:
mean in group North mean in group South 
             791.65               80.10

# 5 IBD

# Load libraries

```{r}
library(ggplot2)
library(dartRverse)
```

# Isolation by distance

Isolation by distance analysis based on a mantel test. Euclidean and genetic distance matrix are calculated (currently only pairwise Fst between population is implemented by default).

```{r}
ibd <- gl.ibd(ak.gen)
```

```{r}
# Save the plot
ggsave("C:/Users/sandr/Documents/Chapter 1 - Github Script/Figures/ibd.pdf", plot = ibd, width = 8, height = 4)
```


# IBD all populations

Distances between all populations via sea were measured with Google Earth Pro and pairwise Fst values entered into a csv file called "IBD_ak.gen".

```{r}
ibd<- read.csv("IBD_ak.gen.csv")
```

```{r}
# plot shows correlation between all populations
ibd.plot <- ggplot(ibd, aes(x=Dist, y=Fst)) +     # fill = Pop
            geom_point(aes(colour = Pop), show.legend = TRUE) +
            scale_color_manual(values=c('blue','brown', 'orange')) +
            geom_smooth(method = "lm", se = FALSE, col = "black", linewidth = 0.6) +
            labs(x = "Geographic distance (km)") +
            labs(y = "Genetic distance (Fst)") +
            theme_classic()

png(file="C:/Users/sandr/OneDrive - James Cook University/DarT/R/Output files/ibdplot_feb25.png",
width=600, height=400)

ibd.plot
```

```{r}
summary(lm(ibd$Fst ~ ibd$Dist))
```

# IBD inter-population

```{r}
ibd_inter<- read.csv("IBD_ak.gen.csv")
```

```{r}
ibd_inter <- filter(ibd_inter, Pop == "Interpopulation")

print(ibd_inter)
```

```{r}
# plot shows correlation between populations between Maggie and adjacent only (Inter-population), and distinguishes between northern and southern bays at Maggie 
ibd_inter.plot <- ggplot(ibd_inter, aes(x=Dist, y=Fst)) +     # fill = Pop
            geom_point(aes(colour = Orientation), show.legend = TRUE) +
            scale_color_manual(values=c('blue', 'orange')) +
            geom_smooth(method = "lm", se = FALSE, col = "black", linewidth = 0.6) +
            labs(x = "Geographic distance (km)") +
            labs(y = expression("Genetic distance (F"["ST"]*")")) +
            theme_classic()

png(file="C:/Users/sandr/OneDrive - James Cook University/DarT/R/Output files/ibd_inter_feb25.png",
width=600, height=400)

ibd_inter.plot
```

```{r}
# plot shows correlation between populations between Maggie and adjacent only (Inter-population), and distinguishes between northern and southern bays at Maggie, and which bays at Maggie correlate to adjacent populations
ibd_inter_bays.plot <- ggplot(ibd_inter, aes(x=Dist, y=Fst)) +
  geom_point(aes(colour = Bay, shape = Orientation), show.legend = TRUE) +  # Add shape aesthetic
  scale_color_manual(values = c('blue', 'orange', "red", "green", "yellow", "black", "brown")) +
  scale_shape_manual(values = c('North' = 22, 'South' = 16)) +  # Adjust shape values here
  geom_smooth(method = "lm", se = FALSE, col = "black", linewidth = 0.6) +
  labs(x = "Geographic distance (km)") +
  labs(y = expression("Genetic distance (F"["ST"]*")")) +
  theme_classic()

png(file="C:/Users/sandr/OneDrive - James Cook University/DarT/R/Output files/ibd_inter_bays_feb25.png",
width=600, height=400)

ibd_inter_bays.plot
```

```{r}
summary(lm(ibd_inter$Fst ~ ibd_inter$Dist))
```


# IBD intra-population Magnetic Island

```{r}
ibd_intra<- read.csv("IBD_ak.gen.csv")
```

```{r}
ibd_intra_mag <- filter(ibd_intra, Pop == "Magnetic")

print(ibd_intra_mag)
```

```{r}
# plot shows correlation between populations at Magnetic Island
ibd_intra_mag.plot <- ggplot(ibd_intra_mag, aes(x=Dist, y=Fst)) +     # fill = Pop
            geom_point(aes()) +
            geom_smooth(method = "lm", se = FALSE, col = "black", linewidth = 0.6) +
            labs(x = "Geographic distance (km)") +
            labs(y = expression("Genetic distance (F"["ST"]*")")) +
            theme_classic()

png(file="C:/Users/sandr/OneDrive - James Cook University/DarT/R/Output files/ibd_intra_mag_feb25.png",
width=600, height=400)

ibd_intra_mag.plot
```

```{r}
# plot shows correlation between populations at Maggie (Magnetic population), and shows correlations between northern only and southern bays only at Maggie, and shows which northern bays at Maggie correlate to southern populations and vice versa (inter).
ibd_intra_mag_bays.plot <- ggplot(ibd_intra_mag, aes(x=Dist, y=Fst)) +
  geom_point(aes(colour = Orientation), show.legend = TRUE) +  # Add shape aesthetic
  scale_color_manual(values = c('darkgreen', "red", "black")) +
  #scale_shape_manual(values = c('North' = 22, 'South' = 16)) +  # Adjust shape values here
  geom_smooth(method = "lm", se = FALSE, col = "black", linewidth = 0.6) +
  labs(x = "Geographic distance (km)") +
  labs(y = expression("Genetic distance (F"["ST"]*")")) +
  theme_classic()

png(file="C:/Users/sandr/OneDrive - James Cook University/DarT/R/Output files/ibd_intra_mag_bays_feb25.png",
width=600, height=400, dpi=600)

ibd_intra_mag_bays.plot
```

```{r}
summary(lm(ibd_inter_mag$Fst ~ ibd_inter_mag$Dist))
```


# IBD intra-population Adjacent

```{r}
ibd_intra<- read.csv("IBD_ak.gen.csv")
```

```{r}
ibd_intra_adj <- filter(ibd_inter, Pop == "Adjacent")

print(ibd_intra_adj)
```

```{r}
# plot shows correlation between populations at adjacent populations (Adjacent population)
ibd_intra_adj.plot <- ggplot(ibd_intra_adj, aes(x=Dist, y=Fst)) + 
            geom_point(aes()) +
            geom_smooth(method = "lm", se = FALSE, col = "black", linewidth = 0.6) +
            labs(x = "Geographic distance (km)") +
            labs(y = expression("Genetic distance (F"["ST"]*")")) +
            theme_classic()

png(file="C:/Users/sandr/OneDrive - James Cook University/DarT/R/Output files/ibd_intra_adj_feb25.png",
width=600, height=400)

ibd_intra_adj.plot
```

# IBD intra-population Maggie and Adjacent

```{r}
ibd_intra<- read.csv("IBD_ak.gen.csv")
```

```{r}
ibd_intra_mag_adj <- filter(ibd_intra, Pop == "Adjacent" | Pop =="Magnetic")

print(ibd_intra_mag_adj)
```

```{r}
# plot shows correlation between populations at adjacent populations only (Adjacent population) between populations at Magnetic populations only (Magnetic population)
ibd_intra_mag_adj.plot <- ggplot(ibd_intra_mag_adj, aes(x=Dist, y=Fst)) +     # fill = Pop
            geom_point(aes(colour = Pop), show.legend = TRUE) +
            scale_color_manual(values=c('blue', 'orange')) +
            geom_smooth(method = "lm", se = FALSE, col = "black", linewidth = 0.6) +
            labs(x = "Geographic distance (km)") +
            labs(y = expression("Genetic distance (F"["ST"]*")")) +
            theme_classic()

png(file="C:/Users/sandr/OneDrive - James Cook University/DarT/R/Output files/ibd_intra_mag_adj_feb25.png",
width=600, height=400)

ibd_intra_mag_adj.plot
```

```{r}
summary(lm(ibd_intra_mag_adj$Fst ~ ibd_intra_mag_adj$Dist))
```


