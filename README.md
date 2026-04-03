# Connectivity and Diversity of Acropora kenti

The dataset analysed in this repository is from this publication (XXX)

Analysis is divided into the following steps

- [Data wrangling and filtering](01.data_wrangling_and_filtering.md)
- [Remove replicated and highly related individuals](02.remove_clones.md)
- [Population Structure](03.population_structure.md)
- [Genetic Diversity](04.genetic_diversity.md)
- [Estimating Effective Population Size](05.ne_estimates.md)
- [Long Runs of Homozygosity](06.LROH)


## Running the code

Rendered versions of all quarto files for this repository are available at the links above.  If you would like to render these yourself and potentially edit/run the code you will need additional data.  Download data and code as follows

```bash
git clone git@github.com:marine-omics/yunbenun-kenti-genetics.git
wget http://data.qld.edu.au/public/Q5999/marine-omics/yunbenun-kenti-genetics/data.tgz
tar -zxvf data.tgz
```

# Raw Data

Unfiltered genotype calls from DaRT sequencing are available in DaRT 2-row format

```bash
wget http://data.qld.edu.au/public/Q5999/marine-omics/yunbenun-kenti-genetics/raw_data/Report_DAc23-8597_5_moreOrders_SNP_2.csv
```

