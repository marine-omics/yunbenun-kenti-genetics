args <- commandArgs(trailingOnly = TRUE)


n <- strtoi(args[1])

library(tidyverse)
library(dartR.base)
library(dartR.popgen)
library(hierfstat)


data <- read_rds("hfstat_nokin.rds")
data_matched <- read_rds("hfstat.rds")

no.so.hfstat <- data[[2]]
adj.mi.hfstat <- data[[3]]

# Rarefaction depth for allelic.richness is fixed here rather than left at its
# default. Two reasons:
#
#  1. Shuffling population labels redistributes missing data between
#     populations, so the default depth (2 * min(ind.count)) drifts from one
#     permutation to the next and the null becomes incomparable to the
#     observed value.
#
#  2. Taking the minimum across this dataset and the all-individuals dataset
#     puts both permutation runs on the same depth, so their p-values are
#     directly comparable. Dropping individuals can only lower the per-locus
#     counts, so this minimum always comes from the kin-removed data (34 and 72
#     alleles for the no/so/ADJ and MI/ADJ groupings respectively).
common_min_n <- function(x, y){
  min(2 * min(ind.count(x), na.rm = TRUE), 2 * min(ind.count(y), na.rm = TRUE))
}

shuffle_hfstat <- function(x, min.n){
  x[,1] <- sample(x[,1],length(x[,1]),replace = FALSE)
  ar <- allelic.richness(x, min.n = min.n)
  bs <- basic.stats(x)
  list("Ar" = ar,"basic_stats" = bs)
}

# Seed from the array task id so individual permutations are reproducible. The
# offset differs from run_permute_all.R so the two nulls stay independent.
set.seed(7919 + n)

no.so.shuff <- shuffle_hfstat(no.so.hfstat, common_min_n(no.so.hfstat, data_matched[[2]]))

adj.mi.shuff <- shuffle_hfstat(adj.mi.hfstat, common_min_n(adj.mi.hfstat, data_matched[[3]]))

out <- list("no.so" = no.so.shuff, "adj.mi" = adj.mi.shuff)

write_rds(out,file = paste("perms_nokin/perm.",n,".rds",sep="",collapse=""))
