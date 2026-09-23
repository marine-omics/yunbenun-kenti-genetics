args <- commandArgs(trailingOnly = TRUE)


n <- strtoi(args[1])

library(tidyverse)
library(dartR.base)
library(dartR.popgen)
library(hierfstat)


data <- read_rds("hfstat.rds")

no.so.hfstat <- data[[2]]
adj.mi.hfstat <- data[[3]]

# Rarefaction depth for allelic.richness is fixed rather than left at its
# default. Shuffling population labels redistributes missing data between
# populations, so the default depth (2 * min(ind.count)) drifts from one
# permutation to the next and the null becomes incomparable to the observed
# value. Computing it once from the unshuffled data keeps every permutation at
# the same depth as the observed statistic.
observed_min_n <- function(x){
  2 * min(ind.count(x), na.rm = TRUE)
}

shuffle_hfstat <- function(x, min.n){
  x[,1] <- sample(x[,1],length(x[,1]),replace = FALSE)
  ar <- allelic.richness(x, min.n = min.n)
  bs <- basic.stats(x)
  list("Ar" = ar,"basic_stats" = bs)
}


no.so.shuff <- shuffle_hfstat(no.so.hfstat, observed_min_n(no.so.hfstat))

adj.mi.shuff <- shuffle_hfstat(adj.mi.hfstat, observed_min_n(adj.mi.hfstat))

out <- list("no.so" = no.so.shuff, "adj.mi" = adj.mi.shuff)

write_rds(out,file = paste("perms/perm.",n,".rds",sep="",collapse=""))

