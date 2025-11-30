args <- commandArgs(trailingOnly = TRUE)


n <- strtoi(args[1])

library(tidyverse)
library(dartR.base)
library(dartR.popgen)
library(hierfstat)


data <- read_rds("hfstat.rds")

no.so.hfstat <- data[[2]]
adj.mi.hfstat <- data[[3]]

shuffle_hfstat <- function(x){
  x[,1] <- sample(x[,1],length(x[,1]),replace = FALSE)
  x
  ar <- allelic.richness(x)
  bs <- basic.stats(x)
  list("Ar" = ar,"basic_stats" = bs)
}


no.so.shuff <- shuffle_hfstat(no.so.hfstat)

adj.mi.shuff <- shuffle_hfstat(adj.mi.hfstat)

out <- list("no.so" = no.so.shuff, "adj.mi" = adj.mi.shuff)

write_rds(out,file = paste("perms/perm.",n,".rds",sep="",collapse=""))

