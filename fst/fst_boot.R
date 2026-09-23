library(tidyverse)
library(dartR.base)
library(dartR.popgen)

start.time <- Sys.time()
options(bspm.sudo = TRUE)

ak.pop.nm <- read_rds("ak.pop.nm.rds")

fst.pop <- gl.fst.pop(ak.pop.nm, nboots=100, percent=95, nclusters = 32)

write_rds(fst.pop,"fst.pop.rds")

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

ak.pop.nm.nosibs <- read_rds("ak.pop.nm.nosibs.rds")

fst.pop.nosibs <- gl.fst.pop(ak.pop.nm.nosibs, nboots=100, percent=95, nclusters = 32)

write_rds(fst.pop.nosibs,"fst.pop.nosibs.rds")

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken