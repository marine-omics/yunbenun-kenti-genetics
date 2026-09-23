args <- commandArgs(trailingOnly = TRUE)


remove_id <- args[1]
focal_pop <- args[2]
iteration <- strtoi(args[3])
neest_path <- getwd()

library(tidyverse)
library(dartR.base)
library(dartR.popgen)

start.time <- Sys.time()

subset_pop <- function(original,remove_id,popi){
  id_pops <- data.frame(id=original$ind.names,pop=original$pop)

  ids <- id_pops %>% 
    filter(pop==popi) %>% 
    filter(id!=remove_id) %>% 
    pull(id)
  
  gl.keep.ind(original,ids,mono.rm = TRUE,recalc = TRUE,verbose = FALSE)
}

ak.chr <- read_rds("ak.chr.rds")

si <- subset_pop(ak.chr,remove_id,focal_pop)

tf <- tempfile(fileext = ".txt")
fd <- dirname(tf)

si.ne <- dartR.popgen::gl.LDNe(si,outfile = basename(tf),outpath = fd, neest.path =  neest_path,critical = c(0, 0.01,0.02,0.05,0.1),plot.out = FALSE,pairing = "separate")

write_rds(si.ne,file = paste("jackknife_sep/",focal_pop,"_",iteration,".rds",collapse = "",sep = ""))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

