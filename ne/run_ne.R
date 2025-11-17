args <- commandArgs(trailingOnly = TRUE)


n <- strtoi(args[1])
focal_pop <- args[2]
iteration <- strtoi(args[3])
neest_path <- getwd()

library(tidyverse)
library(dartR.base)
library(dartR.popgen)

start.time <- Sys.time()

sample_inds <- function(original,sample_size,popi){
  id_pops <- data.frame(id=original$ind.names,pop=original$pop)
  
  ids <- id_pops %>% 
    filter(pop==popi) %>% 
    slice_sample(n=sample_size) %>% 
    pull(id)
  
  gl.keep.ind(original,ids,mono.rm = TRUE,recalc = TRUE,verbose = FALSE)
}

ak.chr.mi.adj <- read_rds("ak.chr.mi.adj.rds")

si <- sample_inds(ak.chr.mi.adj,sample_size = n,focal_pop)

tf <- tempfile(fileext = ".txt")
fd <- dirname(tf)

si.ne <- dartR.popgen::gl.LDNe(si,outfile = basename(tf),outpath = fd, neest.path =  neest_path,critical = c(0, 0.05),plot.out = FALSE,pairing = "separate")

write_rds(si.ne,file = paste(focal_pop,"_",n,"_",iteration,".rds",collapse = "",sep = ""))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

