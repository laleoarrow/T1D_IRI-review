#' Adjust summary data into `fuma` required fromat
library(data.table); library(tidyverse)
res <- fread("/Users/leoarrow/project/iridocyclitis/output/placo/t1d1_iri3_placo.tsv.gz")
res <- res %>% select("SNP","CHR","POS","A1","A2","p.placo") %>% 
  mutate(POS = as.integer(POS))
fwrite(res, "/Users/leoarrow/project/iridocyclitis/output/placo/t1d1_iri3_placo_fuma.tsv.gz", nThread = 16)