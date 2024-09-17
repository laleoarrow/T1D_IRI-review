#! /usr/local/bin/Rscript
# This Rscript is used for FinnGen data to be liftovered from hg38 to hg19

library(MungeSumstats) 
library(tidyverse)
library(liftOver)
library(vroom)


args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {stop("Provide path here")}

finn_path <- args[1]

finn_gwas <- vroom(finn_path) %>%
 dplyr::select(rsids, `#chrom`, pos, alt, ref, af_alt, beta, sebeta, pval, nearest_genes) %>%
 set_names(c("SNP", "CHR", "POS", "A1", "A2", "EAF", "BETA", "SE", "P", "GENE")) %>%
 drop_na(CHR, POS, BETA) %>% # just in case
 mutate(MAF = ifelse(EAF < 0.5, EAF, 1 - EAF))
finn_gwas <- finn_gwas %>% filter(MAF > 0.01)

# Reordering so first three column headers are SNP, CHR and BP in this order.
finn_gwas <- finn_gwas %>% dplyr::rename(BP = POS)
finn_gwas_37 <- MungeSumstats::liftover(
 sumstats_dt = finn_gwas,
 local_chain = "/Users/leoarrow/project/ref/liftover/hg38ToHg19.over.chain",
 ref_genome = "hg38",
 convert_ref_genome = "hg19",
 verbose = TRUE
)


output_path <- sub(".gz", "_hg19.tsv", finn_path)
message(paste0("Writing to >>> ", getwd(), output_path))
MungeSumstats::write_sumstats(
 sumstats_dt = finn_gwas_37, 
 save_path = output_path,
 nThread = 10
)