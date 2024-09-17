#' @Author Lu Ao
#' @description
#' This is a pipline to load your gwas summary statistics and convert its BP using MungeSumstats packages
#'  - `MungeSumstats::lifeover` were used here
#'  - This script is designed for the conversion for any gwas ss data
#'  - Need to be modified for different gwas data
#' @note
#' Tips: use tidyverse to clean the ss data first. This would avoid you with many unexpected problems.
#' for de-bugging: rm(list=ls());gc()
#' detach("package::tidyverse",unload=T)

library(MungeSumstats) # https://neurogenomics.github.io/MungeSumstats/reference/liftover.html
library(tidyverse)
library(liftOver)
library(vroom)

gwas1_path <- "~/project/iridocyclitis/data/iri_finn/ukb_finn/iri_fbmeta_hg38.tsv"
gwas1 <- vroom(gwas1_path)
colnames(gwas1)

devtools::install_github("laleoarrow/MungeSumstats_ModifedLocalChain") # this is a modified version of MungeSumstats with local chain function provided
# https://github.com/neurogenomics/MungeSumstats/blob/master/R/liftover.R
gwas1_38 <- gwas1 %>% dplyr::rename(SNP = `rsid`, CHR = `#CHR`, BP = `POS`) # gwas1 <- gwas1 %>% dplyr::select(-SNP)
gwas1_38 <- gwas1_38 %>% dplyr::select(SNP, CHR, BP, everything()) # smart to re-arrange columns
gwas1_37 <- MungeSumstatsModified::liftover(sumstats_dt = gwas1_38,
                                            local_chain = "/Users/leoarrow/project/ref/liftover/hg38ToHg19.over.chain",
                                            ref_genome = "hg38",
                                            convert_ref_genome="hg19",
                                            verbose = T)
out_path <- gsub("hg38", "hg19_no_maf_filtering", gwas1_path); message(paste0("Save to: ", out_path))
MungeSumstats::write_sumstats(sumstats_dt = gwas1_37, save_path = out_path, nThread = 10)

### ldsc format ###
# MungeSumstats::write_sumstats(sumstats_dt = gwas1_37, save_path = "./data/finn/iri_fbmeta_hg19_ldsc", save_format = "LDSC", nThread = 10)


# The following is for lift type2d data from 37 to 38--------------------------------
t2d_p <- "./data/diabete/2/DIAMANTE-EUR.sumstat.txt" # it is build in 37 (i.e., hg19)
t2d <- vroom(t2d_p) %>% 
  mutate(MAF = ifelse(effect_allele_frequency<0.5, effect_allele_frequency, 1-effect_allele_frequency)) %>% 
  filter(MAF > 0.01)
# Reordering so first three column headers are SNP, CHR and BP in this order.
t2d <- t2d %>% dplyr::rename(SNP = `rsID`, CHR = `chromosome(b37)`, BP = `position(b37)`)
t2d <- t2d %>% dplyr::select(SNP, CHR, BP, everything()) # smart to re-arrange columns
t2d_38 <- MungeSumstatsModified::liftover(sumstats_dt = t2d,
                                            local_chain = "/Users/leoarrow/project/ref/liftover/hg19ToHg38.over.chain",
                                            ref_genome = "hg19",
                                            convert_ref_genome="hg38",
                                            verbose = T)
save_path <- "./data/diabete/2/DIAMANTE-EUR_hg38.sumstat.txt"
MungeSumstats::write_sumstats(sumstats_dt = t2d_38,save_path = save_path, nThread = 10)

# For Iridocyclitis (Gelfman et al.)------------------------------------------------
iri_p <- "./data/iri_nc_large_meta/GCST90295958_withN_rsid.tsv" # it is build in 38
iri <- vroom(iri_p) %>% 
 mutate(MAF = ifelse(effect_allele_frequency<0.5, effect_allele_frequency, 1-effect_allele_frequency)) %>% 
 filter(MAF > 0.01)
# Reordering so first three column headers are SNP, CHR and BP in this order.
iri <- iri %>% dplyr::rename(SNP = `RefSNP_id`, CHR = `chromosome`, BP = `base_pair_location`)
iri <- iri %>% dplyr::select(SNP, CHR, BP, everything()) # smart to re-arrange columns
iri_37 <- MungeSumstatsModified::liftover(sumstats_dt = iri,
                                          local_chain = "/Users/leoarrow/project/ref/liftover/hg38ToHg19.over.chain",
                                          ref_genome = "hg38",
                                          convert_ref_genome="hg19",
                                          verbose = T)
save_path <- "./data/iri_nc_large_meta/GCST90295958_withN_rsid_hg19.tsv"
MungeSumstats::write_sumstats(sumstats_dt = iri_37, save_path = save_path, nThread = 15)
