#' @author Lu Ao, Chongqing, China
#' @description
#' This is a script to convert the summary statistics data to the format that smr can use (i.e., mygwas.ma).
#' Part 1 is to pre-process the summary statistics data.
#' Part 2 is to generate the mygwas.ma file.
#' Part 3 is bash script to run smr.
#' Part 4 is to fdr adjust the smr results.

################# ################# Pre-process for ss data ################# ################# 
library(tidyverse);library(vroom)
gwas1_path <- "~/project/iridocyclitis/data/diabete/1/GCST90014023_buildhg19.tsv" # t1d1
gwas1 <- vroom(gwas1_path) %>% 
  select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF, Samplesize) %>% 
  set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF", "N") %>% 
  mutate(MAF = ifelse(EAF<0.5, EAF, 1-EAF)) %>% 
  filter(MAF > 0.01)

gwas2_path <- "~/project/iridocyclitis/data/iri_finn/ukb_finn/iri_fbmeta_hg19.tsv" # iri3
gwas2 <- vroom(gwas2_path) %>% 
  mutate(N_Cohort = case_when(is.na(FINNGEN_beta) & !is.na(UKBB_beta) ~ "UKB",
                              !is.na(FINNGEN_beta) & is.na(UKBB_beta) ~ "Finn",
                              !is.na(FINNGEN_beta) & !is.na(UKBB_beta) ~ "Finn+UKB",
                              TRUE ~ "Unknown"),
         N = case_when(N_Cohort == "UKB" ~ 2304+413959,
                       N_Cohort == "Finn" ~ 8016+390647,
                       N_Cohort == "Finn+UKB" ~ 2304+413959+8016+390647,
                       TRUE ~ NA),
         EAF = case_when(N_Cohort == "UKB" ~ UKBB_af_alt,
                         N_Cohort == "Finn" ~ FINNGEN_af_alt,
                         N_Cohort == "Finn+UKB" ~ (FINNGEN_af_alt*(8016+390647)+UKBB_af_alt*(2304+413959))/(8016+390647+2304+413959),
                         TRUE ~ NA)) %>% 
  dplyr::select(SNP, CHR, BP, ALT, REF, all_inv_var_meta_p, all_inv_var_meta_beta, all_inv_var_meta_sebeta, EAF, N) %>% 
  set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF", "N") %>%
  mutate(MAF = ifelse(EAF<0.5, EAF, 1-EAF)) %>% 
  filter(MAF > 0.01)

################################## organize files ################################## 
filter_min_p <- function(df) {
  df <- df %>% # mark who is duplicated
    mutate(dup = duplicated(SNP, fromLast = TRUE) | duplicated(SNP, fromLast = FALSE))
  df_dup <- df %>%  # filter duplicated SNP and find the ones with min p
    filter(dup) %>%
    group_by(SNP) %>%
    filter(P == min(P)) %>%
    dplyr::slice(1) %>% # keep the 1 when p is the same
    ungroup()
  df_non_dup <- df %>% # non-duplicated SNP
    filter(!dup)
  rbind(df_non_dup, df_dup) %>% # brilliant!
    select(-dup) # This way only check the duplicated ones
}

type1d.ma <- gwas1 %>% 
  dplyr::select(SNP, A1, A2, EAF, BETA, SE, P, N) %>% 
  drop_na(SNP) %>% 
  filter_min_p() %>% 
  set_names("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
vroom::vroom_write(type1d.ma, './data/smr/t1d1.ma')

iri.ma <- gwas2 %>% 
  dplyr::select(SNP, A1, A2, EAF, BETA, SE, P, N) %>% 
  drop_na(SNP) %>% 
  filter_min_p() %>% 
  set_names("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
vroom::vroom_write(iri.ma, './data/smr/iri3.ma')

################################## run smr ################################## 
#!/bin/zsh
#' # command line for smr
#' Initial param settings
beqtl=/Volumes/LaCie/pQTL/besd/besd-ukbppp/ukbppp.eur
gwas1=/Users/leoarrow/project/iridocyclitis/data/smr/t1d1.ma
gwas2=/Users/leoarrow/project/iridocyclitis/data/smr/iri3.ma
out1=/Users/leoarrow/project/iridocyclitis/output/smr/pqtl/ppp/smr_ppp@t1d1
out2=/Users/leoarrow/project/iridocyclitis/output/smr/pqtl/ppp/smr_ppp@iri3
#' type1dp
smr --bfile /Users/leoarrow/project/ref/1kg.v3/EUR \
--gwas-summary $gwas1 \
--beqtl-summary $beqtl \
--out $out1 \
--thread-num 16
#' uveitis
smr --bfile /Users/leoarrow/project/ref/1kg.v3/EUR \
--gwas-summary $gwas2 \
--beqtl-summary $beqtl \
--out $out2 \
--thread-num 16
# optional
--diff-freq-prop 0.35 \

#' eqtl
beqtl=/Users/leoarrow/project/ref/eqtl/eqtl_gen/cis-eQTL-SMR_20191212/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense
beqtl=/Users/leoarrow/project/ref/eqtl/yanglab/Cage/cage_eqtl_data/CAGE.sparse
beqtl=/Users/leoarrow/project/ref/eqtl/yanglab/Geuvadis/geuvadis_EUR_rsid/geuvadis_EUR_rsid
#' mqtl
beqtl=/Users/leoarrow/project/ref/mqtl/yanglab/US_mQTLS_SMR_format
beqtl=/Users/leoarrow/project/ref/caQTL/Bryois_caQTL_summary/bryois_NatCommun_2018_50kb_cQTLs
####################################### FDR correction ############################################
#' FDR and Bonferroni
#' @param smr_result_path smr result path #'/Users/leoarrow/project/iridocyclitis/output/smr_gtexV8wholeblood_type1d2.smr'
leo_smr_adjust <- function(smr_result_path){
  require(vroom)
  smr_result <- vroom::vroom(smr_result_path, show_col_types = F) %>% 
    mutate(Pass_HEIDI = ifelse(p_HEIDI>=0.05, "Pass", "Fail"),
           FDR = p.adjust(p_SMR, method = "BH"),
           Bonferroni = p.adjust(p_SMR, method = "bonferroni"))
  basename <- basename(smr_result_path) %>% gsub(".smr", ".fdr.smr", .); message(paste0(" - Basename set to : ", basename))
  pathname <- dirname(smr_result_path); message(paste0(" - Path set to : ", pathname))
  writePath <- paste0(pathname, "/", basename); message(paste0(" - Writing to : ", writePath))
  vroom::vroom_write(smr_result, writePath)
}
# 
# leo_smr_adjust("/Users/leoarrow/project/iridocyclitis/output/smr_gtexV8wholeblood_type1d.smr")
# leo_smr_adjust("/Users/leoarrow/project/iridocyclitis/output/smr_gtexV8wholeblood_iri_fkmeta.smr")

#' @param folder_path
#' @example leo_smr_adjust_loop("/Users/leoarrow/project/iridocyclitis/output/GTEx49_sqtl")
library(tidyverse)
library(vroom)
leo_smr_adjust_loop <- function(folder_path) {
  smr_files <- list.files(path = folder_path, pattern = "\\.smr$", full.names = TRUE) # Get list of .smr files in the folder
  for (smr_file in smr_files) { # Loop through each .smr file and apply leo_smr_adjust function
    message(paste0("*Processing: ", smr_file))
    leo_smr_adjust(smr_file)
  }
}
leo_smr_adjust_loop("/Users/leoarrow/project/iridocyclitis/output/smr/eqtl/GTEx49")
# leo_smr_adjust_loop("/home/luao/project/iridocyclitis/output/Hatton_mqtl/")

#' Important note: This part aims to extract sig res from `.fdr.smr` file; and merge them into one ################
#' 
#' @title `leo_smr_sig_extract`: Extracting sig res and send them into one dir (typically `fdr_sig`)
#' @param smr_result_path smr result path #'/Users/leoarrow/project/iridocyclitis/output/smr_gtexV8wholeblood_type1d2.smr'
#' note save with suffix `.sig.smr` in `fdr_sig` dir in the same dir that save `fdr.smr` file
#' 
leo_smr_sig_extract <- function(fdr_path, save_path = NULL){
  require(vroom)
  smr_result <- vroom::vroom(fdr_path,show_col_types = F) %>% 
    filter(Pass_HEIDI == "Pass", FDR < 0.05)
  basename <- basename(fdr_path) %>% gsub(".fdr.smr", ".sig.smr", .); message(paste0(" - Basename set to : ", basename))
  if (is.null(save_path)) {
    save_path <- dirname(fdr_path); message(paste0(" - Path set to : ", save_path))
  } else {
    save_path <- save_path; message(paste0(" - Path set to : ", save_path)) # highly recommended！Set it as "./fdr_sig"
  }
  writePath <- paste0(save_path, "/", basename); message(paste0(" - Writing to : ", writePath))
  vroom::vroom_write(smr_result, writePath)
}

# leo_smr_sig_extract(
#   fdr_path = "/Users/leoarrow/project/iridocyclitis/output/GTEx49_eqtl/Adipose_Subcutaneous_iri_meta.eqtl.fdr.smr",
#   save_path = "/Users/leoarrow/project/iridocyclitis/output/GTEx49_eqtl/fdr_sig/"
# )

#' Title `leo_smr_sig_loop` will preform `leo_smr_sig_extract` for the given dir
#'
#' @param folder_path 
#' @param save_path 
library(tidyverse)
library(vroom) # just in case
leo_smr_sig_loop <- function(folder_path, save_path = NULL) {
  smr_fdr_files <- list.files(path = folder_path, pattern = "\\.fdr.smr$", full.names = TRUE)
  if (!dir.exists(save_path)) {dir.create(save_path)}
  for (smr_fdr_file in smr_fdr_files) {
    message(paste0("*Processing: ", smr_fdr_file))
    leo_smr_sig_extract(smr_fdr_file, save_path = save_path)
  }
}

leo_smr_sig_loop(
  folder_path = "/Users/leoarrow/project/iridocyclitis/output/smr/eqtl/GTEx49",
  save_path = "/Users/leoarrow/project/iridocyclitis/output/smr/eqtl/GTEx49/fdr_sig"
)
# leo_smr_sig_loop(
#   folder_path = "/home/luao/project/iridocyclitis/output/Hatton_mqtl/",
#   save_path = "/home/luao/project/iridocyclitis/output/Hatton_mqtl/fdr_sig"
# )


# Merge and extract positive results
# Now refer to '2.smr_merge_sharedSigGene.R'

# Code below is legacy which we no longer recommend anymore.
# But it does also works.
# leo_smr_sig_merge <- function(folder_path) {
#   smr_sig_files <- list.files(path = folder_path, pattern = "\\.sig.smr$", full.names = TRUE) # Get list of .sig.smr files in the folder
#   smr_sig_files <- lapply(smr_sig_files, vroom::vroom) # Read all .sig.smr files
#   smr_sig_files <- bind_rows(smr_sig_files) # Merge all .sig.smr files
#   writePath <- paste0(folder_path, "/merged.sig.smr"); message(paste0(" - Writing to : ", writePath))
#   vroom::vroom_write(smr_sig_files, writePath)
# }
# leo_smr_sig_merge("/Users/leoarrow/project/iridocyclitis/output/eqtlgen/fdr_sig")