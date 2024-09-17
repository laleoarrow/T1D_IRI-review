#' This file contains both R and Bash code used in translating and formating the ldsc required files.
#' The ldsc is now implemented in the server now.
# Ao Lu, Chongqing, China

# LDSC ------------------------------------------------------------
# Required columns: "SNP", "CHR", "BP", "A1", "A2", "BETA", "P"
# bin/R!
library(vroom); library(data.table); library(tidyverse)
t1d_p1 <- "~/project/iridocyclitis/data/diabete/1/GCST90014023_buildhg19.tsv"
t1d_p2 <- "~/project/iridocyclitis/data/diabete/1_finn/summary_stats_finngen_R10_T1D_hg19.tsv"
t2d_p1 <- "~/project/iridocyclitis/data/diabete/2/DIAMANTE-EUR.sumstat.txt"
t2d_p2 <- "~/project/iridocyclitis/data/diabete/2_finn/finngen_R10_T2D_hg19.tsv"
iri_p1 <- "~/project/iridocyclitis/data/iri_finn/r10/summary_stats_finngen_R10_H7_IRIDOCYCLITIS_hg19.tsv"
iri_p2 <- "~/project/iridocyclitis/data/iri_finn/ukb_finn/iri_fbmeta_hg19.tsv"
iri_p3 <- iri_p2
system(paste0("head -n 3 ", iri_p3))
hm3_path <- "/Users/leoarrow/project/software/ldsc-b站/w_hm3.snplist"
snp.list <- fread(hm3_path) %>% select(SNP)
#' @note This function is used to filter the duplicated SNP and keep the one with min p value.
#' @param df datafame with `SNP` and `P` column.
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
check_duplicate_snps <- function(df) {
  dup_snps <- df[duplicated(df$SNP), "SNP"]
  if (nrow(dup_snps) > 0) {
    message("Found duplicated SNP IDs: ", paste(dup_snps$SNP, collapse = ", "))
  } else {
    message("No duplicated SNP IDs found.")
  }
}

# t1d
t1d_1 <-  vroom(t1d_p1, show_col_types = FALSE) %>% # `T1D (Chiou et al.)`
  select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF, Samplesize) %>% 
  set_names("SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "SE", "EAF", "N") %>% # N for observed sample size
  filter(SNP %in% snp.list$SNP) %>%
  filter_min_p()
t1d_2 <-  vroom(t1d_p2, show_col_types = FALSE) %>% # `T1D (FinnGen)`
  select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF) %>% 
  set_names("SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "SE", "EAF") %>%
  mutate(N = 4320+335112) %>% 
  filter(SNP %in% snp.list$SNP) %>%
  filter_min_p()

# t2d
t2d_1 <- vroom(t2d_p1, show_col_types = FALSE) %>% # `T2D (DIAMANTE)`
  select(rsID, `chromosome(b37)`, `position(b37)`, effect_allele, other_allele, `Fixed-effects_p-value`, `Fixed-effects_beta`, `Fixed-effects_SE`, effect_allele_frequency ) %>% 
  set_names("SNP", "CHR", "BP", "A1", "A2","P", "BETA", "SE", "EAF") %>% 
  mutate(A1 = toupper(A1), A2 = toupper(A2),
         N = 80154+853816) %>% 
  filter(SNP %in% snp.list$SNP) %>%
  filter_min_p()
t2d_2 <- vroom(t2d_p2, show_col_types = FALSE) %>% # `T2D (FinnGen)`
  select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF) %>% 
  set_names("SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "SE", "EAF") %>%
  mutate(N = 65085+335112) %>% 
  filter(SNP %in% snp.list$SNP) %>%
  filter_min_p()

# iri
iri_1 <- vroom(iri_p1, show_col_types = FALSE) %>% # `Iridocyclitis (FinnGen)`
  select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF) %>% 
  set_names("SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "SE", "EAF") %>%
  mutate(N = 8016+390647) %>% 
  filter(SNP %in% snp.list$SNP) %>%
  filter_min_p()
iri_2 <- vroom(iri_p2, show_col_types = FALSE) %>% # `Iridocyclitis (UKB)`
  drop_na(UKBB_beta) %>% 
  select(SNP, CHR, BP, ALT, REF, UKBB_pval, UKBB_beta, UKBB_sebeta, UKBB_af_alt) %>% 
  set_names("SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "SE", "EAF") %>%
  mutate(N = 2304+413959) %>% 
  filter(SNP %in% snp.list$SNP) %>%
  filter_min_p()
iri_3 <- vroom(iri_p3, show_col_types = FALSE) %>% # `Iridocyclitis (FinnGen+UKB)`
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
  select(SNP, CHR, BP, ALT, REF, all_inv_var_meta_p, all_inv_var_meta_beta, all_inv_var_meta_sebeta, EAF, N) %>% 
  set_names("SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "SE", "EAF", "N") %>%
  filter(SNP %in% snp.list$SNP) %>%
  filter_min_p()

# write file in `summary.txt` format
write_ldsc.txt <- function(df, path) {
  df %>% vroom_write(., path, delim = "\t", col_names = T)
}
write_ldsc.txt(t1d_1, "./data/ldsc/light/t1d_1.txt")
write_ldsc.txt(t1d_2, "./data/ldsc/light/t1d_2.txt")
write_ldsc.txt(t2d_1, "./data/ldsc/light/t2d_1.txt")
write_ldsc.txt(t2d_2, "./data/ldsc/light/t2d_2.txt")
write_ldsc.txt(iri_1, "./data/ldsc/light/iri_1.txt")
write_ldsc.txt(iri_2, "./data/ldsc/light/iri_2.txt")
write_ldsc.txt(iri_3, "./data/ldsc/light/iri_3.txt")
gc()

#!/bin/bash------------------------------------------------------------
################### server version ###################
#' @note Munge for a single file >>>>>>
conda activate ldsc # this is important for ldsc
ldsc=/home/luao/software/ldsc/ldsc-master

txt=t2d_1
S1=/home/luao/project/iridocyclitis/data/ldsc/light/${txt}.txt
out1=/home/luao/project/iridocyclitis/data/ldsc/light/${txt}.munged
$ldsc/munge_sumstats.py --sumstats $S1 \
--chunksize 500000 \
--out $out1 \
--merge-alleles /home/luao/software/ldsc/ldsc-master/w_hm3.snplist

#' @note Munge for a whole directory  >>>>>>
conda activate ldsc # this is important for ldsc
ldsc=/home/luao/software/ldsc/ldsc-master
input_dir="/home/luao/project/iridocyclitis/data/ldsc/light"
out="/home/luao/project/iridocyclitis/data/ldsc/light"
for file in $input_dir/*.txt; do
  filename=$(basename $file .txt)
  output_file=${out}/${filename}.munged
  $ldsc/munge_sumstats.py --sumstats "$file" \
    --chunksize 500000 \
    --out "$output_file" \
    --merge-alleles /home/luao/software/ldsc/ldsc-master/w_hm3.snplist
done

#' @note calculate LDSC for two traits >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
conda activate ldsc # this is important for ldsc
ldsc=/home/luao/software/ldsc/ldsc-master
element_exp=t1d_1
element_out=iri_1 iri_2 iri_3
trait1ldsc=/home/luao/project/iridocyclitis/data/ldsc/light/${element_exp}.munged.sumstats.gz
trait2ldsc=/home/luao/project/iridocyclitis/data/ldsc/light/${element_out}.munged.sumstats.gz
out=/home/luao/project/iridocyclitis/output/ldsc/light/${element_exp}@${element_out}
$ldsc/ldsc.py --rg $trait1ldsc,$trait2ldsc \
--ref-ld-chr /home/luao/software/ldsc/eur_w_ld_chr_nomhc/ \
--w-ld-chr /home/luao/software/ldsc/eur_w_ld_chr_nomhc/ \
--out ${out}_nomhc

# loop version
element_exp=("t1d_2" "t2d_2")
element_out=("iri_2")

for exp_element in "${element_exp[@]}"
do
trait1ldsc="/home/luao/project/iridocyclitis/data/ldsc/light/${exp_element}.munged.sumstats.gz"
  for out_element in "${element_out[@]}"
  do
  trait2ldsc="/home/luao/project/iridocyclitis/data/ldsc/light/${out_element}.munged.sumstats.gz"
  out="/home/luao/project/iridocyclitis/output/ldsc/light/${exp_element}@${out_element}"
  echo " >>>>>>>>>>>> Processing: $trait1ldsc and $trait2ldsc >>>>>>>>>>>> "
  
  $ldsc/ldsc.py --rg $trait1ldsc,$trait2ldsc \
  --ref-ld-chr /home/luao/software/ldsc/eur_w_ld_chr_nomhc/ \
  --w-ld-chr /home/luao/software/ldsc/eur_w_ld_chr_nomhc/ \
  --out ${out}_nomhc
  done
done

# extract LDSC result from the log file ------------------------------------------------

#' @title get log file results
extract_results <- function(file_path) {
  lines <- readLines(file_path)
  summary_index <- grep("Summary of Genetic Correlation Results", lines)
  if (length(summary_index) == 0) {return(NULL)}
  col_names <- str_split(lines[summary_index + 1], "\\s+", simplify = TRUE)
  results <- str_split(lines[summary_index + 2], "\\s+", simplify = TRUE)
  result_df <- data.frame(matrix(ncol = length(col_names), nrow = 1))
  colnames(result_df) <- col_names
  for (i in 1:length(col_names)) {
    result_df[1, i] <- results[i]
  }
  return(result_df)
}

# start by put all *.log in one dir
library(dplyr); library(readr)
log_dir <- "/Users/leoarrow/project/iridocyclitis/output/ldsc/"
output_dir <- "/Users/leoarrow/project/iridocyclitis/output/ldsc/"
log_files <- list.files(log_dir, pattern = "\\.log$", full.names = TRUE)
results_list <- lapply(log_files, extract_results)
final_results <- do.call(rbind, results_list)
final_results <- final_results %>% 
  mutate(p1 = gsub(".munged.sumstats.gz", "", basename(p1)),
         p2 = gsub(".munged.sumstats.gz", "", basename(p2)))
final_results

output_dir <- "/Users/leoarrow/project/iridocyclitis/output/ldsc/"
output_file <- paste0(output_dir, "ldsc_results.csv")
write_csv(final_results, output_file)
cat("Summary table created at", output_file, "\n")

# optional: add the Meta data information from Munge step
library(stringr)
library(dplyr)

#' @title get munge log Metadata
extract_metadata <- function(file_path) {
  lines <- readLines(file_path) # file_path = log_files[1]
  metadata_index <- grep("Metadata:", lines)
  
  if (length(metadata_index) == 0) return(NULL)
  
  metadata_lines <- lines[(metadata_index + 1):(metadata_index + 4)]
  valid_lines <- metadata_lines[grepl(" = ", metadata_lines)]
  
  metadata_values <- sapply(valid_lines, function(line) {
    split_line <- strsplit(line, " = ")[[1]]
    as.numeric(trimws(split_line[2]))
  })
  
  trait <- basename(file_path) %>% gsub(".munged.log", "", .)
  
  data.frame(
    trait = trait,
    `Mean chi^2` = metadata_values[1],
    `Lambda GC` = metadata_values[2],
    `Max chi^2` = metadata_values[3]
  )
}

# main
log_files <- list.files(path = "/Users/leoarrow/project/iridocyclitis/output/ldsc/", pattern = "\\.munged.log$", full.names = TRUE)
metadata_results_list <- lapply(log_files, extract_metadata)
final_metadata_results <- do.call(rbind, metadata_results_list)
rownames(final_metadata_results) <- NULL
print(final_metadata_results)

# optional: combine the metadata with the ldsc original results
final_metadata_results_p1 <- final_metadata_results %>%
  rename(
    p1 = trait,
    p1_meta_Mean_chi2 = Mean.chi.2,
    p1_meta_Lambda_GC = Lambda.GC,
    p1_meta_Max_chi2 = Max.chi.2
  )

final_metadata_results_p2 <- final_metadata_results %>%
  rename(
    p2 = trait,
    p2_meta_Mean_chi2 = Mean.chi.2,
    p2_meta_Lambda_GC = Lambda.GC,
    p2_meta_Max_chi2 = Max.chi.2
  )

final_results <- final_results %>%
  left_join(final_metadata_results_p1, by = "p1") %>%
  left_join(final_metadata_results_p2, by = "p2"); final_result
output_file <- gsub(".results", ".with.metadata.results", output_file)
write_csv(final_results, output_file)
cat("Summary table created at", output_file, "\n")

