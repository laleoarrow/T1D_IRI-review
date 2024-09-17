#' This code contains script for GSMR analysis using GCTA analysis from Yanglab.
#' >>> The mr for each exposure contains 4 exposure and 3 outcome
#' >>> All data in hg38 were lifted hg19 before GSMR analysis.
pacman::p_load(vroom, data.table, tidyverse, dplyr, # for importing data
               TwoSampleMR, MendelianRandomization, ieugwasr, # for MR
               ggplot2, forestploter, # for gragh
               pbapply, parallel, tictoc # for parallel
)
#' >>> functions
# This one make sure no duplicated SNP
filter_min_p <- function(df) {
 df <- df %>% # mark who is duplicated
  mutate(dup = duplicated(SNP, fromLast = TRUE) | duplicated(SNP, fromLast = FALSE))
 df_dup <- df %>%  # filter duplicated SNP and find the ones with min p
  filter(dup) %>%
  group_by(SNP) %>%
  filter(p == min(p)) %>%
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
# make `gsmr_exposure.txt` and `gsmr_outcome.txt`
# e.g., (trait path); raw file in GCTA-COJO format (SNP A1 A2 freq b se p N)
# hdl hdl_test.raw
# bmi bmi_test.raw
out_path <- "./data/gcta/"
t1d_p1 <- "./data/diabete/1/GCST90014023_buildhg19.tsv"
t1d_p2 <- "./data/diabete/1_finn/summary_stats_finngen_R10_T1D_hg19.tsv"
t2d_p1 <- "./data/diabete/2/DIAMANTE-EUR.sumstat.txt"
t2d_p2 <- "./data/diabete/2_finn/finngen_R10_T2D_hg19.tsv"
iri_p1 <- "./data/iri_finn/r10/summary_stats_finngen_R10_H7_IRIDOCYCLITIS_hg19.tsv"
iri_p2 <- "~/project/iridocyclitis/data/iri_finn/ukb_finn/iri_fbmeta_hg19.tsv"
iri_p3 <- iri_p2
system(paste0("head -n 3 ", iri_p2))

# re-format to GCTA-COJO format
t1d1 <- vroom(t1d_p1) %>% 
  dplyr::select(SNP, A1, A2, EAF, BETA, SE, P, Samplesize) %>% 
  set_names("SNP", "A1", "A2", "freq", "b", "se", "p", "N") %>% 
  drop_na(SNP) %>% 
  filter_min_p() # check_duplicate_snps(t1d1)
vroom_write(t1d1, file.path(out_path, "T1D_Chiou.raw"), delim = " ")
t1d2 <- vroom(t1d_p2) %>% 
 dplyr::select(SNP, A1, A2, EAF, BETA, SE, P) %>% 
 set_names("SNP", "A1", "A2", "freq", "b", "se", "p") %>% 
 mutate(N = 4320+335112) %>% 
 drop_na(SNP) %>% 
 filter_min_p() # check_duplicate_snps(t1d2)
vroom_write(t1d2, file.path(out_path, "T1D_FinnGen.raw"), delim = " ")
t2d1 <- vroom(t2d_p1) %>% 
 dplyr::select(rsID, effect_allele, other_allele, effect_allele_frequency, `Fixed-effects_beta`, `Fixed-effects_SE`, `Fixed-effects_p-value`) %>% 
 set_names("SNP", "A1", "A2", "freq", "b", "se", "p") %>% 
 mutate(A1 = toupper(A1), A2 = toupper(A2),
        N = 80154+853816) %>% 
 drop_na(SNP) %>% 
 filter_min_p() # check_duplicate_snps(t2d1)
vroom_write(t2d1, file.path(out_path, "T2D_DIAMANTE.raw"), delim = " ")
t2d2 <- vroom(t2d_p2) %>% 
 dplyr::select(SNP, A1, A2, EAF, BETA, SE, P) %>% 
 set_names("SNP", "A1", "A2", "freq", "b", "se", "p") %>% 
 mutate(N = 65085+335112) %>% 
 drop_na(SNP) %>% 
 filter_min_p() # check_duplicate_snps(t2d2)
vroom_write(t2d2, file.path(out_path, "T2D_FinnGen.raw"), delim = " ")
# out
iri1 <- vroom(iri_p1) %>% 
 dplyr::select(SNP, A1, A2, EAF, BETA, SE, P) %>% 
 set_names("SNP", "A1", "A2", "freq", "b", "se", "p") %>% 
 mutate(N = 8016+390647) %>% 
 drop_na(SNP) %>% 
 filter_min_p() # check_duplicate_snps(iri1)
vroom_write(iri1, file.path(out_path, "Iridocyclitis_FinnGen.raw"), delim = " ")
iri2 <- vroom(iri_p2) %>% # "Iridocyclitis (UKB)"
  drop_na(UKBB_beta) %>% 
  dplyr::select(SNP, ALT, REF, UKBB_af_alt, UKBB_beta, UKBB_sebeta, UKBB_pval) %>% 
  set_names("SNP", "A1", "A2", "freq", "b", "se", "p") %>% 
  mutate(N = 2304+413959) %>% 
  filter_min_p()
vroom_write(iri2, file.path(out_path, "Iridocyclitis_UKB.raw"), delim = " ")
iri_3 <- vroom(iri_p3) %>% # `Iridocyclitis (FinnGen+UKB)`
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
  select(SNP, ALT, REF, EAF, all_inv_var_meta_beta, all_inv_var_meta_sebeta, all_inv_var_meta_p, N) %>% 
  set_names("SNP", "A1", "A2", "freq", "b", "se", "p", "N") %>% 
  filter_min_p()
vroom_write(iri_3, file.path(out_path, "Iridocyclitis_FinnGen_UKB.raw"), delim = " ")
system(paste0("head -n 3 ", file.path(out_path, "*.raw"))) # check

# make gsmr_exposure.txt and gsmr_outcome.txt
gcta_path <- "/Users/leoarrow/project/iridocyclitis/data/gcta"
t1d1_path <- file.path(gcta_path, "T1D_Chiou.raw")
t1d2_path <- file.path(gcta_path, "T1D_FinnGen.raw")
t2d1_path <- file.path(gcta_path, "T2D_DIAMANTE.raw")
t2d2_path <- file.path(gcta_path, "T2D_FinnGen.raw")
iri1_path <- file.path(gcta_path, "Iridocyclitis_FinnGen.raw")
iri2_path <- file.path(gcta_path, "Iridocyclitis_UKB.raw")
iri3_path <- file.path(gcta_path, "Iridocyclitis_FinnGen_UKB.raw")
system(paste0("head -n 3 ", t1d1_path))
gsmr_exposure <- c(paste0("T1D_(Chiou_et_al.) ", t1d1_path), # exposure
                   paste0("T1D_(FinnGen) ", t1d2_path),
                   paste0("T2D_(DIAMANTE) ", t2d1_path),
                   paste0("T2D_(FinnGen) ", t2d2_path))
gsmr_outcome <- c(paste0("Iridocyclitis_(FinnGen) ", iri1_path), # outcome
                  paste0("Iridocyclitis_(UKB) ", iri2_path),
                  paste0("Iridocyclitis_(FinnGen+UKB) ", iri3_path))
# write it
writeLines(gsmr_exposure, file.path(out_path, "gsmr_exposure.txt"))
writeLines(gsmr_outcome, file.path(out_path, "gsmr_outcome.txt"))

# others
gmb_path <- "/Users/leoarrow/project/iridocyclitis/data/gmb/GCST90011362_liftover.h19.h.tsv.gz"
gmb1 <- vroom(gmb_path) %>% 
  dplyr::select(SNP, A1, A2, EAF, BETA, SE, P) %>% 
  set_names("SNP", "A1", "A2", "freq", "b", "se", "p") %>% 
  drop_na(SNP) %>% 
  mutate(N = 8956)
  filter_min_p() # check_duplicate_snps(t1d1)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Formal GSMR   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# https://yanglab.westlake.edu.cn/software/gcta/#MendelianRandomisation

# command line for gsmr analysis using gcta
gsmr_exposure=/Users/leoarrow/project/iridocyclitis/data/gcta/gsmr_exposure.txt
gsmr_outcome=/Users/leoarrow/project/iridocyclitis/data/gcta/gsmr_outcome.txt
gcta \
--gsmr2-beta --heidi-thresh 0.01 \
--bfile /Users/leoarrow/project/ref/1kg.v3/EUR \
--gsmr-file $gsmr_exposure $gsmr_outcome \
--gsmr-direction 0 \
--effect-plot \
--thread-num 10 \
--out /Users/leoarrow/project/iridocyclitis/output/gsmr/gsmr2_result

# optional: a single pair gsmr analysis and other command options
gsmr_exposure=/Users/leoarrow/project/iridocyclitis/data/gcta/ukb_demo/gsmr_exposure.txt
gsmr_outcome=/Users/leoarrow/project/iridocyclitis/data/gcta/ukb_demo/gsmr_outcome.txt
gcta \
--bfile /Users/leoarrow/project/ref/1kg.v3/EUR \
--gsmr-file $gsmr_exposure $gsmr_outcome \
--gsmr-direction 2 \
--clump-r2 0.001 \
--gsmr2-beta --heidi-thresh 0.01 \
--threads 10 \
--thread-num 10 \
--out /Users/leoarrow/project/iridocyclitis/data/gcta/ukb_demo/demo_result/ukbdemo_result

# gsmr plot ---------------------------------------------------------------
source("/Users/leoarrow/project/software/gcta-1.94.2-MacOS-ARM-x86_64/gsmr_plot.r")
gsmr_data = read_gsmr_data("/Users/leoarrow/project/iridocyclitis/output/gsmr/gsmr2_result.eff_plot.gz")
gsmr_summary(gsmr_data)      # show a summary of the data
results = gsmr_data[["bxy_result"]] %>% as.data.frame() %>% filter(!(grepl("Finn", Exposure) & grepl("Finn", Outcome)))
results

plot_smr <- function(i, save = T){ # save = T: save each plot locally; F return the plot for patch
  # get param
  exp_str <- results$Exposure[i]
  out_str <- results$Outcome[i]
  p_value <- signif(as.numeric(results$p[i]), 3)
  output_name <- paste0(exp_str, "_VS_", out_str)
  # plot
  if (save) {pdf(paste0("./figure/gsmr/", output_name, ".pdf"), width = 6, height = 6)}
  plot_gsmr_effect(gsmr_data, exp_str, out_str, "#00A087FF")
  # annotate P
  x_pos <- grconvertX(0.2, "npc", "user")
  y_pos <- grconvertY(0.8, "npc", "user")
  text(x = x_pos, y = y_pos, labels = bquote(italic(P) == .(p_value)), cex = 1.2)
  # annotate ABCD..
  if (!save) {
    label_x <- grconvertX(0.05, "npc", "user")
    label_y <- grconvertY(0.95, "npc", "user")
    text(x = label_x, y = label_y, labels = LETTERS[i], cex = 1.5, font = 2)
  }
  if (save) {dev.off()} else {return(recordPlot())}
}
# single
pbapply::pblapply(1:dim(results), function(i) plot_smr(i, save = F))
# all
opar <- par(no.readonly = TRUE)
pdf("./figure/gsmr/all_plots_T1D.pdf", width = 8, height = 8)
par(mfrow = c(2,2))
pbapply::pblapply(1:4, function(i) plot_smr(i, save = FALSE))
dev.off()
