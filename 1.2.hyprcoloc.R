# This file contains both R and Bash code used for coloc analysis.
# Ao Lu, 2024-03-17, Chongqing, China
# rm(list = ls());gc() # setwd("../")
# a=ls();rm(list = setdiff(a, c("T1D", "Iridocyclitis")));gc()

# bin/R!------------------------------------------------------------
source("./code/1.2.leo_coloc_packages.R")
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(vroom))
suppressMessages(library(hyprcoloc))


################## ` Note before actually start #############################
# Please do include all SNPs you have in the region
# Do not trim by significance or MAF or other means 
# (https://chr1swallace.github.io/coloc/articles/a02_data.html)

# * Read my gwas
T1D_path <- "~/project/iridocyclitis/data/diabete/1/GCST90014023_buildhg19.tsv" # t1d1
T1D <- vroom(T1D_path, show_col_types = FALSE) %>% 
  dplyr::select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF, Samplesize) %>% 
  set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF", "N") %>%
  mutate_at("POS", as.integer)

IRI_path <- "~/project/iridocyclitis/data/iri_finn/ukb_finn/iri_fbmeta_hg19.tsv" # iri3
Iridocyclitis <- vroom(IRI_path, show_col_types = FALSE) %>% 
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
  mutate_at("POS", as.integer)

# exclue MHC
T1D <- T1D %>% filter(!(CHR == 6 & POS >= 25000000 & POS <= 34000000))
Iridocyclitis <- Iridocyclitis %>% filter(!(CHR == 6 & POS >= 25000000 & POS <= 34000000))

# placo lead snp
snp_annotated <- data.table::fread("/Users/leoarrow/project/iridocyclitis/output/fuma/placo/FUMA_job493590/[snps].txt")
lead <- data.table::fread("/Users/leoarrow/project/iridocyclitis/output/fuma/placo/FUMA_job493590/【leadSNPs】.txt") %>% 
  dplyr::left_join(snp_annotated %>% dplyr::select(rsID, nearestGene), by = c("rsID" = "rsID")) %>% 
  dplyr::select(rsID, chr, pos, nearestGene) %>%
  set_names("SNP", "CHR", "POS", "nearestGene")
# loci for hyprcoloc
coloc_loci <- fread("/Users/leoarrow/project/iridocyclitis/output/coloc/coloc_results.tsv") %>% 
  mutate(PP34=PP3+PP4) %>% dplyr::filter(PP34>=0.5) %>% arrange(desc(PP34)) %>% dplyr::select(locus, PP3, PP4, gene) %>% 
  mutate(evidence = "PLACO_COLOC",
         chr = sub(":.*", "", locus),
         bp_l = as.numeric(sub(".*:(\\d+)-.*", "\\1", locus)),
         bp_u = as.numeric(sub(".*-(\\d+)", "\\1", locus)),
         probBP = (bp_l+bp_u)/2) %>% 
  distinct(locus, .keep_all = TRUE) %>% 
  dplyr::select(-bp_l, -bp_u) %>% 
  mutate(chr = as.integer(chr),
         Note = sprintf("%.2f/%.2f", PP3, PP4)) %>% 
  dplyr::select(evidence, chr, probBP, gene, Note) %>% 
  set_names("Evidence", "CHR", "BP", "GENE", "Note(PP3/4;qtl_type)")
smr_loci <- fread("/Users/leoarrow/project/iridocyclitis/output/smr/smr_all.figure_nonmhc.tsv") %>% 
  dplyr::select(type, gene,ProbeChr,Probe_bp) %>% 
  mutate(evidence = "SMR") %>% 
  dplyr::select(evidence, ProbeChr, Probe_bp, gene, everything()) %>% 
  distinct(Probe_bp, .keep_all = TRUE) %>% 
  set_names("Evidence", "CHR", "BP", "GENE", "Note(PP3/4;qtl_type)") %>% 
  dplyr::filter(!(CHR = 6 & BP > 25000000 & BP < 34000000))
loci <- bind_rows(coloc_loci, smr_loci)
rm(list = c("snp_annotated", "lead", "coloc_loci", "smr_loci")); gc()

# Reduce the raw data size using loci information
reduce_data <- function(gwas_data, loci, win=500) {
  reduced_data <- do.call(rbind, lapply(1:nrow(loci), function(i) {
    message(paste0("Selecting loci #", i))
    chr <- loci$CHR[i]; bp <- loci$BP[i]
    bp_lower <- as.integer(bp-(win*1000)/2); bp_upper <- as.integer(bp+(win*1000)/2)
    gwas_data %>% dplyr::filter(CHR == chr & POS >= bp_lower & POS <= bp_upper)
  })) %>% distinct() %>% dplyr::filter(!is.na(SNP))
  return(reduced_data)
}
T1D <- reduce_data(T1D, loci, win=500)
Iridocyclitis <- reduce_data(Iridocyclitis, loci, win=500)
gc()

# hyprcoloc 
# 473 ----
#'@note
#' * T1D Iridocyclitis loci already loaded
#' * dir for gut: /Users/leoarrow/project/ref/mbqtl/473/473hg19/*.tsv.gz # gmb for short
pacman::p_load(pbapply, data.table, tidyverse, ggplot2, parallel, future, future.apply)
source("./code/1.2.leo_coloc_packages.R")
gmb_files <- list.files(path = "/Users/leoarrow/project/ref/mbqtl/473/473hg19", pattern = "\\.tsv.gz$", full.names = T)
gwas_list <- list(T1D = T1D, Iridocyclitis = Iridocyclitis); rm(list = c("T1D", "Iridocyclitis")); gc()
# x <- perform_hypr(gmb_files, gwas_list, prefix = "473_", # let 473 be
#              save.plot.dir = "/Users/leoarrow/project/iridocyclitis/figure/hyprcoloc",
#              save.output.dir = "/Users/leoarrow/project/iridocyclitis/output/hyprcoloc",
#              select = c("rsID", "chr", "bp", "beta", "SE", "P.weightedSumZ"))
#----------------------------------main----------------------------------------#
res_stage1s <- list(); save.dir <- "/Users/leoarrow/project/iridocyclitis/figure/hyprcoloc";time1 <- Sys.time()
for (j in length(gmb_files):1) {
  loop_start_time <- Sys.time()
  time_spent <- loop_start_time - time1
  leo_message(paste0("\n# Loop <", j, "> Start Time: ", format(loop_start_time, "%Y-%m-%d %H:%M:%S"), " #"))
  #------------------outer loop-------------------#
  gmb_path <- gmb_files[j]
  gmb_id <- basename(gmb_path) %>% gsub("hg19.tsv.gz","",.)
  gmb_file <- fread(gmb_path) %>%  mutate_at("CHR", as.integer) # %>% reduce_data(loci, win=500)
  leo_message(paste("\n# Processing #", j, ": ", gmb_id))
  #------------------inner loop-------------------#
  res_one_gmbs <- list()
  for (i in 1:nrow(loci)) {
    # initial parameters
    # gwas_list <- list(T1D = T1D, Iridocyclitis = Iridocyclitis); length(gwas_list)
    chr <- loci$CHR[i]
    bp <- loci$BP[i]
    win <- 500 # kb
    gene <- loci$GENE[i]
    Evidence <- ifelse(loci$Evidence[i]=="PLACO_COLOC","PC",loci$Evidence[i])
    reminder <- paste0(Evidence,"_", gmb_id, "@", loci$CHR[i],":",loci$BP[i]-win*1000/2,"-",loci$BP[i]+win*1000/2)
    # format hyprcoloc input
    hypr_loci_dat <- format_hyprc_dat(gwas_list, chr, bp, win=win,
                                      gmb=T, gmb_file=gmb_file, gmb_id=gmb_id)
    beta_matrix <- hypr_loci_dat$beta_matrix
    se_matrix <- hypr_loci_dat$se_matrix
    snp_id = rownames(beta_matrix); leo_message(paste("\nSNP length:", length(snp_id)))
    if (length(snp_id) == 0) {next}
    # hyprcoloc
    binary_traits = c(1,1,0)
    hypr_result <- hyprcoloc::hyprcoloc(beta_matrix, se_matrix,
                                        trait.names = colnames(beta_matrix), snp.id = rownames(beta_matrix),
                                        binary.outcomes = binary_traits)
    reg.res <- hypr_result$results %>%
      mutate(index = reminder, chr = chr, bp = bp, gmb_id = gmb_id, gene = gene) %>%
      select(index, chr, bp, gmb_id, gene, everything())
    # sensitivity
    .check_hyprcoloc_postive <- function(reg.res) {
      coloced_traits <- strsplit(reg.res$traits, ",")[[1]] %>% length()
      logi1 <- nrow(reg.res) == 1
      # logi2 <- reg.res$traits != "None"
      # logi3 <- !is.na(reg.res$posterior_prob) && reg.res$posterior_prob > 0.25
      logi4 <- coloced_traits == 3
      logi = logi1 && logi4
      if (logi) {return(TRUE)} else {return(FALSE)}
    }
    if (.check_hyprcoloc_postive(reg.res)) {
      leo_message(paste0("# Succeed to hyprcoloc for loci: ", reminder, "!!!!!!!!!"), "33")
      sen <- hyprcoloc::sensitivity.plot(beta_matrix, se_matrix,
                                         trait.names = colnames(beta_matrix), snp.id = rownames(beta_matrix),
                                         reg.thresh = c(0.5, 0.6, 0.7), align.thresh = c(0.5, 0.6, 0.7), prior.c = c(0.05, 0.02, 0.01, 0.005),
                                         equal.thresholds = FALSE, similarity.matrix = TRUE)
      sim.mat <- sen[[2]] # for self-defined pheatmap
      save.dir <- save.dir
      title <- reminder %>% gsub(paste0("_",gmb_id),"",.)
      save.path <- paste0(save.dir, "/", reminder, ".pdf")
      p <- self_draw_hypr(sim.mat, title, save = T, save.path = save.path) # grid::grid.draw(p$gtable)
    } else {message(paste0("# Failed to hyprcoloc for loci: ", reminder))}
    # clean
    res_one_gmbs[[i]] <- reg.res
    leo_message("\n# clean cache #");cat(gc())
  }
  #------------------inner loop-------------------#
  res_one_gmb <- do.call("rbind", res_one_gmbs)
  res_stage1s[[j]] <- res_one_gmb
  rm(gmb_file)
  leo_message("\n# clean cache #");cat(gc())
  #------------------outer loop-------------------#
  loop_end_time <- Sys.time()
  last_loop_time <- as.numeric(difftime(loop_end_time, loop_start_time, units = "mins"))
  last_loop_time_formatted <- sprintf("%.2f mins", last_loop_time)
  leo_message(paste0("\n# Last Loop <", j, "> Time: ", last_loop_time_formatted, " #"))
}
res_stage1 <- do.call("rbind", res_stage1s); dim(res_stage1)
res_stage2 <- res_stage1 %>% filter(traits != "None"); dim(res_stage2)
output_path <- "/Users/leoarrow/project/iridocyclitis/output/hyprcoloc/473_res_stage1_0621.tsv"
fwrite(res_stage2, file = output_path, sep = "\t")
    
# 412 ----
gmb_files <- list.files(path = "/Volumes/T9/mbqtl/412DutchMicrobiomeProject/h.hg19", pattern = "\\.f.tsv.gz$", full.names = T)
gwas_list <- list(T1D = T1D, Iridocyclitis = Iridocyclitis); rm(list = c("T1D", "Iridocyclitis")); gc()
#----------------------------------main----------------------------------------#
res <- perform_hypr(gmb_files, gwas_list, prefix = "412_", gsub_id_pattern = ".f.tsv.gz",
                    save.plot.dir = "/Users/leoarrow/project/iridocyclitis/figure/hyprcoloc",
                    save.output.dir = "/Users/leoarrow/project/iridocyclitis/output/hyprcoloc",
                    select = c("variant_id", "chromosome", "base_pair_location", "beta", "standard_error", "p_value"))

# 211 ----
source("./code/1.2.leo_coloc_packages.R")
#' class
# gwas_list <- list(T1D = T1D, Iridocyclitis = Iridocyclitis); rm(list = c("T1D", "Iridocyclitis")); gc()
gmb_files <- list.files(path = "/Volumes/T9/mbqtl/211MiBioGen/MiBioGen_QmbQTL_summary_class", pattern = "\\.txt.gz$", full.names = T)
system(paste("gunzip -c ", gmb_files[1], "| head"))
perform_hypr(gmb_files, gwas_list, prefix = "211_class_",
             save.plot.dir = "/Users/leoarrow/project/iridocyclitis/figure/hyprcoloc",
             save.output.dir = "/Users/leoarrow/project/iridocyclitis/output/hyprcoloc",
             select = c("rsID", "chr", "bp", "beta", "SE", "P.weightedSumZ"))
#' family
gmb_files <- list.files(path = "/Volumes/T9/mbqtl/211MiBioGen/MiBioGen_QmbQTL_summary_family", pattern = "\\.txt.gz$", full.names = T)
system(paste("gunzip -c ", gmb_files[1], "| head"))
perform_hypr(gmb_files, gwas_list, prefix = "211_family_",
             save.plot.dir = "/Users/leoarrow/project/iridocyclitis/figure/hyprcoloc",
             save.output.dir = "/Users/leoarrow/project/iridocyclitis/output/hyprcoloc",
             select = c("rsID", "chr", "bp", "beta", "SE", "P.weightedSumZ"))
#' genus # This is the most specific catalog here for 211
gmb_files <- list.files(path = "/Volumes/T9/mbqtl/211MiBioGen/MiBioGen_QmbQTL_summary_genus", pattern = "\\.txt.gz$", full.names = T)
system(paste("gunzip -c ", gmb_files[1], "| head"))
perform_hypr(gmb_files, gwas_list, prefix = "211_genus_",
             save.plot.dir = "/Users/leoarrow/project/iridocyclitis/figure/hyprcoloc",
             save.output.dir = "/Users/leoarrow/project/iridocyclitis/output/hyprcoloc",
             select = c("rsID", "chr", "bp", "beta", "SE", "P.weightedSumZ"))
#' order
gmb_files <- list.files(path = "/Volumes/T9/mbqtl/211MiBioGen/MiBioGen_QmbQTL_summary_order", pattern = "\\.txt.gz$", full.names = T)
perform_hypr(gmb_files, gwas_list, prefix = "211_order_",
             save.plot.dir = "/Users/leoarrow/project/iridocyclitis/figure/hyprcoloc",
             save.output.dir = "/Users/leoarrow/project/iridocyclitis/output/hyprcoloc",
             select = c("rsID", "chr", "bp", "beta", "SE", "P.weightedSumZ"))
#' phylum
gmb_files <- list.files(path = "/Volumes/T9/mbqtl/211MiBioGen/MiBioGen_QmbQTL_summary_phylum", pattern = "\\.txt.gz$", full.names = T)
perform_hypr(gmb_files, gwas_list, prefix = "211_phylum_",
             save.plot.dir = "/Users/leoarrow/project/iridocyclitis/figure/hyprcoloc",
             save.output.dir = "/Users/leoarrow/project/iridocyclitis/output/hyprcoloc",
             select = c("rsID", "chr", "bp", "beta", "SE", "P.weightedSumZ"))

# 430 ----
gmb_files <- list.files(path = "/Volumes/T9/mbqtl/430German/liftover_hg19", pattern = "\\.tsv.gz$", full.names = T)
perform_hypr(gmb_files, gwas_list, prefix = "430_",
             save.plot.dir = "/Users/leoarrow/project/iridocyclitis/figure/hyprcoloc",
             save.output.dir = "/Users/leoarrow/project/iridocyclitis/output/hyprcoloc",
             select = c("SNP", "CHR", "BP", "BETA", "SE", "P"))
##
system(paste("gunzip -c ", gmb_files[148], "| head"), intern = T) %>% fread(text = ., sep="\t")