# This file contains both R and Bash code used for coloc analysis.
# Ao Lu, Chongqing, China
# rm(list = ls());gc()
# a=ls();rm(list = setdiff(a, c("gwas1", "gwas2")));gc()

# bin/R!------------------------------------------------------------
source("./code/1.2.leo_coloc_packages.R")
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(vroom))
suppressMessages(library(coloc))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(ieugwasr))
suppressMessages(library(plinkbinr))
suppressMessages(library(LDlinkR))

# >>>>>>>>>>>>>>> 2024 >>>>>>>>>>>>>>>>
# Please do include all SNPs you have in the region: do not trim by significance
# or MAF or other means (https://chr1swallace.github.io/coloc/articles/a02_data.html)
# No trim for MAF !!!

# * Read my gwas ----
# system(paste("head -n 2 ", gwas1_path))
gwas1_path <- "~/project/iridocyclitis/data/diabete/1/GCST90014023_buildhg19.tsv" # t1d1
gwas1 <- vroom(gwas1_path) %>% 
  dplyr::select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF, Samplesize) %>% 
  set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF", "N") %>% 
  mutate(Phenotype = "t1d1", VARBETA = SE^2, MAF = ifelse(EAF<0.5, EAF, 1-EAF)) %>% 
  mutate_at("POS", as.integer)

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
  mutate(Phenotype = "iri3", VARBETA = SE^2, MAF = ifelse(EAF<0.5, EAF, 1-EAF)) %>% 
  mutate_at("POS", as.integer)

# exclue MHC
gwas1 <- gwas1 %>% filter(!(CHR == 6 & POS >= 25000000 & POS <= 34000000))
gwas2 <- gwas2 %>% filter(!(CHR == 6 & POS >= 25000000 & POS <= 34000000))

# placo lead snp
snp_annotated <- data.table::fread("/Users/leoarrow/project/iridocyclitis/output/fuma/placo/FUMA_job493590/[snps].txt")
lead <- data.table::fread("/Users/leoarrow/project/iridocyclitis/output/fuma/placo/FUMA_job493590/【leadSNPs】.txt") %>% 
  dplyr::left_join(snp_annotated %>% dplyr::select(rsID, nearestGene), by = c("rsID" = "rsID")) %>% 
  dplyr::select(rsID, chr, pos, nearestGene) %>% 
  set_names("SNP", "CHR", "POS", "nearestGene")

# * Coloc parameter ----
gwas1 = gwas1
gwas2 = gwas2         
type1 <- "cc"         
type2 <- "cc"         
sample_size1 <- 520580
sample_size2 <- 814926
s1 <- 0.03638634      
s2 <- 0.01266373      
coloc_one_side <-  100

### coloc
# refer to `1.2.leo_coloc_packages.R`
source("./code/1.2.leo_coloc_packages.R")
coloc_lead_snps(gwas1, gwas2, lead_snp, window = 500, 
                out_fp = "./figure/coloc", 
                out_p = "./output/coloc")

# summary visualization ----
library(ggplot2)
library(ggchicklet)
library(dplyr)
res <- fread("./output/coloc/coloc_results.tsv") %>% select(gene, PP3, PP4) %>% unique()
res <- res %>%
  mutate(PP_all = PP3+PP4) %>% 
  arrange(desc(PP_all), desc(PP4)) %>%
  mutate(order = 1:nrow(.))

res_long <- res %>% 
  pivot_longer(cols = c(PP3, PP4), names_to = "PP", values_to = "value") 


p <- ggplot(res_long, aes(y = reorder(gene, -order), x = value, fill = PP)) +
  geom_col(color = NA, alpha = 0.8, width = 0.3) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = c("PP4" = "#E64B35FF", "PP3" = "#00A087FF")) +
  labs(x = "Posterior Probabilities", y = "Non-MHC Gene Loci", title = "Colocalization Evidence", fill = "PP Type") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0, size = 20, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black", size = .5),
    axis.ticks = element_line(color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.2),
    legend.title = element_text(size = 14, face = "bold"),
    legend.box.background = element_rect(color = "black", size = 1, fill = NA),
    legend.text = element_text(size = 12)
  ); p
ggsave("./figure/coloc/coloc_all.pdf", plot = p, width = 6, height = 6)