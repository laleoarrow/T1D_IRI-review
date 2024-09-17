#' @title `tsmr` pipeline is designed for diabetes and iridocyclitis causality inference study
#' @author `Lu Ao` (Chongqing, China)
#' @note 
#' $ We conducted the following tsmr analysis
#' $ There are 4 @1-4 exposure and 3 outcome @5-7
#'  - For T1D: 
#'    * `T1D (Chiou et al.)` VS `Iridocyclitis (FinnGen)`; @1 vs @5
#'    * `T1D (Chiou et al.)` VS `Iridocyclitis (UKB)`; @1 vs @6
#'    * `T1D (Chiou et al.)` VS `Iridocyclitis (FinnGen+UKB)`; @1 vs @7
#'    * `T1D (FinnGen)` VS `Iridocyclitis (UKB)`; @2 vs @6
#'  - For T2D:
#'    * `T2D (DIAMANTE)` VS `Iridocyclitis (FinnGen)`; @3 vs @5
#'    * `T2D (DIAMANTE)` VS `Iridocyclitis (UKB)`; @3 vs @6
#'    * `T2D (DIAMANTE)` VS `Iridocyclitis (FinnGen+UKB)`; @3 vs @7
#'    * `T2D (FinnGen)` VS `Iridocyclitis (UKB)`. @4 vs @6
#' $ All data in hg38 were lifted hg19 before MR analysis.
#' $ IV with MAF below 0.01 were excluded
#' 
#' >>> Load packages -----------------------------------------------------------
pacman::p_load(vroom, data.table, tidyverse, TwoSampleMR, ComplexHeatmap,
               MendelianRandomization, ieugwasr, circlize,
               ggplot2, forestploter, RadialMR)
source("./code/1.0.tsmr_packages.R")
# OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJsdWFvQHN0dS5jcW11LmVkdS5jbiIsImlhdCI6MTcxNjQzMjY3NCwiZXhwIjoxNzE3NjQyMjc0fQ.GoJLOrMQEshAcRwTutI41romdmHOaCX3iKDF1gsO_DLgxdte9_0ki7tXNoGliYTb2F9YOP5j9hG4O6LV5fEzPeYKyCss6vVGTECIWUsFODrmCY2U1S1xdUWQj3dJj9UqUU66qK_oEkL8v2zeUdAsR9blgwCnUYLTTXTnSMobG27FfNQEs4gqhpS-ovet9itqk5O1sMob2DcOGx-Na9VUOOgJrzk_Q6gQAfij-mWDQrXOGtTk0UQe4G2ElhXBn7ZSrWJAX721xZQ_862cfEO1nWGVZspVkNAbr2dT8e9AhPyI52djg5uqBSpJtId0to6Eo3CvElxTDN00NNEQUe0ofg"
# Sys.setenv(OPENGWAS_JWT=OPENGWAS_JWT)
# ieugwasr::user()

#' >>> Date import and filter --------------------------------------------------
#' We have `T1D (Chiou et al.)` and `T1D (FinnGen)` @12
#' We have `T2D (DIAMANTE)` and `T2D (FinnGen)` @34
#' We have `Iridocyclitis (FinnGen)` and `Iridocyclitis (UKB)` and `Iridocyclitis (FinnGen+UKB)` @567
t1d_p1 <- "~/project/iridocyclitis/data/diabete/1/GCST90014023_buildhg19.tsv"
t1d_p2 <- "~/project/iridocyclitis/data/diabete/1_finn/summary_stats_finngen_R10_T1D_hg19.tsv"
t2d_p1 <- "~/project/iridocyclitis/data/diabete/2/DIAMANTE-EUR.sumstat.txt"
t2d_p2 <- "~/project/iridocyclitis/data/diabete/2_finn/finngen_R10_T2D_hg19.tsv"
iri_p1 <- "~/project/iridocyclitis/data/iri_finn/r10/summary_stats_finngen_R10_H7_IRIDOCYCLITIS_hg19.tsv"
iri_p2 <- "~/project/iridocyclitis/data/iri_finn/ukb_finn/iri_fbmeta_hg19.tsv"
iri_p3 <- iri_p2
system(paste0("head -n 3 ", iri_p3))

# Exposure ---------------------------------------------------------------------
# T1D
t1d_1 <-  vroom(t1d_p1) %>% # `T1D (Chiou et al.)` # infer Neff as ref (https://github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-Sample-Size-and-Preparing-GWAS-Summary-Statistics)
 select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF, Samplesize) %>% 
 set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF", "N") %>% # N for observed sample size
 mutate(Phenotype = "T1D (Chiou et al.)",
        MAF = ifelse(EAF<0.5, EAF, 1-EAF),
        Neff = 4/((2*MAF*(1-MAF))*SE^2), # inferred Neff
        v = 18942 / (18942+501638), # Case_ratio
        Max_Neff = 4*v*(1-v)*(18942+501638)# ie, the possible max Neff
        ) %>%
  mutate(Neff = ifelse(Neff > 1.1*Max_Neff, 1.1*Max_Neff, Neff)) %>% 
  mutate(Neff = ifelse(Neff < 0.5*Max_Neff, 0.5*Max_Neff, Neff)) %>% 
  select(SNP, CHR, POS, A1, A2, P, BETA, SE, EAF, Phenotype, N, Neff)
t1d_2 <-  vroom(t1d_p2) %>% # `T1D (FinnGen)`
 select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF) %>% 
 set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF") %>%
 mutate(Phenotype = "T1D (FinnGen)",
        N = 4320+335112,
        Neff = TwoSampleMR::effective_n(4320, 335112),
        )
# T2D
t2d_1 <- vroom(t2d_p1) %>% # `T2D (DIAMANTE)`
 select(rsID, `chromosome(b37)`, `position(b37)`, effect_allele, other_allele, `Fixed-effects_p-value`, `Fixed-effects_beta`, `Fixed-effects_SE`, effect_allele_frequency ) %>% 
 set_names("SNP", "CHR", "POS", "A1", "A2","P", "BETA", "SE", "EAF") %>% 
 mutate(A1 = toupper(A1), A2 = toupper(A2),
        Phenotype = "T2D (DIAMANTE)", 
        N = 80154+853816,
        Neff = TwoSampleMR::effective_n(80154, 853816))
t2d_2 <- vroom(t2d_p2) %>% # `T2D (FinnGen)`
 select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF) %>% 
 set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF") %>%
 mutate(Phenotype = "T2D (FinnGen)",
        N = 65085+335112,
        Neff = TwoSampleMR::effective_n(65085, 335112))
# Iri
iri_1 <- vroom(iri_p1) %>% # `Iridocyclitis (FinnGen)`
 select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF) %>% 
 set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF") %>%
 mutate(Phenotype = "Iridocyclitis (FinnGen)", 
        N = 8016+390647,
        Neff = TwoSampleMR::effective_n(8016, 390647))
iri_2 <- vroom(iri_p2) %>% # `Iridocyclitis (UKB)`
 drop_na(UKBB_beta) %>% 
 select(SNP, CHR, BP, ALT, REF, UKBB_pval, UKBB_beta, UKBB_sebeta, UKBB_af_alt) %>% 
 set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF") %>%
 mutate(Phenotype = "Iridocyclitis (UKB)",
        N = 2304+413959,
        Neff = TwoSampleMR::effective_n(2304, 413959))
iri_3 <- vroom(iri_p3) %>% # `Iridocyclitis (FinnGen+UKB)`
 mutate(N_Cohort = case_when(is.na(FINNGEN_beta) & !is.na(UKBB_beta) ~ "UKB",
                             !is.na(FINNGEN_beta) & is.na(UKBB_beta) ~ "Finn",
                             !is.na(FINNGEN_beta) & !is.na(UKBB_beta) ~ "Finn+UKB",
                             TRUE ~ "Unknown"),
        N = case_when(N_Cohort == "UKB" ~ 2304+413959,
                      N_Cohort == "Finn" ~ 8016+390647,
                      N_Cohort == "Finn+UKB" ~ 2304+413959+8016+390647,
                      TRUE ~ NA),
        Neff =  case_when(N_Cohort == "UKB" ~ TwoSampleMR::effective_n(2304, 413959),
                          N_Cohort == "Finn" ~ TwoSampleMR::effective_n(8016, 390647),
                          N_Cohort == "Finn+UKB" ~ 4582.495+15709.64,
                          TRUE ~ NA),
        EAF = case_when(N_Cohort == "UKB" ~ UKBB_af_alt,
                        N_Cohort == "Finn" ~ FINNGEN_af_alt,
                        N_Cohort == "Finn+UKB" ~ (FINNGEN_af_alt*(8016+390647)+UKBB_af_alt*(2304+413959))/(8016+390647+2304+413959),
                        TRUE ~ NA)) %>% 
 select(SNP, CHR, BP, ALT, REF, all_inv_var_meta_p, all_inv_var_meta_beta, all_inv_var_meta_sebeta, EAF, N, Neff) %>% 
 set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF", "N", "Neff") %>%
 mutate(Phenotype = "Iridocyclitis (FinnGen+UKB)")

# extract instruments
p_cutoff <- 5e-8
iv1 <- extract_instruments_local(t1d_1, p = p_cutoff, N = "Neff", bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR", plink_bin = plinkbinr::get_plink_exe())
iv2 <- extract_instruments_local(t1d_2, p = p_cutoff, N = "Neff", bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR", plink_bin = plinkbinr::get_plink_exe())
iv3 <- extract_instruments_local(t2d_1, p = p_cutoff, N = "Neff", bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR", plink_bin = plinkbinr::get_plink_exe())
iv4 <- extract_instruments_local(t2d_2, p = p_cutoff, N = "Neff", bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR", plink_bin = plinkbinr::get_plink_exe())
iv1$chr.exposure <- as.character(iv1$chr.exposure)
iv2$chr.exposure <- as.character(iv2$chr.exposure)
iv3$chr.exposure <- as.character(iv3$chr.exposure)
iv4$chr.exposure <- as.character(iv4$chr.exposure)
iv <- bind_rows(iv1, iv2, iv3, iv4)
rm(iv1,iv2,iv3,iv4);gc()

# R2 & F statistics
# https://mp.weixin.qq.com/s/NYjGglQqQRVDD8O4hnOBwg
# https://www.wolai.com/sharespace/i2dArNtbTo3ZQDswDZRmVq
iv <- iv %>% 
  mutate(R2 = (2*(beta.exposure)^2*(eaf.exposure*(1-eaf.exposure))) / ((se.exposure^2) * samplesize.exposure), # R2 = [2*beta^2*maf*(1-maf)]/[se^2*N]
         F = R2*(samplesize.exposure-2)/(1-R2) # F = (N-2)*(R2/(1-R2)) # i.e., K = 1 for one SNP
         ) %>% 
  filter(F > 10)
# iv maf
iv$maf.exposure <- ifelse(iv$eaf.exposure < 0.5, iv$eaf.exposure, 1-iv$eaf.exposure)
table(iv$maf.exposure > 0.01) # should be all TRUE
iv <- iv %>% subset(iv$maf.exposure > 0.01)
iv$id.exposure <- iv$exposure
iv %>% vroom_write("./output/tsmr/tsmr_iv_date.tsv")

# Outcome ----------------------------------------------------------------------
out_1 <- format_outcome(iri_1, iv$SNP, N = "Neff")
out_2 <- format_outcome(iri_2, iv$SNP, N = "Neff")
out_3 <- format_outcome(iri_3, iv$SNP, N = "Neff")
out_1 <- out_1 %>% mutate(maf.outcome = ifelse(eaf.outcome<0.5, eaf.outcome, 1-eaf.outcome)) %>% subset(maf.outcome>0.01)
out_2 <- out_2 %>% mutate(maf.outcome = ifelse(eaf.outcome<0.5, eaf.outcome, 1-eaf.outcome)) %>% subset(maf.outcome>0.01)
out_3 <- out_3 %>% mutate(maf.outcome = ifelse(eaf.outcome<0.5, eaf.outcome, 1-eaf.outcome)) %>% subset(maf.outcome>0.01)
out = rbind(out_1,out_2,out_3)
rm(out_1,out_2,out_3);gc()

# Harmonise/MR----------------------------------------------------------------
dat_h <- harmonise_data(iv, out) %>% subset(mr_keep)
dat_h$id.exposure <- dat_h$exposure; dat_h$id.outcome <- dat_h$outcome 
dat_h %>% vroom_write("./output/tsmr/tsmr_dat_h_date.tsv")

# unique exposure and outcome pair
outcome_order <- c("Iridocyclitis (FinnGen)", "Iridocyclitis (UKB)", "Iridocyclitis (FinnGen+UKB)")
uni_pair <- dat_h %>% 
  filter(!(grepl("FinnGen", exposure) & grepl("FinnGen", outcome))) %>% 
  select(exposure, outcome) %>% 
  unique() %>% 
  arrange(exposure, factor(outcome, levels = outcome_order))
rownames(uni_pair) <- 1:nrow(uni_pair)

# >>>>>>>>>>>>> formal tsmr analysis # >>>>>>>>>>>>> 
library(parallel); library(pbapply)
nc <- 14
cl <- makeCluster(nc, type = "FORK")
results_list <- pblapply(1:nrow(uni_pair), function(i){
  exp <- uni_pair[i,"exposure"]
  out <- uni_pair[i,"outcome"]
  res <- mr_one_pair(dat_h, exp, out)
  return(res)
  }, cl=cl)
stopCluster(cl);rm(nc);rm(cl);gc()
# The `results_list` is actually a matrix becase the fucntion return a list of two data.frames
ress <- lapply(results_list, function(x) x$res_pair) %>% do.call(rbind, .)
ress_presso <- lapply(results_list, function(x) x$res_pair_presso) %>% do.call(rbind, .)
vroom_write(ress %>% filter(grepl("T1D*", exposure)), "./output/tsmr/T1D_0601.txt")
vroom_write(ress %>% filter(grepl("T2D*", exposure)), "./output/tsmr/T2D_0601.txt")
vroom_write(ress_presso %>% filter(grepl("T1D*", Exposure)), "./output/tsmr/T1D_PRESSO_0601.txt")
vroom_write(ress_presso %>% filter(grepl("T2D*", Exposure)), "./output/tsmr/T2D_PRESSO_0601.txt")

# MR-lap -----------------------------------------------------------------------
library(MRlap) # NOTE that MRlap need observed N
.clean_names <- function(df) { df %>%
    mutate(across(c(exposure, outcome), ~ gsub("[(]", "_", .))) %>%
    mutate(across(c(exposure, outcome), ~ gsub("[)]", "", .))) %>%
    mutate(across(c(exposure, outcome), ~ gsub(" ", "_", .))) %>%
    mutate(across(c(exposure, outcome), ~ gsub("__", "_", .))) %>%
    mutate(across(c(exposure, outcome), ~ gsub("_et_al.", "", .)))
}
t1d_pair <- fread("/Users/leoarrow/project/iridocyclitis/output/tsmr/T1D_0531.txt") %>% 
  select(exposure, outcome) %>% unique() %>% .clean_names() # generated from tsmr pipline
t2d_pair <- fread("/Users/leoarrow/project/iridocyclitis/output/tsmr/T2D_0531.txt") %>% 
  select(exposure, outcome) %>% unique() %>% .clean_names() # generated from tsmr pipline
uni_pair <- rbind(t1d_pair, t2d_pair) %>% 
  mutate(
    element_exp = case_when(
      exposure == "T1D_Chiou" ~ "t1d_1",
      exposure == "T1D_FinnGen" ~ "t1d_2",
      exposure == "T2D_DIAMANTE" ~ "t2d_1",
      exposure == "T2D_FinnGen" ~ "t2d_2",
      TRUE ~ NA_character_
    ),
    element_out = case_when(
      outcome == "Iridocyclitis_FinnGen" ~ "iri_1",
      outcome == "Iridocyclitis_UKB" ~ "iri_2",
      outcome == "Iridocyclitis_FinnGen+UKB" ~ "iri_3",
      TRUE ~ NA_character_
    )
  )
# Note: use N (observed) for MR-lap and LDSC
# >>>>>>>>>>>>> formal mr-lap analysis # >>>>>>>>>>>>> 
library(pbapply) # not cool to use multi-core here
mrlap_res_list <- pblapply(1:nrow(uni_pair), function(i) {
  mrlap_one_pair(exposure_data = get(uni_pair$element_exp[i]),
                 exposure_name = uni_pair$exposure[i],
                 outcome_data = get(uni_pair$element_out[i]),
                 outcome_name = uni_pair$outcome[i],
                 ld_path = "/Users/leoarrow/project/software/ldsc-b站/eur_w_ld_chr_nomhc",
                 hm3_path = "/Users/leoarrow/project/software/ldsc-b站/w_hm3.noMHC.snplist",
                 log = T, log_path = "/Users/leoarrow/project/iridocyclitis/output/tsmr/MRlap"
                 )
}, cl = NULL)
mrlap_res <- do.call(rbind, mrlap_res_list)
vroom_write(mrlap_res, "./output/tsmr/mrlap_0601.tsv")

# others------------------------------------------------------------------------
# sample overlap estimation: 
# https://sb452.shinyapps.io/overlap/
# https://mp.weixin.qq.com/s?__biz=MzIzODExNDk2MQ==&mid=2648747922&idx=1&sn=6e865d1f555d2dbc278ab37bbddf627c&chksm=f12a4e6cc65dc77ae38f7a84572de544f58c86005cd06c4a68b18b7aa829c65591975e0263e9&scene=21#wechat_redirect
# https://mp.weixin.qq.com/s?__biz=MzIzODExNDk2MQ==&mid=2648750119&idx=1&sn=5ece7d926421e714a8fca8ed425d8078&chksm=f0253001e16748c5fd7f4e5be777e8320354e2e7eb1e552aa0b93c15792d46fda9abefb2edd3&scene=132&exptype=timeline_recommend_article_extendread_samebiz&show_related_article=1&subscene=27&scene=132#wechat_redirect

# Figure 1 visualization -------------------------------------------------------
# Aims: To use ComplexHeatmap to present all MR results
# rm(list=ls());gc()
# pre-process data for visualization
res <- fread("/Users/leoarrow/project/iridocyclitis/output/tsmr/tsmr_figure.tsv") %>% 
  mutate(pair = paste0(exposure, " VS ", outcome)) %>%
  mutate(method = gsub("Inverse variance weighted", "IVW", method)) %>% 
  mutate(`N (exp)` = case_when(grepl("Chiou", exposure) ~ 18942+501638,
                               grepl("T1D \\(FinnGen\\)", exposure) ~ 4320+335112,
                               grepl("DIAMANTE", exposure) ~ 80154+853816,
                               grepl("T2D \\(FinnGen\\)", exposure) ~ 65085+335112),
         `N (out)` = case_when(grepl("Iridocyclitis \\(FinnGen\\)", outcome) ~ 8016+390647,
                               grepl("Iridocyclitis \\(UKB\\)", outcome) ~ 2304+413959,
                               grepl("FinnGen\\+UKB", outcome) ~ 2304+413959+8016+390647),
         significance = case_when(pval < 0.001 ~ "***",
                                  pval < 0.01 ~ "**",
                                  pval < 0.05 ~ "*",
                                  pval < 0.1 ~ ".",
                                  TRUE ~ ""),
         lower_ci = b - 1.96 * se,
         upper_ci = b + 1.96 * se)
# res %>% select(pair, N_exp, N_out) %>% unique() %>% View() # check
long_to_wide <- function(data, method_col, pair_col, value_col) {
  wide_data <- data %>%
    select({{method_col}}, {{pair_col}}, {{value_col}}) %>%
    pivot_wider(names_from = {{pair_col}}, values_from = {{value_col}})
  return(wide_data)
}
# change to wide data
betas <- long_to_wide(res, method, pair, b)
significance <- long_to_wide(res, method, pair, significance)
lower_ci <- long_to_wide(res, method, pair, lower_ci)
upper_ci <- long_to_wide(res, method, pair, upper_ci)
# as matrix
betas_matrix <- as.matrix(betas[,-1])
significance_matrix <- as.matrix(significance[,-1])
lower_ci_matrix <- as.matrix(lower_ci[,-1])
upper_ci_matrix <- as.matrix(upper_ci[,-1])
# modify colname
multi_colname <- function(str) {
  str <- gsub(" VS ", " \nOut: ", str)
  str <- gsub("(.+)", "Exp: \\1", str)
  return(str)
}
colnames(betas_matrix) <- colnames(betas_matrix) %>% multi_colname
rownames(betas_matrix) <- betas$method
colnames(significance_matrix) <- colnames(significance_matrix) %>% multi_colname
rownames(significance_matrix) <- betas$method
colnames(lower_ci_matrix) <- colnames(lower_ci_matrix) %>% multi_colname
rownames(lower_ci_matrix) <- betas$method
colnames(upper_ci_matrix) <- colnames(upper_ci_matrix) %>% multi_colname
rownames(upper_ci_matrix) <- betas$method

# 🔥
# list_components()
# seekViewport("annotation_N_2")
library(ComplexHeatmap)
col_fun <- circlize::colorRamp2(
  c(min(betas_matrix),0, max(betas_matrix)),
  c("#4DBBD5FF", "white", "#E64B35FF")
)

N_matrix <- matrix(c(520580,520580,520580,339432,933970,933970,933970,400197,
                     398663,416263,814926,416263,398663,416263,814926,416263), ncol = 2, byrow = F) # 第一列是Exposure； 第二列是Outcome
column_anno <- HeatmapAnnotation(Type = anno_block(gp = gpar(fill = c("#E64B35FF","#4DBBD5FF"), col = "black",alpha = 0.3), # annotate block
                                                   height = unit(0.5,"cm"),
                                                   labels = c("T1D", "T2D"),
                                                   labels_gp = gpar(col = "black", fontsize = 10)),
                                 N = anno_points(N_matrix,
                                                 pch = c(17,19),
                                                 size = unit(2.5, "mm"), # 设置点的大小
                                                 ylim = c(min(N_matrix), max(N_matrix)),
                                                 extend = 0.2,
                                                 gp = gpar(fill = T, alpha = 0.9, col = c(2,"#00A087FF")),
                                                 axis_param = list(
                                                   at = c(400000, 600000, 800000, 1000000), 
                                                   labels = c("400K", "600K", "800K", "1M"))
                                                 ),
                                 border = "black"
)
X <- Heatmap(betas_matrix, cluster_rows = FALSE, cluster_columns = FALSE,
             name = "Effect\nEstimate", # betas_matrix的value是啥
             col = col_fun,
             border = "black", # 小方块border
             rect_gp = gpar(col = "black", lwd = 1), # 小方块border
             heatmap_legend_param = list(legend_height = unit(3.5, "cm")),
             #' >>>>>>>>> annotation >>>>>>>>> 
             top_annotation = column_anno,
             #' >>>>>>>>> cell_fun >>>>>>>>> 
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(paste0(significance_matrix[i, j]), # 显著性标记
                         x, y, gp = gpar(fontsize = 20, fontface = "bold"))
               grid.text(paste0("\n[", round(lower_ci_matrix[i, j], 2), ", ", round(upper_ci_matrix[i, j], 2), "]"), 
                         x, y, gp = gpar(fontsize = 8))  # 置信区间
             },
             #' >>>>>>>>> split >>>>>>>>> 
             column_split = rep(c("T1D", "T2D"), each = 4),
             column_title = NULL,
             column_title_side = "top",
             column_title_gp = gpar(fontsize = 10, fill = "white", just = c("center","center")),
             
             # >>>>>>>>> row col index
             row_names_gp = gpar(fontsize = 10),
             row_names_side = "left",
             column_names_rot = 90, # rotation
             column_names_gp = gpar(fontsize = 10),
             column_names_side = "bottom"
)
pdf("./figure/tsmr/heatmap.pdf", width = 9, height = 5)
draw(X, annotation_legend_list = list(Legend(
  title = "Sample\nSize", type = "points", 
  labels = c("Exposure","Outcome"), pch = c(17,19),
  legend_gp = gpar(col = 2:3)
)))
dev.off()
# legend: https://blog.csdn.net/dxs18459111694/article/details/134520255

#' Abandoned code
# add Neff(max) for each pair
# dat_h <- fread("./output/tsmr/tsmr_dat_h_0601.tsv")
# dat_h %>% select(id.exposure, id.outcome) %>% # check
#   unique() %>% 
#   filter(!(grepl("Finn", id.exposure) & grepl("Finn", id.outcome)))
# get_Neff <- function(dat_h, exp, out, col, type) {
#   exp <- as.character(exp)
#   out <- as.character(out)
#   col <- as.character(col)
#   if (!type %in% c("max", "mean")) {
#     stop("Must give a valid type: 'max' or 'mean'!")
#   }
#   filtered_data <- dat_h %>% filter(exposure == exp, outcome == out) %>% pull(col)
#   result <- case_when(
#     type == "max" ~ max(filtered_data, na.rm = TRUE),
#     type == "mean" ~ mean(filtered_data, na.rm = TRUE)
#   )
#   return(result)
# }
# for (i in 1:nrow(res)) {
#   res[i, "Neff_exp_max"] <- get_Neff(dat_h, res[i, "exposure"], res[i, "outcome"], "samplesize.exposure", type = "max")
#   res[i, "Neff_exp_mean"] <- get_Neff(dat_h, res[i, "exposure"], res[i, "outcome"], "samplesize.exposure", type = "mean")
#   res[i, "Neff_out_max"] <- get_Neff(dat_h, res[i, "exposure"], res[i, "outcome"], "samplesize.outcome", type = "max")
#   res[i, "Neff_out_mean"] <- get_Neff(dat_h, res[i, "exposure"], res[i, "outcome"], "samplesize.outcome", type = "mean")
# }


# In reverse >>>> -------------------------------------------------------------------
# extract instruments
p_cutoff <- 5e-8
iv1 <- extract_instruments_local(iri_1, p = p_cutoff, N = "Neff", bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR", plink_bin = plinkbinr::get_plink_exe())
iv2 <- extract_instruments_local(iri_2, p = p_cutoff, N = "Neff", bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR", plink_bin = plinkbinr::get_plink_exe())
iv3 <- extract_instruments_local(iri_3, p = p_cutoff, N = "Neff", bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR", plink_bin = plinkbinr::get_plink_exe())
iv1$chr.exposure <- as.character(iv1$chr.exposure)
iv2$chr.exposure <- as.character(iv2$chr.exposure)
iv3$chr.exposure <- as.character(iv3$chr.exposure)
iv <- bind_rows(iv1, iv2, iv3)
rm(iv1,iv2,iv3);gc()

# R2 & F statistics
# https://mp.weixin.qq.com/s/NYjGglQqQRVDD8O4hnOBwg
# https://www.wolai.com/sharespace/i2dArNtbTo3ZQDswDZRmVq
iv <- iv %>% 
  mutate(R2 = (2*(beta.exposure)^2*(eaf.exposure*(1-eaf.exposure))) / ((se.exposure^2) * samplesize.exposure), # R2 = [2*beta^2*maf*(1-maf)]/[se^2*N]
         F = R2*(samplesize.exposure-2)/(1-R2) # F = (N-2)*(R2/(1-R2)) # i.e., K = 1 for one SNP
  ) %>% 
  filter(F > 10)
# iv maf
iv$maf.exposure <- ifelse(iv$eaf.exposure < 0.5, iv$eaf.exposure, 1-iv$eaf.exposure)
table(iv$maf.exposure > 0.01) # should be all TRUE
iv <- iv %>% subset(iv$maf.exposure > 0.01)
iv$id.exposure <- iv$exposure
iv %>% vroom_write("./output/tsmr/tsmr_iv_date.tsv")

# Outcome ----------------------------------------------------------------------
out_1 <- format_outcome(t1d_1, iv$SNP, N = "Neff")
out_2 <- format_outcome(t1d_2, iv$SNP, N = "Neff")
out_3 <- format_outcome(t2d_1, iv$SNP, N = "Neff")
out_4 <- format_outcome(t2d_2, iv$SNP, N = "Neff")
out_1 <- out_1 %>% mutate(maf.outcome = ifelse(eaf.outcome<0.5, eaf.outcome, 1-eaf.outcome)) %>% subset(maf.outcome>0.01)
out_2 <- out_2 %>% mutate(maf.outcome = ifelse(eaf.outcome<0.5, eaf.outcome, 1-eaf.outcome)) %>% subset(maf.outcome>0.01)
out_3 <- out_3 %>% mutate(maf.outcome = ifelse(eaf.outcome<0.5, eaf.outcome, 1-eaf.outcome)) %>% subset(maf.outcome>0.01)
out_4 <- out_4 %>% mutate(maf.outcome = ifelse(eaf.outcome<0.5, eaf.outcome, 1-eaf.outcome)) %>% subset(maf.outcome>0.01)
out = rbind(out_1,out_2,out_3,out_4)
rm(out_1,out_2,out_3,out_4);gc()

# Harmonise/MR----------------------------------------------------------------
dat_h <- harmonise_data(iv, out) %>% subset(mr_keep)
dat_h$id.exposure <- dat_h$exposure; dat_h$id.outcome <- dat_h$outcome 
table(dat_h$id.exposure); table(dat_h$id.outcome)
# dat_h %>% vroom_write("./output/tsmr/rev_tsmr_dat_h_date.tsv")

# unique exposure and outcome pair
exporsure_order <- c("Iridocyclitis (FinnGen)", "Iridocyclitis (UKB)", "Iridocyclitis (FinnGen+UKB)")
# outcome_order <- c("T1D (Chiou et al.)", "T1D (FinnGen)", "T2D (DIAMANTE)", "T2D (FinnGen)")
uni_pair <- dat_h %>% 
  filter(!(grepl("FinnGen", exposure) & grepl("FinnGen", outcome))) %>% 
  select(exposure, outcome) %>% 
  unique() %>% 
  arrange(outcome, factor(exposure, levels = exporsure_order))
rownames(uni_pair) <- 1:nrow(uni_pair)

x <- mr(uni_pair)

# >>>>>>>>>>>>> formal tsmr analysis # >>>>>>>>>>>>> 
library(parallel); library(pbapply)
nc <- 8
cl <- makeCluster(nc, type = "FORK")
results_list <- pblapply(1:nrow(uni_pair), function(i){
  exp <- uni_pair[i,"exposure"]
  out <- uni_pair[i,"outcome"]
  res <- mr_one_pair(dat_h, exp, out, save_plot = F, res_dir= "./output/tsmr/reverse", fig_dir="./figure/tsmr/reverse")
  return(res)
}, cl=cl)
stopCluster(cl);rm(nc);rm(cl);gc()
# The `results_list` is actually a matrix becase the fucntion return a list of two data.frames
ress <- lapply(results_list, function(x) x$res_pair) %>% do.call(rbind, .)
ress_presso <- lapply(results_list, function(x) x$res_pair_presso) %>% do.call(rbind, .)
vroom_write(ress, "./output/tsmr/reverse/res_all.txt")
vroom_write(ress_presso, "./output/tsmr/reverse/res_presso_all.txt")

