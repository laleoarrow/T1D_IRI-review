# Multivarible MR for glycemic traits on iridocyclitis
# >>> load packages & function
pacman::p_load(vroom, data.table, tidyverse, ggplot2, forestploter,
               TwoSampleMR, MendelianRandomization, MRPRESSO, MVMR, ieugwasr)
source("./code/1.0.tsmr_packages.R")

# https://api.opengwas.io/
OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJsdWFvQHN0dS5jcW11LmVkdS5jbiIsImlhdCI6MTcyMDg2OTUwNiwiZXhwIjoxNzIyMDc5MTA2fQ.QSU5blzZexib-rOjfBV9e8wSSRSv5LeisjcbhO9uJQ5osmUYnVj-VRxdniL07C6_b3c9wzwS0gmbp6zGLw6dV0TniqdL32E-mMc1CnF5jiKext9Jc3i5dYmaMqjKigkFjZNlevl2Ri_pqJpt4u4ushPITcuNKznhu5ZUcCp0WDxe6DXRqaxgRZgv7fVzujgSscxMt58bF6yXCgwZ-acUuqPkom_aq03tBVXXgdimPK7btVP8eF09OXBjaNXE_436PA57xyiBwaRuSFhIh39FxgbGjF1FKEAwrThXFvehOh3t_kbkOX6BW8Rv-uYJXbBftIkaTQGffbMOlqDJlo5BIA"
Sys.setenv(OPENGWAS_JWT=OPENGWAS_JWT)
ieugwasr::user()
system("curl cip.cc")

# exposure
gwas_id <- c("ebi-a-GCST90002232", # Fasting glucose
             "ebi-a-GCST90002244", # Glycated hemoglobin levels
             "ebi-a-GCST90002238", # Fasting insulin
             "ebi-a-GCST90002227") # Two-hour glucose
exposure_dat <- mv_extract_exposures(gwas_id)
# Clumping 1, 181 variants, using EUR population reference
# Removing 49 of 181 variants due to LD with other variants or absence from LD reference panel
# Extracting data for 132 SNP(s) from 4 GWAS(s)
# Finding proxies for 1 SNPs in outcome ebi-a-GCST90002232
# Extracting data for 1 SNP(s) from 1 GWAS(s)
# Finding proxies for 1 SNPs in outcome ebi-a-GCST90002238
# Extracting data for 1 SNP(s) from 1 GWAS(s)
# Harmonising Fasting glucose || id:ebi-a-GCST90002232 (ebi-a-GCST90002232) and Two-hour glucose || id:ebi-a-GCST90002227 (ebi-a-GCST90002227)
# Removing the following SNPs for being palindromic with intermediate allele frequencies:
#   rs10487796
# Harmonising Fasting glucose || id:ebi-a-GCST90002232 (ebi-a-GCST90002232) and Fasting insulin || id:ebi-a-GCST90002238 (ebi-a-GCST90002238)
# Removing the following SNPs for being palindromic with intermediate allele frequencies:
#   rs10487796
# Harmonising Fasting glucose || id:ebi-a-GCST90002232 (ebi-a-GCST90002232) and Glycated hemoglobin levels || id:ebi-a-GCST90002244 (ebi-a-GCST90002244)
# Removing the following SNPs for being palindromic with intermediate allele frequencies:
#   rs10487796
iv <- exposure_dat
# iv maf
iv$maf.exposure <- ifelse(iv$eaf.exposure < 0.5, iv$eaf.exposure, 1-iv$eaf.exposure)
table(iv$maf.exposure > 0.01) # should be all TRUE
iv <- iv %>% subset(iv$maf.exposure > 0.01)
iv$id.exposure <- iv$exposure
n_iv <- dim(iv)[1]/length(unique(iv$id.exposure)) # how many iv for the mvmr

# outcome
iri_p3 <- "~/project/iridocyclitis/data/iri_finn/ukb_finn/iri_fbmeta_hg19.tsv"
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
out_3 <- format_outcome(iri_3, iv$SNP, N = "Neff") # neff is useless for outcome data though

# harmonize data
mvdat <- mv_harmonise_data(iv, out_3)
beta <- mvdat$exposure_beta; se <- mvdat$exposure_se; out_b <- as.matrix(mvdat$outcome_beta); out_se <- as.matrix(mvdat$outcome_se)
colnames(beta) <- paste0("b_", colnames(beta)); colnames(se) <-  paste0("se_", colnames(se))
colnames(out_b) <- "out_b"; colnames(out_se) <- "out_se"
ivs <- cbind(beta, se, out_b, out_se) %>% as.data.frame() %>% mutate(SNP = rownames(.))
vroom::vroom_write(ivs, "./output/tsmr/bloodsugar/mvmr_iv_0716.tsv")
n_iv <- dim(mvdat$exposure_beta)[1] # how many iv for the mvmr

# >>>>>>>>>>>>>>> twosamplemr——mvmr >>>>>>>>>>>>>>> #
res_IVW_tsmr <- mv_multiple(mvdat, intercept = F)
res_IVW_tsmr_OR <- generate_odds_ratios(res_IVW_tsmr$result); res_IVW_tsmr_OR
res_IVW_tsmr_OR <- res_IVW_tsmr_OR %>% mutate(Method = "MVMR-IVW", N_iv = n_iv) %>% select(-id.exposure,-id.outcome)

res_RUR_tsmr <- mv_residual(mvdat, intercept = F) # ru for residualised unweighted regression (Burgess et al. 2015)
res_RUR_tsmr_OR <- generate_odds_ratios(res_RUR_tsmr$result); res_RUR_tsmr_OR
res_RUR_tsmr_OR <- res_RUR_tsmr_OR %>% mutate(Method = "MVMR-RUR", N_iv = n_iv) %>% select(-id.exposure,-id.outcome)

twosamplemr_res <- rbind(res_IVW_tsmr_OR, res_RUR_tsmr_OR)
twosamplemr_res %>% vroom::vroom_write("./output/tsmr/bloodsugar/mvmr_tsmr_IVW_RUR_0716.tsv")

# >>>>>>>>>>>>>>> presso——mvmr >>>>>>>>>>>>>>> #
presso_summary <- cbind(mvdat[["outcome_beta"]],
                        mvdat[["exposure_beta"]][,1],
                        mvdat[["exposure_beta"]][,2],
                        mvdat[["exposure_beta"]][,3],
                        mvdat[["exposure_beta"]][,4],
                        mvdat[["exposure_se"]][,1],
                        mvdat[["exposure_se"]][,2],
                        mvdat[["exposure_se"]][,3],
                        mvdat[["exposure_se"]][,4],
                        mvdat[["outcome_se"]]) %>% as.data.frame()

library(MRPRESSO)
res_mvmr_presso <- mr_presso(BetaOutcome = "V1",
                             BetaExposure = c("V2","V3","V4","V5"),
                             SdOutcome = "V10",
                             SdExposure = c("V6","V7", "V8", "V9"),
                             OUTLIERtest = TRUE,
                             DISTORTIONtest = TRUE,
                             data = presso_summary,
                             NbDistribution = ifelse(nrow(presso_summary)/0.05<1000, 1000, nrow(presso_summary)/0.05),
                             SignifThreshold = 0.05)
# global test
res_presso$`MR-PRESSO results` %>% as.data.frame() %>% vroom::vroom_write("./output/tsmr/bloodsugar/mvmr_presso_global_0716.tsv")
res_presso_result <- res_presso$`Main MR results` %>% 
  as.data.frame() %>% 
  mutate(Exposure =  rep(colnames(mvdat$exposure_beta),2)) %>% 
  dplyr::select(-`T-stat`) %>% 
  dplyr::rename(b = `Causal Estimate`, se = Sd, pval = `P-value`) %>%
  TwoSampleMR::generate_odds_ratios()
res_presso_result %>% vroom::vroom_write("./output/tsmr/bloodsugar/mvmr_presso_result_0716.tsv")

# >>>>>>>>>>>>>>> MendelianRandomization——mvmr >>>>>>>>>>>>>>> #
library(MendelianRandomization)
MRMVInputObject <- mr_mvinput(bx = cbind(mvdat[["exposure_beta"]][,1], mvdat[["exposure_beta"]][,2], mvdat[["exposure_beta"]][,3], mvdat[["exposure_beta"]][,4]),
                              bxse = cbind(mvdat[["exposure_se"]][,1], mvdat[["exposure_se"]][,2], mvdat[["exposure_se"]][,3], mvdat[["exposure_se"]][,4]),
                              by = mvdat[["outcome_beta"]],
                              byse = mvdat[["outcome_se"]])
res_IVW_MRObject <- mr_mvivw(MRMVInputObject, 
                             model = "default",
                             correl = FALSE,
                             distribution = "normal",
                             alpha = 0.05)
res_egger_MRObject <- mr_mvegger(MRMVInputObject,
                                 orientate = 1,
                                 correl = FALSE,
                                 distribution = "normal",
                                 alpha = 0.05)
res_median_MRObject <- mr_mvmedian(MRMVInputObject,
                                   distribution = "normal",
                                   alpha = 0.05,
                                   iterations = 10000,
                                   seed = 314159265)
get_result_from_MRObject <- function(MRObject, mvdat){
  exp <- mvdat$exposure_beta %>% colnames()
  out <- mvdat$outname$outcome
  res <- data.frame(
    exposure = exp,
    outcome = rep(out,length(exp)),
    b = MRObject@Estimate,
    se = tryCatch(MRObject@StdError, error = function(e) MRObject@StdError.Est),
    pval = tryCatch(MRObject@Pvalue, error = function(e) MRObject@Pvalue.Est),
    Method = rep(class(MRObject)[[1]],length(exp)),
    N_iv = MRObject@SNPs
  )
  res <- res %>% TwoSampleMR::generate_odds_ratios() %>% 
    dplyr::select(exposure, outcome, b, se, pval, lo_ci, up_ci, or, or_lci95, or_uci95, Method, N_iv)
  return(res)
}
mv_ivw <- get_result_from_MRObject(res_IVW_MRObject, mvdat)
mv_egger <- get_result_from_MRObject(res_egger_MRObject, mvdat)
mv_lasso <- get_result_from_MRObject(res_lasso_MRObject, mvdat)
mv_median <- get_result_from_MRObject(res_median_MRObject, mvdat)
mv_MendelianRandomization_res <- rbind(mv_ivw, mv_egger, mv_lasso, mv_median)
mv_MendelianRandomization_res %>% vroom::vroom_write("./output/tsmr/bloodsugar/mvmr_mvmrres_0716.tsv")

get_het_pleiotropy_int_from_MRObject <- function(IVW=res_IVW_MRObject, Egger=res_egger_MRObject){
  # IVW & Egger for the heterogeneity 
  # Egger intercept for horizontal pleiotropy
  res = data.frame(
    MVMR_IVW_Heterogeneity_Q = res_IVW_MRObject@Heter.Stat[1],
    MVMR_IVW_Heterogeneity_Q_P = res_IVW_MRObject@Heter.Stat[2],
    MVMR_Egger_Heterogeneity_Q = res_egger_MRObject@Heter.Stat[1],
    MVMR_Egger_Heterogeneity_Q_P = res_egger_MRObject@Heter.Stat[2],
    MVMR_Egger_pleiotropy_intercept = res_egger_MRObject@Intercept,
    MVMR_Egger_pleiotropy_intercept_P = res_egger_MRObject@Pvalue.Int
  )
  return(res)
}
h_p_res_MVMR <- get_het_pleiotropy_int_from_MRObject() %>% t %>% as.data.frame()
colnames(h_p_res_MVMR) = "MVMR Value"; h_p_res_MVMR
h_p_res_MVMR$`MVMR Value` <- format(h_p_res_MVMR$`MVMR Value`, scientific = FALSE)
h_p_res_MVMR <- h_p_res_MVMR %>% rownames_to_column(var = "Metric") %>% dplyr::select(Metric, `MVMR Value`)
h_p_res_MVMR %>% vroom::vroom_write("./output/tsmr/bloodsugar/mvmr_het_pleiotropy_0716.tsv")

# >>>>>>>>>>>>>>> MVMR——mvmr >>>>>>>>>>>>>>> #
# calculate F statistics
library(MVMR)
r_input <- format_mvmr(BXGs = cbind(mvdat[["exposure_beta"]][,1], mvdat[["exposure_beta"]][,2], mvdat[["exposure_beta"]][,3], mvdat[["exposure_beta"]][,4]), # beta-coefficient values for the exposure
                       BYG = mvdat[["outcome_beta"]],
                       seBXGs = cbind(mvdat[["exposure_se"]][,1], mvdat[["exposure_se"]][,2], mvdat[["exposure_se"]][,3], mvdat[["exposure_se"]][,4]),
                       seBYG = mvdat[["outcome_se"]],
                       RSID = rownames(mvdat$exposure_beta)); head(r_input); names(r_input); class(r_input)
strength_F <- strength_mvmr(r_input = r_input, gencov = 0)
colnames(strength_F) <- dimnames(mvdat$exposure_beta)[[2]]
strength_F %>% vroom::vroom_write("./output/tsmr/bloodsugar/mvmr_strength_0718.tsv")
pres <- pleiotropy_mvmr(r_input = r_input, gencov = 0)
res <- ivw_mvmr(r_input = r_input)