#' Give Messages with my color ----
#'
#' @param color prefered color
#' @param msg messages to show; it should be a word pasted pharse.
#'
#' @examples
#' leo_message("This is a pink message")
#' leo_message("This is a green message","32")
#' leo_message("This is a blue message","34")
#' leo_message("This is a light purple message","95")
leo_message <- function(msg, color = "31") {
 message(paste0("\033[", color, "m", msg, "\033[0m\n"))
}

#' Clump data locally----
#'
#' @param dat a datafame with `SNP`, `pval.exposure` and `id.exposure` col
#' @param pop Super-population to use as reference panel. Options are: AFR, AMR, EAS, EUR, SAS
#' @param bfile If this is provided then will go ld locally
#'
#' @return a subset of input dat with independant SNP
#' @description
#' Here is the reference:
#' https://github.com/MRCIEU/TwoSampleMR/issues/173  
#' https://blog.csdn.net/xiaozheng1213/article/details/126269969
#' library(ieugwasr) 
#' library(plinkbinr) # devtools::install_github("explodecomputer/plinkbinr")
#' plinkbinr::get_plink_exe()
clump_data_local <- function(dat, pop = NULL, bfile = NULL, plink_bin = NULL) {
 require(ieugwasr); require(plinkbinr)
 leo_message(" - Clumping data locally")
 leo_message(paste0(" - Reminder: bfile located: ", "/Users/leoarrow/project/ref/1kg.v3/EUR"))
 if (is.null(pop) & is.null(bfile)) { stop("Must indicate LD method") }
 if (!is.null(pop) & is.null(bfile)) { leo_message("Performing LD online") }
 if (is.null(pop) & !is.null(bfile)) { leo_message("Performing LD locally") }
 # if (!is.null(pop) & !is.null(bfile)) { stop("Only one LD method can be used") }
 dat1 <- ieugwasr::ld_clump(
  dat = dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),
  clump_kb = 10000,
  clump_r2 = 0.001, 
  pop = pop,
  bfile = bfile,
  plink_bin = plinkbinr::get_plink_exe()
 )
 dat2<- subset(dat, SNP %in% dat1$rsid)
 return(dat2)
}

#' Extract instruments locally
#'
#' @param dat: Summary statistics needs: SNP, CHR, POS, A1, A2, EAF, BETA, SE, P, Phenotype, N
#' @param p: P cut-off value
#' @param pop: param for clump_data_local
#' @param bfile: param for clump_data_local
#' 
#' @return a clumped tsmr format dataframe
#' @details
#' P column name is `P`
#' 
extract_instruments_local <- function(dat, p = 5e-08, N = "Neff", pop = "EUR", bfile = NULL, plink_bin = NULL) {
 instruments <- subset(dat, P < p) # 调试 dat = gwas1[1:1000,]
 leo_message(paste0(" - A total of <",dim(instruments)[1], "> SNP passed the P-value threshold"))
 leo_message(" - Using effective N as default!!!")
 instruments <- format_data(
  instruments, 
  type = "exposure",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  chr_col = "CHR",
  pos_col = "POS",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "EAF",
  beta_col = "BETA",
  se_col = "SE",
  pval_col = "P",
  samplesize_col = N,
  min_pval = 1e-200
 )
 instruments <- clump_data_local(instruments, pop = pop, bfile = bfile) # F2中执行了F1（本地Clump）
 return(instruments)
}

#' format outcome data
#'
#' @param dat a dataframe for outcome with SNP, CHR, POS, A1, A2, EAF, BETA, SE, P, Phenotype, samplesize columns
#' @param snp a str vector out of iv$SNP
#' @param N N column name for sample size (effective or observed)
#'
#' @return a tsmr format outcome dataframe
#' 
format_outcome <- function(dat, snp = iv$SNP, N = "Neff") {
 leo_message("You should check samplesize column (effective or observed)")
 leo_message(paste0("Now using <", N, "> for samplesize indice !!!!!!!!"))
 out <- format_data(
  dat, 
  snps = snp,
  type = "outcome",
  # defaut col names
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  chr_col = "CHR",
  pos_col = "POS",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "EAF",
  beta_col = "BETA",
  se_col = "SE",
  pval_col = "P",
  samplesize_col = N,
  min_pval = 1e-200
 )
 return(out)
}

#' `perform_mr_for_one_pair` perform a MR for one pair of exp and out
#' 
#' @param dat_h harmonized data
#' @param exp a str indicating the exposure in dat_h
#' @param out a str indicating the outcome in dat_h
#' @param res_dir dir path where the result fo MR analysis stored
#' @param fig_dir 
#' @param save_plot only save the plot if T, defaut T.
#' 
#' @return list(res_pair=res_pair, res_pair_presso=res_pair_presso)
#' @example
#' clusterEvalQ(cl, {
#' library(vroom)
#' library(tidyverse)
#' library(TwoSampleMR)
#' library(ggplot2)
#' library(ggsci)
#' source("./code/1.0.tsmr_packages.R")
#' })
#' clusterExport(cl, varlist = c("uni_pair","dat_h"))
mr_one_pair <- function(dat_h, exp = "", out = "", save_plot = T, res_dir= "./output/tsmr", fig_dir="./figure/tsmr") {
  # Initialize
  pairname <- paste0(exp, " VS ", out); leo_message(paste0(" - MR Pair: ", pairname))
  if (!dir.exists(res_dir)) {dir.create(res_dir, recursive = T)}; leo_message(paste0(" - Setting Check: `res_dir` using ", res_dir))
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive = T)}; leo_message(paste0(" - Setting Check: `fig_dir` using ", fig_dir))
  
  # MR for one specific pair of exposure and outcome
  dat_h_pair <- dat_h %>% filter(exposure == exp, outcome == out)
  dat_h_pair <- subset(dat_h_pair, mr_keep); rownames(dat_h_pair) <- NULL
  if (nrow(dat_h_pair) < 2) {
    res_pair <- mr(dat_h_pair, method_list =  c("mr_wald_ratio"))
  } else {
    res_pair <- mr(dat_h_pair, method_list =  c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  }
  res_pair <- generate_odds_ratios(res_pair)
  res_pair <- res_pair %>% 
    mutate(R2_total = scales::percent(sum(dat_h_pair$R2), accuracy = 0.01),
           Neff_max_exp = as.integer(max(dat_h_pair$samplesize.exposure[1])),
           Neff_max_out = as.integer(max(dat_h_pair$samplesize.outcome[1])))
  
  # mr_heterogeneity & pleiotropy test
  heterogeneity_pair <- mr_heterogeneity(dat_h_pair)
  pleiotropy_res <- mr_pleiotropy_test(dat_h_pair)
  presso_res <- tryCatch(expr = {run_mr_presso(dat_h_pair, SignifThreshold = 0.05, NbDistribution = ifelse(nrow(dat_h_pair)/0.05<1000, 1000, nrow(dat_h_pair)/0.05))},
                         error = function(e) {message("Error in MR-PRESSO: ", e$message);return(NULL)})
  
  # integrate results into `res_pair`
  res_pair$heterogeneity_IVW_Q <- ifelse(!is.null(heterogeneity_pair$Q[which(heterogeneity_pair$method == "Inverse variance weighted")]), heterogeneity_pair$Q[which(heterogeneity_pair$method == "Inverse variance weighted")], "NA")
  res_pair$heterogeneity_IVW_Q_df <- ifelse(!is.null(heterogeneity_pair$Q_df[which(heterogeneity_pair$method == "Inverse variance weighted")]), heterogeneity_pair$Q_df[which(heterogeneity_pair$method == "Inverse variance weighted")], "NA")
  res_pair$heterogeneity_IVW_Q_pval <- ifelse(!is.null(heterogeneity_pair$Q_pval[which(heterogeneity_pair$method == "Inverse variance weighted")]), heterogeneity_pair$Q_pval[which(heterogeneity_pair$method == "Inverse variance weighted")], "NA")
  res_pair$heterogeneity_Egger_Q <- ifelse(!is.null(heterogeneity_pair$Q[which(heterogeneity_pair$method == "MR Egger")]), heterogeneity_pair$Q[which(heterogeneity_pair$method == "MR Egger")], "NA")
  res_pair$heterogeneity_Egger_Q_df <- ifelse(!is.null(heterogeneity_pair$Q_df[which(heterogeneity_pair$method == "MR Egger")]), heterogeneity_pair$Q_df[which(heterogeneity_pair$method == "MR Egger")], "NA")
  res_pair$heterogeneity_Egger_Q_pval <- ifelse(!is.null(heterogeneity_pair$Q_pval[which(heterogeneity_pair$method == "MR Egger")]), heterogeneity_pair$Q_pval[which(heterogeneity_pair$method == "MR Egger")], "NA")
  res_pair$pleiotropy_Egger_intercept <- ifelse(!is.null(pleiotropy_res$egger_intercept), pleiotropy_res$egger_intercept, "NA")
  res_pair$pleiotropy_Egger_intercept_se <- ifelse(!is.null(pleiotropy_res$se), pleiotropy_res$se, "NA")
  res_pair$pleiotropy_Egger_intercept_pval <- ifelse(!is.null(pleiotropy_res$pval), pleiotropy_res$pval, "NA")
  
  # presso results
  if (!is.null(presso_res)) {
    res_pair_presso <- presso_res[[1]]$`Main MR results` %>% mutate(Exposure = exp, Outcome = out) %>% select(Exposure, Outcome, everything()) 
    outlier_detected <- ifelse(is.na(res_pair_presso$`Causal Estimate`[2]), "NO", "YES")
    if (outlier_detected == "YES") {
      res_pair_presso$outlier_detected <- outlier_detected
      res_pair_presso$Global_RSSobs <- presso_res[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs
      res_pair_presso$Global_Pvalue <- presso_res[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
      res_pair_presso$Distortion_Coefficient <- presso_res[[1]]$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`
      res_pair_presso$Distortion_Pvalue <- presso_res[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
    } else {
      res_pair_presso$outlier_detected <- c(outlier_detected, "NA")
      res_pair_presso$Global_RSSobs <- c(presso_res[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs, NA)
      res_pair_presso$Global_Pvalue <- c(presso_res[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue, NA)
      res_pair_presso$Distortion_Coefficient <- NA
      res_pair_presso$Distortion_Pvalue <- NA
    }
  } else {
    res_pair_presso <- NULL
  }
  
  # res_path <- file.path(res_dir, paste0("tsmr_", exp, "_VS_", out, ".tsv")); leo_message(paste0(" >>> Save result to: ", res_path))
  # vroom::vroom_write(res_pair, res_path)
  # res_path <- file.path(res_dir, paste0("tsmr_presso_", exp, "_VS_", out, ".tsv")); leo_message(paste0(" >>> Save result to: ", res_path))
  # vroom::vroom_write(res_pair_presso, res_path)
  
  if (save_plot) {
    #' >>>>>>>>>>>>>  draw and save draw   >>>>>>>>>>>>>  
    #' p1: Scatter plot
    #' p2: Forest plot
    #' p3: Funnel plot
    #' p4: LOO plot
    #' For Forest: width = 7, height = 8 (T1D)
    #' For Forest: width = 7, height = 10 (T2D)
    #' For others: width = 7, height = 6
    #' fig_dir <- "./figure/tsmr"
    # >>>>>>>>>>>>> Scatter Plot >>>>>>>>>>>>>>>>>>>
    # p1 <- mr_scatter_plot_modified(mr_results = res_pair, dat = dat_h_pair); p1[[1]]
    p1 <- mr_scatter_plot(mr_results = res_pair, dat = dat_h_pair)
    p1[[1]] <- p1[[1]]+ggsci::scale_color_npg()+ggsci::scale_fill_npg()+ggplot2::theme_classic()+
      ggplot2::theme(legend.position = "top",legend.direction = "vertical",
                     axis.title.x = element_text(size = 16, colour = "black"),  # X label
                     axis.title.y = element_text(size = 16, colour = "black"),  # Y label
                     axis.text = element_text(size = 14, colour = "black"), # x/y axis size
                     legend.title = element_text(size = 14, face = "italic"), # legend title
                     legend.text = element_text(size = 12),
                     legend.background = element_blank()
      )+ggplot2::guides(colour = ggplot2::guide_legend(ncol = 3))
    p1[[1]][["layers"]][[1]][["aes_params"]]$colour <- "black" # error bar (|)
    p1[[1]][["layers"]][[1]][["aes_params"]]$alpha <- 0.75 # error bar (|)
    p1[[1]][["layers"]][[2]][["aes_params"]]$colour <- "black" # error bar (-)
    p1[[1]][["layers"]][[2]][["aes_params"]]$alpha <- 0.75 # error bar (-)
    p1[[1]][["layers"]][[3]][["aes_params"]]$colour <- "black" # snp point
    p1[[1]][["layers"]][[4]][["aes_params"]]$size <- 1 # geom_abline line
    p1[[1]][["layers"]][[4]][["aes_params"]]$alpha <- 0.75 # geom_abline line
    fig_path <- file.path(fig_dir, paste0("tsmr_scatter_", exp, "_VS_", out, ".pdf")); leo_message(paste0(" >>> Save figure to: ", fig_path))
    ggsave(fig_path, plot = p1[[1]], width = 7, height = 6)
    # >>>>>>>>>>>>> mr_forest_plot / mr_funnel_plot >>>>>>>>>>>>>>>>>>>
    res_single_pair <- mr_singlesnp(dat_h_pair)
    p2 <- mr_forest_plot(res_single_pair)
    p3 <- mr_funnel_plot(res_single_pair)
    p2[[1]] <- p2[[1]]+ggsci::scale_color_npg()+ggsci::scale_fill_npg()+ggplot2::theme_classic()+
      ggplot2::theme(legend.position = "none", 
                     axis.title.y = element_text(colour = "black"), 
                     axis.title.x = element_text(colour = "black"),
                     axis.text = element_text(colour = "black"))
    p3[[1]] <- p3[[1]]+ggsci::scale_color_npg()+ggsci::scale_fill_npg()+ggplot2::theme_classic()+
      ggplot2::theme(legend.title = element_text(size = 14, face = "italic"), 
                     legend.text = element_text(size = 12),
                     legend.position = "top",legend.direction = "horizontal",
                     axis.title.y = element_text(size = 16,colour = "black"), 
                     axis.title.x = element_text(size = 16,colour = "black"),
                     axis.text = element_text(size = 14,colour = "black"))
    fig_path <- file.path(fig_dir, paste0("tsmr_forrest_", exp, "_VS_", out, ".pdf")); leo_message(paste0(" >>> Save figure to: ", fig_path))
    if (grepl("T2D*", exp)) {
      ggsave(fig_path, plot = p2[[1]], width = 7, height = 14)
    } else {
      ggsave(fig_path, plot = p2[[1]], width = 7, height = 8)}
    fig_path <- file.path(fig_dir, paste0("tsmr_funnel_", exp, "_VS_", out, ".pdf")); leo_message(paste0(" >>> Save figure to: ", fig_path))
    ggsave(fig_path, plot = p3[[1]], width = 7, height = 6)
    # >>>>>>>>>>>>> LOO Plot >>>>>>>>>>>>>>>>>>>
    loo_pair = mr_leaveoneout(dat_h_pair)
    p4 <- mr_leaveoneout_plot(loo_pair)
    p4[[1]] <- p4[[1]]+ggsci::scale_color_npg()+ggsci::scale_fill_npg()+ggplot2::theme_classic()+
      ggplot2::theme(legend.position = "none", 
                     axis.title.y = element_text(colour = "black"), 
                     axis.title.x = element_text(colour = "black"),
                     axis.text = element_text(colour = "black"))
    fig_path <- file.path(fig_dir, paste0("tsmr_loo_", exp, "_VS_", out, ".pdf")); leo_message(paste0(" >>> Save figure to: ", fig_path))
    if (grepl("T2D*", exp)) {
      ggsave(fig_path, plot = p4[[1]], width = 7, height = 14)
    } else {
      ggsave(fig_path, plot = p4[[1]], width = 7, height = 8)}
  }
  
  # reture `res-pair` and `presso` results as a list
  return(list(
    res_pair=res_pair %>% select(-id.exposure,-id.outcome) %>% select(exposure, outcome, everything()), 
    res_pair_presso=res_pair_presso
    ))
}


#' The `mr_scatter_plot` in the `TwoSampleMR` package could be better.
#' This is a modified version of `mr_scatter_plot` in the `TwoSampleMR` package
#' When the slope of two method is really close, the scatter plot may not properly plot them!
#' Change the line type here here then.
#' 
#' @param mr_results same as TwoSampleMR::mr_scatter_plot
#' @param dat same as TwoSampleMR::mr_scatter_plot
#' @example 
#' p1 <- mr_scatter_plot_modified(mr_results = res_pair, dat = dat_h_pair)
#' print(p1[[1]])
mr_scatter_plot_modified <- function(mr_results, dat) {
  #' library(ggplot2);library(ggsci);library(TwoSampleMR);library(dplyr)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d) {
    d <- plyr::mutate(d)
    if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & 
                      id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if ("MR Egger" %in% mrres$method) {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, 
                                  d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    if ("MR Egger (bootstrap)" %in% mrres$method) {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, 
                                            d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    # debugging
    # print(mrres)
    
    p <- ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, y = beta.outcome)) + 
      ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome), 
                             colour = "black", alpha= 0.75, width = 0) + 
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), 
                              colour = "black", alpha= 0.75, height = 0) + 
      ggplot2::geom_point(alpha= 0.75) + 
      ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, slope = b, colour = method, linetype = method), 
                           show.legend = TRUE, size = 1.5, alpha = 0.8) + 
      # ggplot2::scale_colour_manual(values = c("Inverse variance weighted" = "red", 
      #                                         "MR Egger" = "blue",
      #                                         "Weighted median" = "green",
      #                                         "IVW radial" = "purple")) + 
      ggplot2::scale_linetype_manual(values = c("Inverse variance weighted" = "solid", # change setting here!!!!!
                                                "MR Egger" = "solid",
                                                "Weighted median" = "solid",
                                                "IVW radial" = "dotted")) + 
      ggsci::scale_color_npg()+ggsci::scale_fill_npg()+ggplot2::theme_classic()+
      ggplot2::labs(colour = "MR Test", linetype = "MR Test", x = paste("SNP effect on", d$exposure[1]), 
                    y = paste("SNP effect on", d$outcome[1])) + 
      ggplot2::theme(legend.position = "top", legend.direction = "vertical") + 
      ggplot2::theme(legend.position = "top",legend.direction = "vertical",
                     axis.title.x = element_text(size = 16, colour = "black"),  # X label
                     axis.title.y = element_text(size = 16, colour = "black"),  # Y label
                     axis.text = element_text(size = 14, colour = "black"), # x/y axis size
                     legend.title = element_text(size = 14, face = "italic"), # legend title
                     legend.text = element_text(size = 12),
                     legend.background = element_blank()
      )+ggplot2::guides(colour = ggplot2::guide_legend(ncol = 4), linetype = ggplot2::guide_legend(ncol = 4))
    return(p)
  })
  mrres
}

#' `mrlap_one_pair` perform the MR-lap for one pair of exposure and outcome
#' https://github.com/n-mounier/MRlap?tab=readme-ov-file
#'
#' @param exposure_name exposure_name
#' @param outcome_name outcome_name
#' @param exposure_data exposure_data
#' @param outcome_data outcome_data
#' @param ld_path ldsc required file
#' @param hm3_path ldsc required file; tutorial use nomhc version list
#' @param log logi, if T, would save ldsc log to log_path, defaut F
#' @param log_path path to where you wanna store the ldsc log
#' 
mrlap_one_pair <- function(exposure_data, outcome_data, exposure_name, outcome_name, ld_path, hm3_path, log = F, log_path = ".") {
  leo_message(paste0(" - MR-lap for: ", exposure_name, " VS ", outcome_name))
  if (log) {leo_message(paste0("log file would be stored at >>>"), log_path)}
  # locate iv_mrlap snps
  iv_mrlap <- extract_instruments_local(exposure_data, p = 5e-08, N = "Neff", bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR", plink_bin = plinkbinr::get_plink_exe())
  iv_mrlap <- iv_mrlap %>% 
    mutate(R2 = (2*(beta.exposure)^2*(eaf.exposure*(1-eaf.exposure)))/((se.exposure^2) * samplesize.exposure),
           F = R2*(samplesize.exposure-2)/(1-R2)) %>% 
    filter(F > 10) %>% 
    mutate(maf.exposure = ifelse(eaf.exposure<0.5, eaf.exposure, 1-eaf.exposure)) %>%
    filter(maf.exposure > 0.01)
  iv_mrlap$id.exposure <- iv_mrlap$exposure
  out_mrlap <- format_outcome(outcome_data, iv_mrlap$SNP)
  out_mrlap <- out_mrlap %>%
    mutate(maf.outcome = ifelse(eaf.outcome<0.5, eaf.outcome, 1-eaf.outcome)) %>% 
    subset(maf.outcome>0.01) %>% 
    mutate(id.outcome = outcome)
  dat_h_mrlap <- harmonise_data(iv_mrlap, out_mrlap) %>% subset(mr_keep)
  
  # merge exposure and outcome with w_hm3.xxx.snplist (with/without MHC)
  hm3_snp.list <- fread(hm3_path) %>% select(SNP)
  iv_snp.list <- dat_h_mrlap %>% select(SNP)
  snp.list <- rbind(hm3_snp.list, iv_snp.list) %>% unique()
  exposure_data <- exposure_data %>% filter(SNP %in% snp.list$SNP) %>% select(-Neff)
  outcome_data <- outcome_data %>% filter(SNP %in% snp.list$SNP) %>% select(-Neff)
  
  # MR-lap
  wd <- getwd()
  if (log) {setwd(log_path)}
  A <- MRlap(exposure = exposure_data,
             exposure_name = exposure_name,
             outcome = outcome_data,
             outcome_name = outcome_name,
             ld = ld_path,
             hm3 = hm3_path,
             save_logfiles = log,
             do_pruning = FALSE,
             user_SNPsToKeep = dat_h_mrlap$SNP
             )
  setwd(wd)
  
  # Results
  res_lap_pair <- data.frame(
    exposure = exposure_name,
    outcome = outcome_name,
    m_IVs = A[[1]]$m_IVs,
    observed_effect = A[[1]]$observed_effect,
    observed_se = A[[1]]$observed_effect_se,
    observed_p = A[[1]]$observed_effect_p,
    corrected_effect = A[[1]]$corrected_effect,
    corrected_effect_se = A[[1]]$corrected_effect_se,
    corrected_effect_p = A[[1]]$corrected_effect_p,
    test_difference = A[[1]]$test_difference,
    p_difference = A[[1]]$p_difference
  )
  return(res_lap_pair)
}