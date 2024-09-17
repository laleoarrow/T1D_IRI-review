# Coloc ---------
#' Give Messages with my color ----
#'
#' @param color prefered color
#' @param msg messages to show
#'
#' @examples
#' leo_message("This is a pink message")
#' leo_message("This is a green message","32")
#' leo_message("This is a yellow message","33")
#' leo_message("This is a blue message","34")
#' leo_message("This is a light purple message","95")
leo_message <- function(msg, color = "31") {
  message(paste0("\033[", color, "m", msg, "\033[0m\n"))
}

#' Clump data locally----
#'
#' @param dat: a datafame with `SNP`, `pval.exposure` and `id.exposure` col
#'
#' @return a subset of input dat
#' @description
#' Here is the reference:
#' - https://github.com/MRCIEU/TwoSampleMR/issues/173  
#' - https://blog.csdn.net/xiaozheng1213/article/details/126269969
clump_data_local <- function(dat) {
  require(ieugwasr); require(plinkbinr)
  leo_message(" - Clumping data locally")
  dat1 <- ieugwasr::ld_clump(
    dat = dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),
    clump_kb = 10000,
    clump_r2 = 0.001,
    bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR",
    plink_bin = plinkbinr::get_plink_exe()
  )
  dat2<- subset(dat,SNP %in% dat1$rsid)
  return(dat2)
}

#' Extract instruments locally
#'
#' @param dat: Summary statistics needs: SNP, CHR, POS, A1, A2, EAF, BETA, SE, P, Phenotype
#' @param p: P cut-off value
#'
#' @return a clumped tsmr format dataframe
extract_instruments_local <- function(dat, p = 5e-5) {
  instruments <- subset(dat, P < p)
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
    min_pval = 1e-400
  )
  instruments <- clump_data_local(instruments)
  return(instruments)
}


#' COLOC region indexed by placo lead SNP----
#' @note This is designed for coloc using lead snp from placo lead snp
#'       This is for two case/control gwas (type is set to be "cc")
#'
#' @param gwas1 gwas full data for one gwas; do not pre-prune it via maf
#' @param gwas2 gwas full data for another gwas; do not pre-prune it via maf
#' @param lead_snp self-defined df with at least `SNP`/`CHR`/`POS` column
#' @param out_fp output path for figure
#' @param out_p output path
#' @param window coloc windows size setting (unit kb)
coloc_lead_snps <- function(gwas1, gwas2, lead_snp, window = 500, 
                            out_fp = "./figure/coloc", 
                            out_p = "./output/coloc"){
  # param setting
  one_side <- window/2; leo_message(paste0(" - Coloc window <", window, "> kb"))
  length_lead_snp <- nrow(lead_snp); leo_message(paste0(" - Given <", length_lead_snp, "> lead snps"))
  # each locus
  require(pbapply);require(parallel)
  cl <- makeCluster(16, type = "FORK")
  result_list <- pbapply::pblapply(1:length_lead_snp, cl=cl, FUN = function(i) {
    leo_message(paste0(" - Selecting the No.", i, " COLOC region"))
    lead_chr = lead_snp$CHR[i]; message(paste0("lead_chr is: ", lead_chr))
    lead_pos = lead_snp$POS[i]; message(paste0("lead_pos is: ", lead_pos))
    lower_pos = lead_snp$POS[i] - one_side*1000; message(paste0("lower_pos is: ", lower_pos))
    higher_pos = lead_snp$POS[i] + one_side*1000; message(paste0("higher_pos is: ", higher_pos))
    locus <- paste0(lead_chr, ":", lower_pos, "-", higher_pos); message("locus: ",locus,"\n")
    gene <- lead_snp$nearestGene[i]; message("gene: ",gene,"\n")
      
    # subset common SNP
    tab1 <- gwas1 %>%
      dplyr::filter(CHR == lead_chr, POS > lower_pos, POS < higher_pos) %>%
      distinct(SNP, .keep_all = TRUE)
    tab2 <- gwas2 %>%
      dplyr::filter(CHR == lead_chr, POS > lower_pos, POS < higher_pos) %>%
      distinct(SNP, .keep_all = TRUE)
    commonsnps <- intersect(tab1$SNP, tab2$SNP)
    tab1 <- tab1 %>% dplyr::filter(SNP %in% commonsnps)
    tab2 <- tab2 %>% dplyr::filter(SNP %in% commonsnps)
    stopifnot(length(tab1$SNP) == length(tab2$SNP)) # just in case
    
    # foramt for harmonize
    tab1 <- tab1 %>% # give exposure name
      rename(SNP = SNP,
             beta.exposure = BETA,
             se.exposure = SE,
             effect_allele.exposure = A1,
             other_allele.exposure = A2,
             eaf.exposure = EAF,
             exposure = Phenotype) %>% 
      mutate(id.exposure = paste0(exposure, "_", i))
    tab2 <- tab2 %>% # give outcome name
      rename(SNP = SNP,
             beta.outcome = BETA,
             se.outcome = SE,
             effect_allele.outcome = A1,
             other_allele.outcome = A2,
             eaf.outcome = EAF,
             outcome = Phenotype) %>% 
      mutate(id.outcome = paste0(outcome, "_", i))
    
    # harmonize the snp
    dat_h <- TwoSampleMR::harmonise_data(tab1, tab2)
    dat_h <- dat_h %>% subset(mr_keep) %>% drop_na(SNP)
    n_snp_in_loci <- nrow(dat_h); leo_message(paste0(" - There are <", n_snp_in_loci, "> SNP in the coloc region"))
    if (n_snp_in_loci <= 1000) {leo_message("Warning: Not enough common SNP found in the coloc region. Please consider broaden the region.")}
    
    # coloc data list
    ## List of 5
    ##  $ snp     : chr [1:500] "s1" "s2" "s3" "s4" ...
    ##  $ position: int [1:500] 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ beta    : num [1:500] 0.337 0.211 0.257 0.267 0.247 ...
    ##  $ varbeta : num [1:500] 0.01634 0.00532 0.00748 0.01339 0.00664 ...
    ##  $ type    : chr "cc"
    coloc1 <- dat_h %>% {list(snp = .$SNP,
                              position = .$POS.y,
                              beta = .$beta.exposure,
                              varbeta = .$se.exposure^2,
                              pvalues = .$P.y,  # only for potential plot
                              type = "cc")} # plot_dataset(coloc1)
    coloc2 <- dat_h %>% {list(snp = .$SNP,
                              position = .$POS.x,
                              beta = .$beta.outcome,
                              varbeta = .$se.outcome^2,
                              pvalues = .$P.x, # only for potential plot
                              type = "cc")} # plot_dataset(coloc2)
    
    fmap1 <- finemap.abf(dataset = coloc1)
    fsnp1 <- fmap1[-nrow(fmap1),]; fp1 <- tail(fmap1,n = 1)
    
    fmap2 <- finemap.abf(dataset = coloc2)
    fsnp2 <- fmap2[-nrow(fmap2),]; fp2 <- tail(fmap2,n = 1)
    
    my.res <- coloc.abf(dataset1 = coloc1, dataset2 = coloc2)
    if (i == length_lead_snp) {system("echo 'The r is for coloc round' > ./figure/coloc/note.txt")}
    # fig
    figure_name <- file.path(out_fp, paste0("sen_r", i, "_", locus, ".pdf"))
    leo_message(paste0(">>> Save sensitivity figure to: ",figure_name))
    pdf(figure_name, width = 7, height = 5)
    sensitivity(my.res, rule="H4 > 0.5")
    dev.off()
    
    # results
    coloc_res_df <- tibble(
      region = paste0("region-", i),
      locus = locus,
      gene = gene,
      n_snps = my.res[["summary"]][["nsnps"]], # or n_snp_in_loci
      PP3    = my.res[["summary"]][["PP.H3.abf"]],
      PP4    = my.res[["summary"]][["PP.H4.abf"]],
      hit_snp = my.res$results %>% dplyr::filter(SNP.PP.H4 == max(SNP.PP.H4)) %>% pull(snp),
      hit_snp_beta1 = dat_h %>% dplyr::filter(SNP == hit_snp) %>% pull(beta.exposure),
      hit_snp_beta2 = dat_h %>% dplyr::filter(SNP == hit_snp) %>% pull(beta.outcome),
      hit_snp_P1 = dat_h %>% dplyr::filter(SNP == hit_snp) %>% pull(P.y),
      hit_snp_P2 = dat_h %>% dplyr::filter(SNP == hit_snp) %>% pull(P.x),
      hit_snp_lABF1 = my.res$results %>% dplyr::filter(snp == hit_snp) %>% pull(lABF.df1),
      hit_snp_lABF2 = my.res$results %>% dplyr::filter(snp == hit_snp) %>% pull(lABF.df2),
      hit_snp_sum.lABF = my.res$results %>% dplyr::filter(snp == hit_snp) %>% pull(internal.sum.lABF),
      hit_snp_pp = max(my.res$results$SNP.PP.H4),
      # fine map 1
      finemap.abf_P1 = fp1$SNP.PP,
      # causality1 = ifelse(finemap.abf_P1<0.05, "Supportive", "Negative"),
      hit1_snp = fsnp1 %>% dplyr::filter(SNP.PP == max(SNP.PP)) %>% pull(snp),
      hit1_snp_lABF = fsnp1 %>% dplyr::filter(SNP.PP == max(SNP.PP)) %>% pull(lABF.),
      hit1_snp_pp = max(fsnp1$SNP.PP),
      # fine map 1
      finemap.abf_P2 = fp2$SNP.PP,
      # causality2 = ifelse(finemap.abf_P2<0.05, "Supportive", "Negative"),
      hit2_snp = fsnp2 %>% dplyr::filter(SNP.PP == max(SNP.PP)) %>% pull(snp),
      hit2_snp_lABF = fsnp2 %>% dplyr::filter(SNP.PP == max(SNP.PP)) %>% pull(lABF.),
      hit2_snp_pp = max(fsnp2$SNP.PP),
      # coloc
      prior_p1 = my.res$priors["p1"],
      prior_p2 = my.res$priors["p2"],
      prior_p12 = my.res$priors["p12"],
      PP1    = my.res[["summary"]][["PP.H1.abf"]],
      PP2    = my.res[["summary"]][["PP.H2.abf"]]
    )
    return(coloc_res_df)
  }); stopCluster(cl)
  
  # save out
  all_coloc_results <- do.call(rbind, result_list)
  all_coloc_results %>% data.table::fwrite("./output/coloc/coloc_results.tsv", sep = "\t", col.names = T)
  rownames(all_coloc_results) <- NULL
  return(all_coloc_results)
}


# hyprcoloc ------------
.filter_min_p <- function(df) {
  df <- df %>% # mark who is duplicated
    mutate(dup = duplicated(SNP, fromLast = TRUE) | duplicated(SNP, fromLast = FALSE))
  df_dup <- df %>%  # filter duplicated SNP and find the ones with min p
    dplyr::filter(dup) %>%
    dplyr::group_by(SNP) %>%
    dplyr::filter(P == min(P)) %>%
    dplyr::slice(1) %>% # keep the 1 when p is the same
    dplyr::ungroup()
  df_non_dup <- df %>% # non-duplicated SNP
    dplyr::filter(!dup)
  rbind(df_non_dup, df_dup) %>% # brilliant!
    dplyr::select(-dup) # This way only check the duplicated ones
}
.format_hyprc_dat1 <- function(gwas, chr, bp_lower, bp_upper) {
  extracted_data <- gwas %>%
    dplyr::filter(CHR == chr, POS >= bp_lower, POS <= bp_upper) %>%
    .filter_min_p() %>%
    dplyr::select(SNP, BETA, SE) %>% 
    drop_na(SNP)
  if (nrow(extracted_data) == 0) {warning("No data extracted for the specified region.")}
  return(extracted_data)
}
.format_hyprc_dat2 <- function(gwas, chr, bp_lower, bp_upper) {
  extracted_data <- gwas %>%
    dplyr::filter(CHR == chr, BP >= bp_lower, BP <= bp_upper) %>%
    .filter_min_p() %>%
    dplyr::select(SNP, BETA, SE) %>% 
    drop_na(SNP)
  if (nrow(extracted_data) == 0) {warning("No data extracted for the specified region.")}
  return(extracted_data)
}

#' Format multiple GWAS datasets for hyprcoloc ----
#'
#' @param gwas_lists A list of GWAS data.frame. The name of the list is the trait name
#' @param chr The lead chromosome number to extract.
#' @param bp The lead bp number to extract
#' @param win hyprcoloc window setting (unit kb,defaut 500kb)
#'
#' @return A list containing beta matrix and standard error matrix.
#' @examples
#' format_hyprc_dat(list("gwas1", "gwas2"), 1, 100000, 500000)
format_hyprc_dat <- function(gwas_lists, chr, bp, win=500, 
                             gmb = F, gmb_file, gmb_id) {
  # param
  bp_lower = as.integer(bp - (win*1000)/2)
  bp_upper = as.integer(bp + (win*1000)/2)
  leo_message(paste0("\n - hyprcoloc region >>> ", chr, ":", bp_lower, "-", bp_upper))
  # extract
  extracted_data_list <- lapply(gwas_lists, function(gwas) {
    .format_hyprc_dat1(gwas, chr, bp_lower, bp_upper)
  })
  ## mbQTL
  if (gmb) {
    gmb_file_tmp <- .format_hyprc_dat2(gmb_file, chr, bp_lower, bp_upper) %>% as_tibble()
    extracted_data_list[[length(extracted_data_list)+1]] <- gmb_file_tmp
  }
  # common snps
  common_snps <- Reduce(intersect, lapply(extracted_data_list, function(x) x$SNP))
  
  # common snp and beta/se matrix
  beta_matrix <- extracted_data_list %>%
    lapply(function(x) {x %>% 
        dplyr::filter(SNP %in% common_snps) %>%
        arrange(factor(SNP, levels = common_snps)) %>%
        dplyr::select(BETA)
    }) %>%
    bind_cols() %>%
    as.matrix()
  se_matrix <- extracted_data_list %>%
    lapply(function(x) {
      x %>%
        dplyr::filter(SNP %in% common_snps) %>%
        arrange(factor(SNP, levels = common_snps)) %>%
        dplyr::select(SE)
    }) %>%
    bind_cols() %>%
    as.matrix()
  
  if (!gmb) {
    colnames(beta_matrix) <- names(gwas_list)
    colnames(se_matrix) <- names(gwas_list)
    rownames(beta_matrix) <- common_snps
    rownames(se_matrix) <- common_snps
  } else {
    colnames(beta_matrix) <- c(names(gwas_list), gmb_id)
    colnames(se_matrix) <- c(names(gwas_list), gmb_id)
    rownames(beta_matrix) <- common_snps
    rownames(se_matrix) <- common_snps
  }
  
  return(list(beta_matrix = beta_matrix, se_matrix = se_matrix))
}


#' Plot the sensitivity pheatmap for hypercoloc analysis
#'
#' @param sim.mat output from hyprcoloc::sensitivity.plot(..., similarity.matrix = TRUE)
#' @param title pheatmap title
#' @param save logi, if T, it would save the plot into save.path
#' @param save.path path to save the plot (recommended as pdf)
#'
#' @return pheatmap object
#' @export a pheatmap for sensitivity analysis
#'
#' @examples 
#' self_draw_hypr(sim.mat, title, save = T, save.path = save.path)
#' p <- self_draw_hypr(sim.mat, title, save = F)
self_draw_hypr <- function(sim.mat, title, save = F, save.path="./figure/hyprcoloc/demo.pdf"){
  library(pheatmap);library(ggplot2);library(RColorBrewer);library(grid)
  if (rownames(sim.mat)[3] %>% nchar() > 15) {
    rownames(sim.mat)[3] <- "mbQTL*"
  }
  if (colnames(sim.mat)[3] %>% nchar() > 15) {
    colnames(sim.mat)[3] <- "mbQTL*"
  }
  # number in each cell
  sim.mat.percent <- sim.mat * 100
  display_numbers <- matrix(paste0(sprintf("%.2f", round(sim.mat.percent, 2)), "%"), nrow = nrow(sim.mat))
  # self-defined color
  breaksList <- seq(0, 1, by = .02)
  colors <- colorRampPalette(c("#FFFFFF", "#4DBBD5FF"))(length(breaksList))
  # colors <- colorRampPalette(brewer.pal(9, "RdPu"))(length(breaksList))
  p <- pheatmap::pheatmap(sim.mat,
                          cluster_rows = FALSE, 
                          cluster_cols = FALSE, 
                          display_numbers = display_numbers,
                          breaks = breaksList,
                          border_color = "black", 
                          show_colnames = TRUE, 
                          show_rownames = TRUE, 
                          drop_levels = TRUE, 
                          number_color = "black",
                          color = colors, 
                          fontsize = 14,
                          fontsize_row = 14,
                          fontsize_col = 14,
                          main = title)
  if (save && !is.null(save.path)) {
    .save_pheatmap_pdf <- function(plot, save.path, width=5, height=5) {
      pdf(save.path, width=width, height=height)
      leo_message(paste0("Saving to >>> ", save.path), "93")
      grid::grid.newpage()
      grid::grid.draw(plot$gtable)
      dev.off()
    }
    .save_pheatmap_pdf(p, save.path, width = 6, height = 3)
  }
  return(plot)
}

#' perform_hypr
#' @note This function performs hyprcoloc for 2 GWAS and 1 QTL
#'
#' @param gmb_files gut mbQTL data files, it can also be other QTL data
#' @param gwas_list a reduced gwas list containing 2 gwas
#' @param prefix mbQTL prefix which is used as an identifier for sensitivity plot output
#' @param save.plot.dir path to save the sensitivity plot
#' @param save.output.dir path to save the hyprcoloc output
#' @param select columns to select for gmb data
#' @param gsub_id_pattern pattern to gsub the id for gmb/smb data
#'
#' @return
#' @export
#'
#' @examples
#' gmb_files <- list.files(path = "/Users/leoarrow/project/ref/mbqtl/skin", pattern = "\\_rsid.tsv.gz$", full.names = T)
#' gwas_list <- list(T1D = T1D, Iridocyclitis = Iridocyclitis); rm(list = c("T1D", "Iridocyclitis")); gc()
#' perform_hypr(gmb_files, gwas_list)
perform_hypr <- function (gmb_files, gwas_list, prefix = "", gsub_id_pattern = ".gz",
                          save.plot.dir = "/Users/leoarrow/project/iridocyclitis/figure/hyprcoloc",
                          save.output.dir = "/Users/leoarrow/project/iridocyclitis/output/hyprcoloc",
                          select = c("hm_rsid", "hm_chrom", "hm_pos", "hm_beta", "standard_error", "p_value")) {
 #----------------------------------main----------------------------------------#
  res_stage1s <- list() # only use for the 1st time
  prefix <- prefix; save.dir <- save.plot.dir; time1 <- Sys.time()
  for (j in length(gmb_files):1) { # length(gmb_files)
    loop_start_time <- Sys.time()
    time_spent <- loop_start_time - time1
    leo_message(paste0("\n# Loop <", j, "> Start Time: ", format(loop_start_time, "%Y-%m-%d %H:%M:%S"), " #"))
    #------------------outer loop-------------------#
    gmb_path <- gmb_files[j]
    gmb_id <- basename(gmb_path) %>% gsub(gsub_id_pattern,"",.)
    gmb_file <- fread(gmb_path, select = select, nThread = 1) %>% 
      set_names("SNP", "CHR", "BP", "BETA", "SE", "P") %>% 
      mutate_at("CHR", as.integer) %>% 
      drop_na(BETA)
    leo_message(paste0("\n# Processing <", j, ">: ", gmb_id))
    #------------------inner loop-------------------#
    res_one_gmbs <- list()
    for (i in 1:nrow(loci)) {
      # initial parameters
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
        dplyr::select(index, chr, bp, gmb_id, gene, everything())
      # sensitivity
      .check_hyprcoloc_postive <- function(reg.res) {
        if (nrow(reg.res) == 0) {return(FALSE)}
        lengths <- sapply(reg.res$traits, function(x) length(strsplit(x, ",")[[1]])) # hyprcoloced trait number
        logi <- any(lengths == 3)
        return(logi)
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
        save.path <- paste0(save.dir, "/", prefix, "_", reminder, ".pdf")
        p <- self_draw_hypr(sim.mat, title, save = T, save.path = save.path)
      } else {message(paste0("# Failed to hyprcoloc for loci: ", reminder))}
      # clean
      res_one_gmbs[[i]] <- reg.res
      leo_message("\n# clean cache #");cat(gc())
    }
    #------------------inner loop-------------------#
    res_one_gmb <- do.call("rbind", res_one_gmbs)
    rm(gmb_file)
    leo_message("\n# clean cache #");cat(gc())
    #------------------outer loop-------------------#
    loop_end_time <- Sys.time()
    last_loop_time <- as.numeric(difftime(loop_end_time, loop_start_time, units = "mins"))
    last_loop_time_formatted <- sprintf("%.2f mins", last_loop_time)
    leo_message(paste0("\n# Last Loop <", j, "> Time: ", last_loop_time_formatted, " #"))
    # return `res_one_gmb`
    res_stage1s[[j]] <- res_one_gmb
  }
  res_stage1 <- do.call("rbind", res_stage1s)
  res_stage2 <- res_stage1 %>% dplyr::filter(traits != "None") # %>% dplyr::filter(grepl("*id*", traits))
  save.output.dir <- save.output.dir
  output_path <- paste0(save.output.dir, "/", prefix, "res_stage1.tsv")
  fwrite(res_stage2, file = output_path, sep = "\t")
  leo_message(paste0("\n# All Done! #"))
  return(res_stage2)
}

#' perform_hypr with multi-core CPU using `pbmcapply`
#' @note This function performs hyprcoloc for 2 GWAS and 1 QTL
#'
#' @param gmb_files gut mbQTL data files, it can also be other QTL data
#' @param gwas_list a reduced gwas list containing 2 gwas
#' @param prefix mbQTL prefix which is used as an identifier for sensitivity plot output
#' @param save.plot.dir path to save the sensitivity plot
#' @param save.output.dir path to save the hyprcoloc output
#' @param select columns to select for gmb data
#' @param gsub_id_pattern pattern to gsub the id for gmb/smb data
#'
#' @return
#' @export
#'
#' @examples
#' gmb_files <- list.files(path = "/Volumes/T9/mbqtl/412DutchMicrobiomeProject/source_data", pattern = "\\h.tsv.gz$", full.names = T)
#' gwas_list <- list(T1D = T1D, Iridocyclitis = Iridocyclitis); rm(list = c("T1D", "Iridocyclitis")); gc()
#' perform_hypr(gmb_files, gwas_list)
mc_perform_hypr <- function (gmb_files, gwas_list, prefix = "", gsub_id_pattern = ".gz",
                             save.plot.dir = "/Users/leoarrow/project/iridocyclitis/figure/hyprcoloc",
                             save.output.dir = "/Users/leoarrow/project/iridocyclitis/output/hyprcoloc",
                             select = c("hm_rsid", "hm_chrom", "hm_pos", "hm_beta", "standard_error", "p_value"),
                             mc = 14) {
  #----------------------------------main----------------------------------------#
  prefix <- prefix; save.dir <- save.plot.dir
  res_stage1s <- pbmcapply::pbmclapply(1:length(gmb_files), function(j) {
    #------------------outer loop-------------------#
    gmb_path <- gmb_files[j]
    # gmb_path <- smb_files[j] # bugging
    gmb_id <- basename(gmb_path) %>% gsub(gsub_id_pattern,"",.)
    gmb_file <- fread(gmb_path, select = select, nThread = 1) %>% 
      set_names("SNP", "CHR", "BP", "BETA", "SE", "P") %>% 
      mutate_at("CHR", as.integer) %>% 
      drop_na(BETA)
    leo_message(paste0("\n# Processing <", j, ">: ", gmb_id))
    #------------------inner loop-------------------#
    res_one_gmbs <- list()
    for (i in 1:nrow(loci)) {
      # initial parameters
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
        dplyr::select(index, chr, bp, gmb_id, gene, everything())
      # sensitivity
      .check_hyprcoloc_postive <- function(reg.res) {
        if (nrow(reg.res) == 0) {return(FALSE)}
        lengths <- sapply(reg.res$traits, function(x) length(strsplit(x, ",")[[1]])) # hyprcoloced trait number
        logi <- any(lengths == 3)
        return(logi)
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
        save.path <- paste0(save.dir, "/", prefix, "_", reminder, ".pdf")
        p <- self_draw_hypr(sim.mat, title, save = T, save.path = save.path)
      } else {message(paste0("# Failed to hyprcoloc for loci: ", reminder))}
      # clean
      res_one_gmbs[[i]] <- reg.res
      leo_message("\n# clean cache #");cat(gc())
    }
    #------------------inner loop-------------------#
    res_one_gmb <- do.call("rbind", res_one_gmbs)
    rm(gmb_file)
    leo_message("\n# clean cache #");cat(gc())
    #------------------outer loop-------------------#
    res_one_gmb$j <- j
    return(res_one_gmb)
  }, mc.cores = mc)
  res_stage1 <- do.call("rbind", res_stage1s) %>% as.data.frame()
  res_stage2 <- res_stage1 %>% dplyr::filter(traits != "None") # %>% dplyr::filter(grepl("*GCST*", traits))
  save.output.dir <- save.output.dir
  output_path <- paste0(save.output.dir, "/", prefix, "res_stage1.tsv")
  fwrite(res_stage2, file = output_path, sep = "\t")
  leo_message(paste0("\n# All Done! #"))
  return(res_stage2)
}