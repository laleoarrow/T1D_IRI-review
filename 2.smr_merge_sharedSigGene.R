#' This is used for locating shared genes identified in SMR from eQTL to GWAS.
#' 
#' @param fdr_res_path1 Path to the first SMR results file.
#' @param fdr_res_path2 Path to the second SMR results file.
#' @param gwas1_name Name of the first GWAS.
#' @param gwas2_name Name of the second GWAS.
#'
#' @return A combined CSV file with shared genes.
#' 
#' @note For debugging, you may use the following example:
#' fdr_res_path1 = "/Users/leoarrow/project/iridocyclitis/output/GTEx49_eqtl/fdr_sig/Artery_Tibial-type1d.eqtl.sig.smr"
#' fdr_res_path2 = "/Users/leoarrow/project/iridocyclitis/output/GTEx49_eqtl/fdr_sig/Artery_Tibial-iri_meta.eqtl.sig.smr"

library(tidyverse)
library(vroom)
compare_and_merge <- function(fdr_res_path1, fdr_res_path2, gwas1_name, gwas2_name) {
  smr1 <- vroom(fdr_res_path1, show_col_types = FALSE, progress = FALSE)
  smr2 <- vroom(fdr_res_path2, show_col_types = FALSE, progress = FALSE)
  # common_genes <- intersect(smr1$Gene, smr2$Gene) ; n <- length(common_genes) # for the if(dim(smr_all)[1] != n) below
  common_probe <- intersect(smr1$probeID, smr2$probeID); n <- length(common_probe) 
  smr1 <- smr1 %>% filter(probeID %in% common_probe) %>% arrange(probeID) 
  smr2 <- smr2 %>% filter(probeID %in% common_probe) %>% arrange(probeID) 
  # merge them with suffixes
  tissue <- strsplit(basename(fdr_res_path1),"@")[[1]][1]
  smr_all <- merge(smr1, smr2, by = "probeID", all = TRUE)
  if (dim(smr_all)[1] != n) {stop(paste0("Error: The number of shared probeID is not equal to ", n))} # just in case; Decrepited now
  # build the header columns  
  # n <- dim(smr_all)[1]
  res <- data.frame(
    tissue = rep(tissue, n),
    gwas1 = rep(gwas1_name, n),
    gwas2 = rep(gwas2_name, n),
    n_shared_genes = rep(n, n),
    probeID = common_probe
  ) %>% arrange(probeID)
  # merge the header columns with results; and make sure the order is satisfying please
  m <- merge(res, smr_all, by = "probeID", all = TRUE) # m for merged_data
  return(m)
}

####################################### main #######################################
# The first part here is made for GTEx; if it is for other source, please refer to the part 2.
# 
# Part 1: for GTEx  ############################
# Specify parameters and execute the function
gwas_names <- c("t1d1", "iri3")
fdr_path <- "/Users/leoarrow/project/iridocyclitis/output/smr/eqtl/GTEx49/fdr_sig" # To be more precise, it's actually the `fdr_sig` dir
fdr_ress <- list.files(path = fdr_path, pattern = "\\.sig.smr$", full.names = TRUE)
# Go looooooop! GTEx ###########################
# final_df <- data.frame()
# looped_tissues <- c() # to avoid duplicated tissues loop
# # Important Note: 需要先手动改名让tissue_gwas变为tissue@gwas的格式!!!!!!!!!!!!!!!!!!
# # For checking in one iteration: fdr_res = fdr_ress[1]
# for (fdr_res in fdr_ress) {
#   fdr_basename <- basename(fdr_res)
#   dir_path <- dirname(fdr_res)
#   t <- strsplit(basename(fdr_basename),"@")[[1]][1] # t for tissue
#   if (t %in% looped_tissues) {next} # skip if already looped
#   message(paste0("* Processing tissue: ", t))
#   
#   gwas1_name <- gwas_names[1]
#   gwas2_name <- gwas_names[2]
#   
#   fdr_res_path1 <- paste0(dir_path, "/", t, "@", gwas1_name, ".sig.smr") # This part sometimes need to be adjusted!
#   fdr_res_path2 <- paste0(dir_path, "/", t, "@", gwas2_name, ".sig.smr") # This part sometimes need to be adjusted!
#   
#   # Chunk for mannual check the data (optional)
#   # if (file.exists(fdr_res_path1) & file.exists(fdr_res_path2)) {stop("Error: One of the files does not exist. Please check the file path.")}
#   if (xor(file.exists(fdr_res_path1), file.exists(fdr_res_path2))) {next} # xor: only one of them exists # happen when one smr has now sig res.
#   
#   tmp_df <- compare_and_merge(fdr_res_path1, fdr_res_path2, gwas1_name, gwas2_name)
#   looped_tissues <- c(looped_tissues, t)
#   final_df <- rbind(final_df, tmp_df)
# }
# # write the results to a CSV file
# output_path <- paste0(fdr_path, "/smr_combined.csv"); message(paste0(" - Writing to ", output_path))
# write.csv(final_df, output_path, row.names = FALSE)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Part 2: for other sources ############################
gwas_names <- c("t1d1", "iri3")
fdr_path <- "/Users/leoarrow/project/iridocyclitis/output/smr/pqtl/ppp/fdr_sig" # To be more precise, it's actually the `fdr_sig` dir
fdr_ress <- list.files(path = fdr_path, pattern = "\\.sig.smr$", full.names = TRUE)
# Go looooooop! eqtlgen ###########################
final_df <- data.frame()
looped_tissues <- c() # to avoid duplicated tissues loop
for (fdr_res in fdr_ress) {
  fdr_basename <- basename(fdr_res)
  dir_path <- dirname(fdr_res)
  t <- strsplit(basename(fdr_basename),"@")[[1]][1] # t for tissue
  if (t %in% looped_tissues) {next} # skip if already looped
  message(paste0("* Processing: ", t))
  
  gwas1_name <- gwas_names[1]
  gwas2_name <- gwas_names[2]
  
  fdr_res_path1 <- paste0(dir_path, "/", t, "@", gwas1_name, ".sig.smr")
  fdr_res_path2 <- paste0(dir_path, "/", t, "@", gwas2_name, ".sig.smr")
  tmp_df <- compare_and_merge(fdr_res_path1, fdr_res_path2, gwas1_name, gwas2_name)
  looped_tissues <- c(looped_tissues, t)
  final_df <- rbind(final_df, tmp_df)
}

# write the results to a CSV file
t <- sub(".*/([^/]+)/[^/]+/?$", "\\1", fdr_path)
output_path <- paste0(fdr_path, "/", t, "@combined.csv"); message(paste0(" - Writing to ", output_path))
write.csv(final_df %>% arrange(ProbeChr.x, Probe_bp.x), output_path, row.names = FALSE)

# Two neccessary function is `2.smr_mygwas.ma.R`
path <- "/Users/leoarrow/project/iridocyclitis/output/smr/pqtl/ppp"; fdr_path <- file.path(path,"fdr_sig"); message(fdr_path)
leo_smr_adjust_loop(path) # 保存
leo_smr_sig_loop(folder_path = path, save_path = fdr_path)

##########Additional for eqtlgen for transformation from esemble id to gene id########
#' @title BioMart
if (!requireNamespace("biomaRt", quietly = TRUE)) {install.packages("biomaRt")}
library(biomaRt)
x1=vroom::vroom("/Users/leoarrow/project/iridocyclitis/output/smr/eqtl/eqtlgen/fdr_sig/smr_eqtlgen@iri3.sig.smr") 
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
gene_list <- x1$Gene # Gene list
genes <- getBM(attributes = c('ensembl_gene_id', 'gene_biotype', 'external_gene_name'),
               filters = 'ensembl_gene_id',
               values = gene_list,
               mart = ensembl)
final_results <- do.call(rbind, results_list)
#' @title clusterProfiler
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
library(vroom)
x1=vroom("/Users/leoarrow/project/iridocyclitis/output/smr/eqtl/eqtlgen/fdr_sig/smr_eqtlgen@t1d1.sig.smr") 
x2=vroom("/Users/leoarrow/project/iridocyclitis/output/smr/eqtl/eqtlgen/fdr_sig/smr_eqtlgen@iri3.sig.smr") 
# gene_list <- x1$Gene # Gene list
# gene_list <- gsub("\\..*", "",  gene_list) # remove version number "XXX.X"
ENSEMBL_to_SYMBOL <- function(gene_list){
  gene.symbol <- bitr(geneID = gene_list, 
                      fromType = "ENSEMBL",
                      toType = c("SYMBOL"),
                      OrgDb = org.Hs.eg.db)
}
gene1 <- ENSEMBL_to_SYMBOL(x1 %>% pull(Gene))
gene2 <- ENSEMBL_to_SYMBOL(x2 %>% pull(Gene))
intersect(gene1, gene2)
