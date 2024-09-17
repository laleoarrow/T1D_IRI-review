#' This is a toolkit for converting the `#CHR`/`POS` to `rsid` in the GWAS summary statistics.
# CHR:POS to RSID ----
library(data.table)
gwas_path <- "/Users/leoarrow/project/iridocyclitis/data/bloodsugar/GCST90002238_buildGRCh37.tsv.gz"
dat <- fread(gwas_path)
# dat <- setDT(gwas2) # mannually define
colnames(dat)

# 转rs id
library("SNPlocs.Hsapiens.dbSNP155.GRCh37")
# library("SNPlocs.Hsapiens.dbSNP155.GRCh38")
library(BSgenome)
# setDT(dat)
dat[, ranges := paste(`chromosome`, ":", base_pair_location, "-", base_pair_location, sep = "")]
# try
# dat_head <- head(dat)
# dat_head[, ranges := sprintf("%s:%s-%s", `#CHR`, POS, POS)]

snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37
# snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
snp.res <- snpsByOverlaps(snps, GRanges(dat$ranges))


library(stringr)
library(tidyverse)
snp.res.dt <- as.data.table(snp.res)
snp.res.dt[, ranges := stringr::str_c(seqnames, ":", pos, "-", pos)]
trans.dat <- merge(snp.res.dt, dat, by = "ranges") # colnames(trans.dat)
columns_to_remove <- c("strand", "alleles_as_ambig", "ranges", "seqnames", "pos"); trans.dat[, (columns_to_remove) := NULL]
trans.dat <- trans.dat %>% drop_na(RefSNP_id) # rm(list = ls()[ls() != "trans.dat"]);gc()
vroom::vroom_write(trans.dat, gsub("*Ch37*","Ch37_rsid", gwas_path)) # 保存


#' This is a toolkit for converting the `rsid` to `#CHR`/`POS` in the GWAS summary statistics
# RSID to CHR:POS ----
library(BSgenome)
library(dplyr)
library(magrittr)
library(data.table)
# example: my_rsids <- c("rs10458597", "rs12565286", "rs7553394")
my_rsids_df <- fread("/Users/leoarrow/project/software/ldsc-b站/eur_w_ld_chr/w_hm3.snplist") %>% select(1)
my_rsids <- as.vector(my_rsids_df[[1]])
id <- seq_along(my_rsids)
dat <- data.frame(my_rsids,id)

get_chrpos<-function(dat,snp_col=NULL,GRCh=NULL){
  rsids <- data.table::setDT(data.frame(snpsById(GRCh, my_rsids, ifnotfound = "drop"))) %>%
    dplyr::select(seqnames,pos,RefSNP_id) %>%
    dplyr::rename(CHR=seqnames,POS=pos,SNP=RefSNP_id)
  final_res<-merge(rsids,dat,by.x="SNP",by.y=snp_col,all.y=TRUE)
  message(paste0("Cannot find following snp:\n",final_res[is.na(final_res$CHR),]$SNP,"\n"))
  return(final_res)
}

# Load ref and run
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
dat <- get_chrpos(dat,
                  snp_col = "my_rsids", # colname for snp
                  GRCh = SNPlocs.Hsapiens.dbSNP155.GRCh37::SNPlocs.Hsapiens.dbSNP155.GRCh37)
dat %>% filter(CHR == 6 & POS >= 25000000 & POS <= 34000000) %>% View()

