#' This is a R script built for mbQTL data without RS_id
# CHR:POS to RSID ----
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {stop("Need Path")}
mb_path <- args[1] # mb_path for mbQTL path
# mb_path <- "/Volumes/T9/ref/mbqtl/473/GCST90032348_buildGRCh38.tsv.gz"

dat <- fread(mb_path)
# dat <- setDT(gwas2) # mannually define
# colnames(dat)


library("SNPlocs.Hsapiens.dbSNP155.GRCh38")
library(BSgenome)
# setDT(dat)
dat[, ranges := paste(`chromosome`, ":", base_pair_location, "-", base_pair_location, sep = "")]
# try
# dat_head <- head(dat)
# dat_head[, ranges := sprintf("%s:%s-%s", `#CHR`, POS, POS)]

snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
snp.res <- snpsByOverlaps(snps, GRanges(dat$ranges))

### merge data
library(stringr)
library(tidyverse)
snp.res.dt <- as.data.table(snp.res)
snp.res.dt[, ranges := stringr::str_c(seqnames, ":", pos, "-", pos)]
trans.dat <- merge(snp.res.dt, dat, by = "ranges") # colnames(trans.dat)
columns_to_remove <- c("strand", "alleles_as_ambig", "ranges", "seqnames", "pos"); trans.dat[, (columns_to_remove) := NULL]
trans.dat <- trans.dat %>% drop_na(RefSNP_id) # rm(list = ls()[ls() != "trans.dat"]);gc()
fwrite(trans.dat, mb_path, nThread =15)

