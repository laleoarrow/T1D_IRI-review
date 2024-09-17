#' This code aims to perform Placo analysis----------------
#' Source code and Toy example by Placo author were attached after mail codes
# Preparation ----
library(vroom);library(data.table);library(tidyverse)
gwas1_path <- "~/project/iridocyclitis/data/diabete/1/GCST90014023_buildhg19.tsv" # t1d1
gwas1 <- vroom(gwas1_path) %>% 
 select(SNP, CHR, BP, A1, A2, P, BETA, SE, EAF) %>% 
 set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF") %>% 
 drop_na(CHR, POS) %>%
 mutate(Phenotype = "T1D (Chiou et al.)", MAF = ifelse(EAF<0.5, EAF, 1-EAF)) %>% 
 filter(MAF > 0.01)

gwas2_path <- "~/project/iridocyclitis/data/iri_finn/ukb_finn/iri_fbmeta_hg19.tsv" # iri3
gwas2 <- vroom(gwas2_path) %>% 
  mutate(N_Cohort = case_when(is.na(FINNGEN_beta) & !is.na(UKBB_beta) ~ "UKB",
                              !is.na(FINNGEN_beta) & is.na(UKBB_beta) ~ "Finn",
                              !is.na(FINNGEN_beta) & !is.na(UKBB_beta) ~ "Finn+UKB",
                              TRUE ~ "Unknown"),
         EAF = case_when(N_Cohort == "UKB" ~ UKBB_af_alt,
                         N_Cohort == "Finn" ~ FINNGEN_af_alt,
                         N_Cohort == "Finn+UKB" ~ (FINNGEN_af_alt*(8016+390647)+UKBB_af_alt*(2304+413959))/(8016+390647+2304+413959),
                         TRUE ~ NA)) %>% 
  select(SNP, CHR, BP, ALT, REF, all_inv_var_meta_p, all_inv_var_meta_beta, all_inv_var_meta_sebeta, EAF) %>% 
  set_names("SNP", "CHR", "POS", "A1", "A2", "P", "BETA", "SE", "EAF") %>%
  mutate(Phenotype = "Iridocyclitis (FinnGen+UKB)", MAF = ifelse(EAF<0.5, EAF, 1-EAF)) %>% 
  filter(MAF > 0.01)


# harmonize these two datasets using TwoSampleMR
# pre_reduce the row numbers in each gwas if possible 
# e.g., delete snp with too large Z^2; or locate the intersect SNP in advance
# library(TwoSampleMR)
gwas1 <- gwas1 %>% # give exposure name
 rename(SNP = SNP,
        beta.exposure = BETA,
        se.exposure = SE,
        effect_allele.exposure = A1,
        other_allele.exposure = A2,
        eaf.exposure = EAF,
        exposure = Phenotype) %>% 
  mutate(id.exposure = outcome)
gwas2 <- gwas2 %>% # give outcome name
  rename(SNP = SNP,
        beta.outcome = BETA,
        se.outcome = SE,
        effect_allele.outcome = A1,
        other_allele.outcome = A2,
        eaf.outcome = EAF,
        outcome = Phenotype) %>% 
  mutate(id.outcome = outcome)
# reduce the row numbers if possible
.pre_filter_z <- function(df, beta, se){
  n1 = nrow(df);df <- df %>% mutate(Z = {{beta}}/{{se}}) %>% filter(Z*Z < 80)
  n2 = nrow(df);message(paste0("Delete <", n2-n1, "> SNP with too large Z^2"))
  return(df)
}
.filter_min_p <- function(df) {
  df <- df %>% # mark who is duplicated
    mutate(dup = duplicated(SNP, fromLast = TRUE) | duplicated(SNP, fromLast = FALSE))
  df_dup <- df %>%  # filter duplicated SNP and find the ones with min p
    filter(dup) %>%
    group_by(SNP) %>%
    filter(P == min(P)) %>% # P name needs to be adjusted
    dplyr::slice(1) %>% # keep the 1 when p is the same
    ungroup()
  df_non_dup <- df %>% # non-duplicated SNP
    filter(!dup)
  rbind(df_non_dup, df_dup) %>% # brilliant!
    select(-dup) # This way only check the duplicated ones
}
intersect_snp <- intersect(gwas1$SNP, gwas2$SNP); message(paste0("A total of <", length(intersect_snp), "> SNPs"))
gwas1 <- gwas1 %>% .pre_filter_z(., beta.exposure, se.exposure) %>% filter(SNP %in% intersect_snp) %>% .filter_min_p()
gwas2 <- gwas2 %>% .pre_filter_z(., beta.outcome, se.outcome) %>% filter(SNP %in% intersect_snp) %>% .filter_min_p()
rm(intersect_snp); gc()

# harmonise_data
h_gwas <- TwoSampleMR::harmonise_data(gwas1, gwas2) # For debug: h_gwas_bak <- h_gwas
h_gwas <- h_gwas %>%  
  subset(mr_keep) %>% 
  mutate(Z1 = beta.exposure/se.exposure, Z2 = beta.outcome/se.outcome,
         ZZ1 = Z.y, ZZ2 = Z.x, # just in case tsmr change the source code one day
         P1 = P.y, P2 = P.x) %>% 
  filter(!(Z1*Z1 > 80 | Z2*Z2 > 80)) # just in case
# !! This one is really tricky, after the harmonise via tsmr, *.y is exposure and *.x is the outcome!
# !! Clue in souce code of TwoSampleMR::harmonise_data: res.tab <- merge(outcome_dat, exposure_dat, by = "SNP")
# !! So make sure the Z and P is corresponds to one !! Be carefull about this !!
if (any(h_gwas$Z1 != h_gwas$ZZ1)) { # so just in case tsmr change the source code one day here
  stop("Error: Column Z1 and ZZ1 contain different values! Check if Z1 should be Z.y!!")
} else {
  h_gwas <- h_gwas %>% select(-ZZ1, -ZZ2)
}


# build Z and P matrix for placo
Z.matrix <- h_gwas %>% 
 select(Z1, Z2) %>% 
 as.matrix()
P.matrix <- h_gwas %>%
 select(P1, P2) %>%
 as.matrix()

# PLACO ----
# rm(gwas1);rm(gwas2);gc()
require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.1.1.R?raw=TRUE") # credits goes to the author of PLACO
## If the traits are dependent or correlated, we suggest decorrelating the Z-scores (only once), then apply Steps 1 and 2 on the decorrelated Z-scores
## Step 0a: Obtain the correlation matrix of Z-scores----
R <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)
## Step 0b: Decorrelate the matrix of Z-scores----
k <- 2 # number of traits
"%^%" <- function(x, pow) # function for raising matrix to any power
 with(eigen(x), vectors %*% (values^pow * t(vectors)))
Z.matrix.decor <- Z.matrix %*% (R %^% (-0.5))
colnames(Z.matrix.decor) <- paste("Z",1:k,sep="")
## Step 1: Obtain the variance parameter estimates (only once)----
VarZ <- var.placo(Z.matrix.decor, P.matrix, p.threshold=1e-4)
## Step 2: Apply test of pleiotropy for each variant----
k <- 2 # number of traits
p <- dim(Z.matrix.decor)[1] # number of variants
### multicore placo----
library(parallel); library(pbapply)
nc <- detectCores()-6
cl <- makeCluster(nc)
clusterExport(cl, varlist = c(
 "VarZ", "Z.matrix.decor",
 "var.placo", "placo",".p.bessel",'.pdfx')
)
# out <- sapply(1:p, function(i) placo(Z=Z.matrix[i,], VarZ=VarZ))
out <- pbsapply(1:p, function(i) placo(Z=Z.matrix.decor[i,], VarZ=VarZ), cl=cl)
stopCluster(cl);rm(cl);gc()

# get placo results
t_out = t(out) %>% as.data.frame()
h_gwas$T.placo <-  t_out$T.placo
h_gwas$p.placo <- t_out$p.placo # format(object.size(h_gwas), units = "GB")
h_gwas <- h_gwas %>% 
  select(SNP, CHR.x, POS.x, effect_allele.exposure,other_allele.exposure, 
         beta.exposure, beta.outcome, se.exposure, se.outcome,
         eaf.exposure, eaf.outcome, outcome, exposure,
         Z1, Z2, P1, P2, T.placo, p.placo) %>% 
  rename(CHR=CHR.x, POS=POS.x, A1=effect_allele.exposure, A2=other_allele.exposure)
  
# save placo results
outpath = "/Users/leoarrow/project/iridocyclitis/output/placo"
fwrite(x = h_gwas, file = file.path(outpath,'t1d1_iri3_placo.tsv.gz'), sep="\t", nThread = 16)
# plcao significant snps
h_gwas_sig <- h_gwas %>% 
  filter(p.placo < 5e-8)
fwrite(x = h_gwas_sig, file = file.path(outpath,'t1d1_iri3_placo_sig.tsv'), sep="\t", nThread = 16)

# Clump placo data-----
library(data.table);library(tidyverse);library(TwoSampleMR);library(ieugwasr)
placo <- fread("/Users/leoarrow/project/iridocyclitis/output/placo/t1d1_iri3_placo.tsv.gz", nThread = 16)
placo_lead <- placo %>% 
 filter(p.placo < 5e-8) %>% 
 arrange(p.placo)

clump_data_local <- function(dat, clump_kb = 1000, clump_r2 = 0.2) {
 # Note: `clump_kb = 1000, clump_r2 = 0.2` is choice of placo author
 require(ieugwasr); require(plinkbinr)
 message(" - Clumping data locally")
 dat1 <- ieugwasr::ld_clump(
  dat = dplyr::tibble(rsid=dat$SNP, pval=dat$p.placo, id=dat$id.exposure), # Need indicate SNP here
  clump_kb = clump_kb, 
  clump_r2 = clump_r2,
  clump_p = 5e-8,
  bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR",
  plink_bin = plinkbinr::get_plink_exe()
 )
 dat2<- subset(dat, SNP %in% dat1$rsid)
 return(dat2)
}
clumped <- clump_data_local(placo_lead)
clumped_strict <- clump_data_local(placo_lead, clump_kb = 10000, clump_r2 = 0.001)
fwrite(clumped, file = "./output/placo/t1d1_iri3_placo_1000kb0.2_clumped.tsv", sep="\t", nThread = 16)
fwrite(clumped_strict, file = "./output/placo/t1d1_iri3_placo_10000kb0.001_clumped.tsv", sep="\t", nThread = 16)

# Manhattan for placo----
library(data.table);library(tidyverse)
placo <- data.table::fread("/Users/leoarrow/project/iridocyclitis/output/placo/t1d1_iri3_placo.tsv.gz", nThread = 16)
data_draw <- placo %>% # colnames(placo)
 mutate(`Effect Direction` = ifelse(sign(beta.exposure) == sign(beta.outcome), "Consistent effects", "Opposite effects"),
        category = case_when(P1 < 5e-8 & P2 < 5e-8 ~ "P1/2 < 5e-8", # both sig
                             P1 < 5e-8 & P2 >= 5e-8 & P2 <= 1e-5 ~ "P1/2 < 5e-8, P2/1 in [5e-8, 1e-5]", # one sig one suggestive
                             P2 < 5e-8 & P1 >= 5e-8 & P1 <= 1e-5 ~ "P1/2 < 5e-8, P2/1 in [5e-8, 1e-5]", # one sig one suggestive
                             P1 >= 5e-8 & P1 <= 1e-5 & P2 >= 5e-8 & P2 <= 1e-5 ~ "P1/2 in [5e-8, 1e-5]", # both suggestive
                             TRUE ~ "P1/2 > 1e-5")) %>% # both null
 select(CHR = CHR, BP = POS, P = p.placo,
        SNP, `Effect Direction`, category) # opposite effect or not
data_draw$category <- factor(data_draw$category, levels = c("P1/2 < 5e-8", # both sig
                                                            "P1/2 < 5e-8, P2/1 in [5e-8, 1e-5]", # one sig one suggestive
                                                            "P1/2 in [5e-8, 1e-5]", # both suggestive
                                                            "P1/2 > 1e-5")) # both null
# color
nature_colors <- ggsci::pal_npg()(6)[c(1, 2, 3, 5, 4, 6)]
```
"#E64B35FF" "#4DBBD5FF" "#00A087FF" "#F39B7FFF" "#3C5488FF" "#8491B4FF"
```
chr_colors <- rep(c(nature_colors[5], nature_colors[6]), length.out = 23); names(chr_colors) <- as.character(1:23) # base color
data_draw$chr_color <- chr_colors[as.character(data_draw$CHR)] # base color
data_draw$final_color <- case_when(
 data_draw$category == "P1/2 < 5e-8" ~ "#E64B35FF",
 data_draw$category == "P1/2 < 5e-8, P2/1 in [5e-8, 1e-5]" ~ "#00A087FF",
 data_draw$category == "P1/2 in [5e-8, 1e-5]" ~ "#4DBBD5FF",
 TRUE ~ data_draw$chr_color
)
data_draw$final_category <- case_when(
 data_draw$final_color == "#E64B35FF" ~ "P1/2 < 5e-8",
 data_draw$final_color == "#00A087FF" ~ "P1/2 < 5e-8, P2/1 in [5e-8, 1e-5]",
 data_draw$final_color == "#4DBBD5FF" ~ "P1/2 in [5e-8, 1e-5]",
 data_draw$final_color == nature_colors[5] ~ "P1/2 > 1e-5",
 data_draw$final_color == nature_colors[6] ~ "P2/1 > 1e-5"
)
color_labels <- setNames(c("#E64B35FF","#00A087FF","#4DBBD5FF", nature_colors[5],nature_colors[6]),
                         c("P1/2 < 5e-8", # label在后
                           "P1/2 < 5e-8, P2/1 in [5e-8, 1e-5]",
                           "P1/2 in [5e-8, 1e-5]",
                           "P1/2 > 1e-5",
                           "P2/1 > 1e-5"))
# legend annotation
cat_percent <- data_draw %>%
 group_by(category) %>%
 summarise(count = n()) %>%
 mutate(percent = count / sum(count) * 100,
        annotation = paste0(category, ": ",count, " (", round(percent, 4), "%)"))


chr_len <- data_draw %>%
 group_by(CHR) %>% 
 summarise(chr_len = max(BP))
chr_pos <- chr_len  %>%
 mutate(total = cumsum(chr_len)-chr_len) %>%
 select(-chr_len)
Snp_pos <- chr_pos %>%
 left_join(data_draw, ., by="CHR") %>% # sweet note: data_draw added in here
 arrange(CHR, BP) %>%
 mutate(BPcum = BP+total)
X_axis <-  Snp_pos %>%
 group_by(CHR) %>% 
 summarize(center = (max(BPcum)+min(BPcum))/2) 


adjust_dat <- Snp_pos %>% filter(P<5e-8)
library(ggplot2); library(ggrastr); library(ggsci)
x1 = ggplot(adjust_dat, aes(x = BPcum, y = -log10(P), 
                            color = final_category,
                            shape = `Effect Direction`)) +
 rasterise(geom_point(size = 3, alpha = 1), dpi = 500) +
 scale_shape_manual(values = c("Consistent effects" = 19, "Opposite effects" = 17), name = "Effect Direction") +
 scale_color_manual(values = color_labels, name = "SNP Category") +
 scale_x_continuous(label = X_axis$CHR, breaks= X_axis$center) +
 geom_hline(yintercept = c(-log10(1e-05), -log10(5e-08)), color = c('grey', 'black'), linetype = c("dashed", "dashed")) + # http://www.sthda.com/english/wiki/ggplot2-line-types-how-to-change-line-types-of-a-graph-in-r-software
 theme_classic() +
 labs(title = "PLACO Manhattan Plot", x = "Chromosome Position", y = "-log10(P-Value)") +
 theme(plot.title = element_text(hjust = 0, size = 20, face = "bold"),
       axis.title.x = element_text(size = 16), 
       axis.title.y = element_text(size = 16), 
       axis.text = element_text(size = 14),
       legend.position = c(0.8, 0.8),
       legend.title = element_text(size = 14, face = "bold"),
       legend.text = element_text(size = 12),
       legend.background = element_blank())
ggsave("./figure/manhatton/placo_t1d1_iri3.pdf", plot = x1, width = 15, height = 6)  

# Annotation using FUMA
snp_annotated <- data.table::fread("/Users/leoarrow/project/iridocyclitis/output/fuma/placo/FUMA_job493590/[snps].txt")
lead <- data.table::fread("/Users/leoarrow/project/iridocyclitis/output/fuma/placo/FUMA_job493590/【leadSNPs】.txt") %>% 
  left_join(snp_annotated %>% select(rsID, nearestGene), by = c("rsID" = "rsID"))
IndSigSNPs <- data.table::fread("/Users/leoarrow/project/iridocyclitis/output/fuma/placo/FUMA_job493590/【IndSigSNPs】.txt") %>% 
  left_join(snp_annotated %>% select(rsID, nearestGene), by = c("rsID" = "rsID"))
# lead %>% fwrite("/Users/leoarrow/project/iridocyclitis/output/fuma/placo/FUMA_job493590/【lead】addGene.txt", sep = "\t")
# IndSigSNPs %>% fwrite("/Users/leoarrow/project/iridocyclitis/output/fuma/placo/FUMA_job493590/【IndSigSNPs】addGene.txt", sep = "\t")
IndSigSNPs <- IndSigSNPs %>% mutate(lead = ifelse(rsID %in% lead$rsID, "Lead", "IndSig")) %>% 
  select(rsID, nearestGene, lead)
IndSigSNPs <- IndSigSNPs %>% left_join(Snp_pos %>% select(SNP, BPcum, P), by = c("rsID" = "SNP"))
IndSigSNPs <- IndSigSNPs %>% left_join(Snp_pos %>% select(SNP, final_category), by = c("rsID" = "SNP"))
IndSigSNPs <- IndSigSNPs %>% left_join(Snp_pos %>% select(SNP, `Effect Direction`), by = c("rsID" = "SNP"))
IndSigSNPs$nearestGene = ifelse(IndSigSNPs$lead=="Lead", paste0(IndSigSNPs$nearestGene,"*"), IndSigSNPs$nearestGene)
# ggplot with annotation
x2 = ggplot(Snp_pos, aes(x = BPcum, y = -log10(P), color = final_category, shape = `Effect Direction`)) +
  rasterise(geom_point(size = 3, alpha = 1), dpi = 500) +
  ggrepel::geom_text_repel(data = IndSigSNPs, aes(x=BPcum, y=-log10(P), label=nearestGene),
                           fontface = "italic", box.padding = 0.5, point.padding = 0.5,
                           force_pull = 0,
                           force = 1,
                           max.overlaps = Inf,
                           angle = 90, 
                           direction = "x",
                           nudge_x = 0, 
                           nudge_y = 5,  
                           size = 4,
                           hjust = 0,
                           segment.color = "black", color = "black") +
  scale_shape_manual(values = c("Consistent effects" = 19, "Opposite effects" = 17), name = "Effect Direction") +
  scale_color_manual(values = color_labels, name = "SNP Category") +
  scale_x_continuous(label = X_axis$CHR, breaks= X_axis$center) +
  geom_hline(yintercept = c(-log10(1e-05), -log10(5e-08)), color = c('grey', 'black'), linetype = c("dashed", "dashed")) + # http://www.sthda.com/english/wiki/ggplot2-line-types-how-to-change-line-types-of-a-graph-in-r-software
  theme_classic() +
  labs(title = "PLACO Manhattan Plot", x = "Chromosome", y = "-log10(PLACO P)") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.background = element_blank())
ggsave("./figure/manhatton/placo_t1d1_iri3_anno_bak.pdf", plot = x2, width = 15, height = 6) 



# Source code of placo; In case network in R died; Here is the source code for the PLACO function ----
.pdfx<-function(x) besselK(x=abs(x),nu=0)/pi
.p.bessel<-function(z, varz, AbsTol=1e-13){
 p1<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[1])),Inf, abs.tol=AbsTol)$value)
 p2<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[2])),Inf, abs.tol=AbsTol)$value)
 p0<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]),Inf, abs.tol=AbsTol)$value)
 pval.compnull<-p1+p2-p0
 return(pval.compnull)
}
var.placo<-function(Z.matrix, P.matrix, p.threshold=1e-4){
 # Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits 
 # Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
 # p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
 # checks
 k<-ncol(Z.matrix)
 if(k!=2) stop("This method is meant for 2 traits only. Columns correspond to traits.")
 ZP<-cbind(Z.matrix,P.matrix)
 ZP<-na.omit(ZP)
 
 rows.alt<-which(ZP[,3]<p.threshold & ZP[,4]<p.threshold)
 if(length(rows.alt)>0){
  ZP<-ZP[-rows.alt,]
  if(nrow(ZP)==0) stop(paste("No 'null' variant left at p-value threshold",p.threshold))
  if(nrow(ZP)<30) warning(paste("Too few 'null' variants at p-value threshold",p.threshold))
 }
 varz<-diag(var(ZP[,c(1,2)]))
 return(varz)
}

#---------------- Function for estimating correlation matrix of the Z's
cor.pearson<-function(Z.matrix, P.matrix, p.threshold=1e-4){
 # Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits 
 # Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
 # p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
 # checks
 k<-ncol(Z.matrix)
 if(k!=2) stop("This method is meant for 2 traits only.")
 # estimating correlation
 row.exclude<-which( apply(P.matrix, MARGIN = 1, function(x) any(x < p.threshold)) == TRUE )
 if(length(row.exclude)>0) Z.matrix<-Z.matrix[-row.exclude,]
 R<-cor(Z.matrix)
 return(R)
}

placo<-function(Z, VarZ, AbsTol=.Machine$double.eps^0.8){
 # Z: vector of Z-scores of size k=2 (i.e., collection of Z-scores of a particular SNP for k=2 traits)
 # VarZ: vector of variances of Z-scores (covariance assumed 0; so need to be independent traits)
 # AbsTol: absolute tolerance (accuracy paramater) for numerical integration.
 # checks				
 k<-length(Z)		
 if(k!=2) stop("This method is meant for 2 traits only.")
 if(length(VarZ)!=k) stop("Provide variance estimates for 2 traits as obtained using var.placo() function.")
 
 # test of pleiotropy: PLACO
 pvalue.b=.p.bessel(z=Z, varz=VarZ, AbsTol=AbsTol)
 return(list(T.placo=prod(Z), p.placo=pvalue.b))
}


# Toy illustration by placo author ------------
set.seed(1)
## For an example, let's first simulate a toy set of GWAS summary 
## statistics on 2 uncorrelated traits and 1000 variants
require(MASS)
k <- 2
p <- 1000
Z.matrix <- mvrnorm(n=p, mu=rep(0,k), Sigma=diag(1,k))
P.matrix <- matrix(NA, nrow=p, ncol=k)
for(j in 1:k){
 P.matrix[,j] <- sapply(1:nrow(Z.matrix), 
                        function(i) pchisq(Z.matrix[i,j]^2,df=1,ncp=0,lower.tail=F))
}	
colnames(Z.matrix) <- paste("Z",1:k,sep="")
colnames(P.matrix) <- paste("P",1:k,sep="")

## Steps to implementing PLACO
# Step 1: Obtain the variance parameter estimates (only once)
VarZ <- var.placo(Z.matrix, P.matrix, p.threshold=1e-4)
# Step 2: Apply test of pleiotropy for each variant
library(parallel)
library(pbapply)
nc <- detectCores()-1
cl <- makeCluster(nc)
clusterExport(cl, varlist = c(ls(),".p.bessel",'.pdfx'))
system.time(
 out <- pbsapply(1:p, function(i) placo(Z=Z.matrix[i,], VarZ=VarZ), cl=cl)
)
stopCluster(cl)
# Check the output for say variant 100
dim(out)
out[,100]$T.placo
out[,100]$p.placo

## If the traits are dependent or correlated, we suggest
## decorrelating the Z-scores (only once), then apply Steps 1 and 2
## on the decorrelated Z-scores

# Step 0a: Obtain the correlation matrix of Z-scores
R <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)
# Step 0b: Decorrelate the matrix of Z-scores
# function for raising matrix to any power
"%^%" <- function(x, pow)
 with(eigen(x), vectors %*% (values^pow * t(vectors)))
Z.matrix.decor <- Z.matrix %*% (R %^% (-0.5))
colnames(Z.matrix.decor) <- paste("Z",1:k,sep="")










