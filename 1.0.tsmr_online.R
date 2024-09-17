# TSMR online pipeline
# Lu Ao, Chongqing, China
library("MendelianRandomization")
library("TwoSampleMR")
library("data.table")
library("vroom")
# https://api.opengwas.io/
# OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJsdWFvQHN0dS5jcW11LmVkdS5jbiIsImlhdCI6MTcyMDg2OTUwNiwiZXhwIjoxNzIyMDc5MTA2fQ.QSU5blzZexib-rOjfBV9e8wSSRSv5LeisjcbhO9uJQ5osmUYnVj-VRxdniL07C6_b3c9wzwS0gmbp6zGLw6dV0TniqdL32E-mMc1CnF5jiKext9Jc3i5dYmaMqjKigkFjZNlevl2Ri_pqJpt4u4ushPITcuNKznhu5ZUcCp0WDxe6DXRqaxgRZgv7fVzujgSscxMt58bF6yXCgwZ-acUuqPkom_aq03tBVXXgdimPK7btVP8eF09OXBjaNXE_436PA57xyiBwaRuSFhIh39FxgbGjF1FKEAwrThXFvehOh3t_kbkOX6BW8Rv-uYJXbBftIkaTQGffbMOlqDJlo5BIA"
# Sys.setenv(OPENGWAS_JWT=OPENGWAS_JWT)
# ieugwasr::user()
# system("curl cip.cc")

# Exposure----
gwas_id <- c("ebi-a-GCST90002232", # fasting glucose
             "ebi-a-GCST90002244", # Glycated hemoglobin levels
             "ebi-a-GCST90002238", # fasting insulin
             "ebi-a-GCST90002227") # Two-hour glucose
gwas <- extract_instruments(outcomes = gwas_id,
                            p1 = 5e-08,
                            clump = F,
                            p2 = 5e-08,
                            r2 = 0.001,
                            kb = 10000,
                            opengwas_jwt = ieugwasr::get_opengwas_jwt())
# gwas_c <- clump_data(gwas,
#                      clump_kb = 10000,
#                      clump_r2 = 0.001,
#                      clump_p1 = 1,
#                      clump_p2 = 1, pop = "EUR")
gwas_c <- clump_data_local(gwas, 
                           plink_bin = plinkbinr::get_plink_exe(),
                           bfile = "/Users/leoarrow/project/ref/1kg.v3/EUR"
                           )

iv <- gwas_c
iv <- iv %>% 
  mutate(R2 = (2*(beta.exposure)^2*(eaf.exposure*(1-eaf.exposure))) / ((se.exposure^2) * samplesize.exposure), # R2 = [2*beta^2*maf*(1-maf)]/[se^2*N]
         F = R2*(samplesize.exposure-2)/(1-R2) # F = (N-2)*(R2/(1-R2)) # i.e., K = 1 for one SNP
  ) %>% 
  filter(F > 10) # filter out SNPs with F < 10
iv$maf.exposure <- ifelse(iv$eaf.exposure < 0.5, iv$eaf.exposure, 1-iv$eaf.exposure)
table(iv$maf.exposure > 0.01) # should be all TRUE
iv <- iv %>% subset(iv$maf.exposure > 0.01)
iv$id.exposure <- iv$exposure
iv %>% vroom_write("./output/tsmr/bloodsugar/xxx.tsv")

# Outcome----
# read local outcome in 1.0.tsmr_di.R
source("./code/1.0.tsmr_packages.R")
out_3 <- format_outcome(iri_3, iv$SNP, N = "Neff")
out_3 <- out_3 %>% mutate(maf.outcome = ifelse(eaf.outcome<0.5, eaf.outcome, 1-eaf.outcome)) %>% subset(maf.outcome>0.01)
out = out_3
rm(out_1,out_2,out_3);gc()
# online
# outcome <- extract_outcome_data(
#   snps = iv$SNP,
#   outcomes = c('finn-b-H7_IRIDOACUTE', 'finn-b-H7_IRIDOCYCLITIS')
# )

# Harmonise data
dat_h <- harmonise_data(iv, out) %>% subset(mr_keep)
dat_h$id.exposure <- dat_h$exposure; dat_h$id.outcome <- dat_h$outcome 
dat_h %>% vroom_write("./output/tsmr/bloodsugar/tsmr_dat_h_bloodsugar_0717.tsv")

# mr
# unique exposure and outcome pair
outcome_order <- c("Fasting glucose || id:ebi-a-GCST90002232", "Fasting insulin || id:ebi-a-GCST90002238", "Glycated hemoglobin levels || id:ebi-a-GCST90002244", "Two-hour glucose || id:ebi-a-GCST90002227")
uni_pair <- dat_h %>% 
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
  res <- mr_one_pair(dat_h, exp, out, res_dir= "./output/tsmr/bloodsugar/", fig_dir="./figure/tsmr/bloodsugar/")
  return(res)
}, cl=cl)
stopCluster(cl);rm(nc);rm(cl);gc()
# The `results_list` is actually a matrix becase the fucntion return a list of two data.frames
ress <- lapply(results_list, function(x) x$res_pair) %>% do.call(rbind, .)
ress_presso <- lapply(results_list, function(x) x$res_pair_presso) %>% do.call(rbind, .)
vroom_write(ress, "./output/tsmr/bloodsugar/tsmr_bloodsugar_tsmrres_0717.txt")
vroom_write(ress_presso, "./output/tsmr/bloodsugar/tsmr_bloodsugar_PRESSO_0717.txt")