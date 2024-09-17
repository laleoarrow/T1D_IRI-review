library(ggplot2)
library(ggsci)
library(vroom)
library(tidyverse)
library(data.table)
library(dplyr)
library(tidyr)

# >>>>>>>>>>>>>>>>>>>>>> overall >>>>>>>>>>>>>>>>>>>>>> 
file_path_t1d <- "/Users/leoarrow/project/iridocyclitis/data/diabete/1_hla_finn/hla_summary_stats_T1D"
file_path_iri <- "/Users/leoarrow/project/iridocyclitis/data/iri_hla/hla_summary_stats_H7_IRIDOCYCLITIS"

data_t1d <- vroom(file_path_t1d)
data_iri <- vroom(file_path_iri)
merge_data <- merge(data_t1d, data_iri, by = "id", suffixes = c("_T1D", "_IRI"))
n <- nrow(merge_data); n 
  
# Both significant alleles & Alleles with consistent effects
sig_t1d <- data_t1d %>% filter(pval < 5e-10); dim(sig_t1d)
sig_iri <- data_iri %>% filter(pval < 5e-10); dim(sig_iri)
common_sig <- merge(sig_t1d, sig_iri, by = "id", suffixes = c("_T1D", "_IRI")) ; dim(common_sig)
consistent_direction <- common_sig %>% filter(beta_T1D * beta_IRI > 0) ; dim(consistent_direction)

# annotation
merge_data <- merge_data %>% mutate(`Both significance` = ifelse(id %in% common_sig$id, "1", "0"),
                                    `Consistent Direction` = ifelse(id %in% consistent_direction$id, "1", "0"),
                                    annotation = case_when(`Both significance` == "1" & `Consistent Direction` == "1" ~ "**",
                                                           `Both significance` == "1" & `Consistent Direction` == "0" ~ "*",
                                                           TRUE ~ ""))
# count "**" and "*"
counts <- merge_data %>%
  count(annotation) %>% set_names("Annotation", "N") %>%
  # filter(annotation %in% c("**", "*")) %>% 
  mutate(Percentage = sprintf("%.2f%%", (N / nrow(merge_data)) * 100)) %>% 
  mutate(Type = case_when(
    Annotation == "" ~ paste0("Not all significant (P > 5e-10)"),
    Annotation == "*" ~  paste0("Both significant (P < 5e-10)"),
    Annotation == "**" ~ "Both significance with consistent effect"
  )) %>% 
  select(Type, everything())

# change from wide to long
merge_data_long <- tidyr::pivot_longer(merge_data,
                                       cols = c(beta_T1D, sebeta_T1D, beta_IRI, sebeta_IRI, pval_T1D, pval_IRI),
                                       names_to = c(".value", "Phenotype"), names_sep = "_")
# alpha setting
merge_data_long <- merge_data_long %>%
  mutate(`-log10(P)` = -log10(pval),
         alpha_label = factor(
           case_when(
             `-log10(P)` >= -log10(5e-10) ~ ">5e-10",
             `-log10(P)` >= -log10(5e-8)  ~ ">5e-8",
             `-log10(P)` >= -log10(1e-5)  ~ ">1e-5",
             TRUE                         ~ "<1e-5"
           ), levels = c("<1e-5", ">1e-5", ">5e-8", ">5e-10"))
         )
head(merge_data_long, n = 4) %>% as.data.frame()
vertical_lines <- merge_data_long %>% filter(annotation %in% c("*","**"))

# illustration
p <- ggplot(merge_data_long, aes(x = id, y = beta, ymin = beta - 1.96*sebeta, ymax = beta + 1.96*sebeta)) +
  # geom_vline(data = vertical_lines, aes(xintercept = as.numeric(factor(id)), linetype = "dashed", color = "grey", alpha = 0.4)) +
  geom_pointrange(aes(color = Phenotype, alpha = alpha_label)) +
  scale_alpha_manual(values = c(0.25, 0.5, 0.75, 1.0), labels = c("<1e-5", ">1e-5", ">5e-8", ">5e-10")) +
  geom_hline(yintercept = 0) +
  # coord_flip() +
  theme_classic() +
  labs(
    x = "MHC Allele ID",
    y = "Beta Value",
    title = "Comparison of MHC Allele Effects on T1D and Iridocyclitis",
    alpha = expression("-log"[10] * italic("P")) # smart
  ) + 
  ggsci::scale_color_npg()+
  theme(
    plot.title = element_text(hjust = 0, size = 20, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.8),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.background = element_blank()
    ); p
ggsave("./figure/mhc/mhc_all.pdf", plot = p, width = 15, height = 6)

# add table
library(grid)
library(gridExtra)
table_theme <- ttheme_minimal(
  base_size = 12,
  colhead = list(
    padding = unit(c(10, 3, 10, 3), "mm"),
    fg_params = list(hjust = 0, x = 0.1, fontface = "bold"),
    bg_params = list(fill = "grey90", col = NA)
  ),
  core = list(
    padding = unit(c(10, 3, 10, 3), "mm"),
    fg_params = list(hjust = 0, x = 0.1),
    bg_params = list(fill = NA, col = NA)
  )
)
table_grob <- gridExtra::tableGrob(counts %>% select(-Annotation), rows = NULL, theme = table_theme)
grid.newpage()
grid.draw(table_grob)

# plot figure and table in it
pdf("./figure/mhc/mhc_all_table.pdf", width = 15, height = 6)
grid.newpage()
vp1 <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # main
vp2 <- viewport(width = 0.25, height = 0.2, x = 0.1, y = 0.9, just = c("left", "top"))
# 1
pushViewport(vp1)
print(p, newpage = FALSE)
popViewport()
# 2
pushViewport(vp2)
grid.draw(table_grob)
popViewport()
dev.off()

# >>>>>>>>>>>>>>>>>>>>>> both sig >>>>>>>>>>>>>>>>>>>>>> 
rm(list = ls()); gc()
file_path_t1d <- "/Users/leoarrow/project/iridocyclitis/data/diabete/1_hla_finn/hla_summary_stats_T1D"
file_path_iri <- "/Users/leoarrow/project/iridocyclitis/data/iri_hla/hla_summary_stats_H7_IRIDOCYCLITIS"

data_t1d <- vroom(file_path_t1d)
data_iri <- vroom(file_path_iri)
merge_data <- merge(data_t1d, data_iri, by = "id", suffixes = c("_T1D", "_IRI"))
n <- nrow(merge_data); n 

# 
sig_t1d <- data_t1d %>% filter(pval < 5e-10); dim(sig_t1d)
sig_iri <- data_iri %>% filter(pval < 5e-10); dim(sig_iri)
common_sig <- merge(sig_t1d, sig_iri, by = "id", suffixes = c("_T1D", "_IRI")) ; dim(common_sig)
n_both_sig <- nrow(common_sig); n_both_sig

# Separating alleles with consistent and inconsistent effects.
consistent_direction <- common_sig %>% filter(beta_T1D * beta_IRI > 0) %>% arrange(beta_T1D)
inconsistent_direction <- common_sig %>% filter(beta_T1D * beta_IRI <= 0) %>% arrange(beta_T1D)

# add direction
consistent_direction <- consistent_direction %>% mutate(direction = "Consistent") %>% mutate(order = row_number())
inconsistent_direction <- inconsistent_direction %>% mutate(direction = "Inconsistent") %>% mutate(order = row_number() + nrow(consistent_direction))

# merge
ordered_data <- bind_rows(consistent_direction, inconsistent_direction)

# wide to long
merge_data_long <- tidyr::pivot_longer(ordered_data,
                                       cols = c(beta_T1D, sebeta_T1D, beta_IRI, sebeta_IRI, pval_T1D, pval_IRI),
                                       names_to = c(".value", "Phenotype"), names_sep = "_")
merge_data_long <- merge_data_long %>% arrange(order)

p <- ggplot(merge_data_long, aes(x = reorder(id, -order), y = beta, ymin = beta - 1.96*sebeta, ymax = beta + 1.96*sebeta)) +
  # geom_vline(data = vertical_lines, aes(xintercept = as.numeric(factor(id)), linetype = "dashed", color = "grey", alpha = 0.4)) +
  geom_pointrange(aes(color = Phenotype)) +
  scale_alpha_manual(values = c(0.25, 0.5, 0.75, 1.0), labels = c("<1e-5", ">1e-5", ">5e-8", ">5e-10")) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  theme_bw() +
  labs(
    x = "MHC Allele ID",
    y = "Beta Value (95% CI)",
    title = "MHC Allele Effects Comparison",
    alpha = expression("-log"[10] * italic("P")) # smart
  ) + 
  # scale_color_manual(values = c("#F39B7FFF","#00A087FF"))+
  scale_color_manual(values = c("#00A087FF","#E64B35FF"))+
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
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
    # legend.background = element_blank()
  ); p
ggsave("./figure/mhc/mhc_sig.pdf", plot = p, width = 6, height = 6)

# plot.title = element_text(hjust = 0, size = 20, face = "bold"),
# axis.title.x = element_text(size = 16),
# axis.title.y = element_text(size = 16),
# axis.text = element_text(size = 14), 
# legend.position = c(0.8, 0.8),
# legend.title = element_text(size = 14, face = "bold"),
# legend.text = element_text(size = 12),
# legend.background = element_blank()







