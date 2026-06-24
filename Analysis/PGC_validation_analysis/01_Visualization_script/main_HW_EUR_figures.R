# =================================================================
# Main figure: HW + EUR, 8 major cell types
# Supplementary:  HW + EUR, 20 subtypes
# PGC Cell-2025 Figure 3 style: multi-color by significant-method pattern
# raw P < 0.05 threshold
# =================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

# Group palette (PGC Cell-2025 Figure 3 inspired)
GROUP_COLORS <- c(
  "All three"        = "#cc3a36",
  "MAGMA+LDSC"       = "#fa8072",
  "MAGMA+scPagwas"   = "#ee7e60",
  "LDSC+scPagwas"    = "#f4a261",
  "MAGMA only"       = "#6cba73",
  "LDSC only"        = "#4c79b6",
  "scPagwas only"    = "#b698c5",
  "None"             = "#b7b7b7"
)
GROUP_LEVELS <- names(GROUP_COLORS)

classify <- function(m, l, s) {
  flags <- c(MAGMA = !is.na(m) & m < 0.05,
             LDSC  = !is.na(l) & l < 0.05,
             scPagwas = !is.na(s) & s < 0.05)
  sigs  <- names(flags)[flags]
  if (length(sigs) == 0)      "None"
  else if (length(sigs) == 3) "All three"
  else if (length(sigs) == 2) paste(sigs, collapse = "+")
  else                        paste(sigs, "only")
}

# =================================================================
# Helper: read & merge for one (GWAS, level)
# =================================================================
read_magma <- function(path) {
  d <- read.table(path, header = TRUE, comment.char = "#",
                  stringsAsFactors = FALSE)
  data.frame(cell_type = d$VARIABLE, MAGMA_P = d$P)
}
read_ldsc <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  data.frame(cell_type = d$cell_type, LDSC_P = d$Coef_p)
}
read_scpagwas <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  data.frame(cell_type = d$celltype, scPagwas_P = d$pvalue)
}

# Canonical name normalizer: spaces / underscores -> dots
norm_ct <- function(x) gsub("[ _]", ".", x)

build_panel <- function(gwas_tag, magma_path, ldsc_path, scp_path, ct_order) {
  M <- read_magma(magma_path); M$cell_type <- norm_ct(M$cell_type)
  L <- read_ldsc(ldsc_path);   L$cell_type <- norm_ct(L$cell_type)
  S <- read_scpagwas(scp_path);S$cell_type <- norm_ct(S$cell_type)

  df <- M %>% left_join(L, by = "cell_type") %>%
              left_join(S, by = "cell_type")
  # Restrict to canonical order
  df <- df[df$cell_type %in% ct_order, ]
  # Guarantee all expected types appear
  miss <- setdiff(ct_order, df$cell_type)
  if (length(miss)) {
    df <- bind_rows(df,
      data.frame(cell_type = miss, MAGMA_P = NA, LDSC_P = NA, scPagwas_P = NA))
  }
  df <- df %>%
    rowwise() %>%
    mutate(group = classify(MAGMA_P, LDSC_P, scPagwas_P),
           mean_neglogP = mean(-log10(c(MAGMA_P, LDSC_P, scPagwas_P)),
                                na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(GWAS = gwas_tag,
           cell_type = factor(cell_type, levels = rev(ct_order)),
           group = factor(group, levels = GROUP_LEVELS))
  df
}

# =================================================================
# MAJOR — 8 cell types
# =================================================================
CT_MAJOR <- c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons",
              "Endothelial.cells","Oligodendrocytes","OPCs",
              "Microglia","Astrocytes")

major_HW <- build_panel("HW",
  "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/HW_anno_covar_v2.gsa.out",
  file.path(DAT, "ldsc_v2_HW_major.csv"),
  "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_HW_paperparams_20260615/cell_type_pvalues_paperparams.csv",
  CT_MAJOR)
major_EUR <- build_panel("EUR",
  "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/EUR_anno_covar_v2.gsa.out",
  file.path(DAT, "ldsc_v2_EUR_major.csv"),
  "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/cell_type_pvalues_v2.csv",
  CT_MAJOR)

major_df <- bind_rows(major_HW, major_EUR) %>%
  mutate(GWAS = factor(GWAS, levels = c("HW","EUR")))

write.csv(major_df, file.path(DAT, "main_HW_EUR_major.csv"), row.names = FALSE)

p_major <- ggplot(major_df, aes(x = mean_neglogP, y = cell_type, fill = group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", mean_neglogP)),
            hjust = -0.2, size = 3.1) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "gray50") +
  facet_wrap(~ GWAS, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = GROUP_COLORS, name = "Significance\n(raw P < 0.05)",
                    drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = expression(Mean~-log[10](P)~"across methods"),
       y = NULL,
       title = "Integrated cell-type enrichment — scPagwas + MAGMA + LDSC",
       subtitle = "PGC Cell-2025 Figure 3 style. Raw P < 0.05 per method.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 9))

ggsave(file.path(FIG, "main_HW_EUR_major.pdf"), p_major,
       width = 10, height = 4.8)
ggsave(file.path(FIG, "main_HW_EUR_major.png"), p_major,
       width = 10, height = 4.8, dpi = 200, bg = "white")

cat("\n=========== MAIN (major, 8 cell types) ===========\n")
print(major_df %>% arrange(GWAS, desc(mean_neglogP)), row.names = FALSE)

# =================================================================
# SUBTYPE — 20 subtypes (Supplementary)
# =================================================================
CT_SUBTYPE <- c("Ex-mix","Ex-L2.3","Ex-L2.4","Ex-L4.6","Ex-L5","Ex-L5.6","Ex-L6","Ex-NRGN",
                "In-CALM1","In-LAMP5","In-PVALB","In-SHANK2","In-SST","In-VIP",
                "Purkinje.neurons","Endothelial.cells","Oligodendrocytes","OPCs",
                "Microglia","Astrocytes")

subtype_HW <- build_panel("HW",
  "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/HW_final_anno_covar_v2.gsa.out",
  file.path(DAT, "ldsc_v2_HW_subtype.csv"),
  "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_HW_paperparams_20260615/cell_type_pvalues_paperparams_subtype.csv",
  CT_SUBTYPE)
subtype_EUR <- build_panel("EUR",
  "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/EUR_final_anno_covar_v2.gsa.out",
  file.path(DAT, "ldsc_v2_EUR_subtype.csv"),
  "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/cell_type_pvalues_v2_subtype.csv",
  CT_SUBTYPE)

subtype_df <- bind_rows(subtype_HW, subtype_EUR) %>%
  mutate(GWAS = factor(GWAS, levels = c("HW","EUR")))

write.csv(subtype_df, file.path(DAT, "supp_HW_EUR_subtype.csv"), row.names = FALSE)

p_subtype <- ggplot(subtype_df, aes(x = mean_neglogP, y = cell_type, fill = group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", mean_neglogP)),
            hjust = -0.2, size = 2.8) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "gray50") +
  facet_wrap(~ GWAS, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = GROUP_COLORS, name = "Significance\n(raw P < 0.05)",
                    drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = expression(Mean~-log[10](P)~"across methods"),
       y = NULL,
       title = "Supplementary: subtype-level cell-type enrichment (20 subtypes)",
       subtitle = "PGC Cell-2025 Figure 3 style. Raw P < 0.05 per method.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 9))

ggsave(file.path(FIG, "supp_HW_EUR_subtype.pdf"), p_subtype,
       width = 11, height = 7)
ggsave(file.path(FIG, "supp_HW_EUR_subtype.png"), p_subtype,
       width = 11, height = 7, dpi = 200, bg = "white")

cat("\n=========== SUPP (subtype, 20 types) — top 10 each =============\n")
for (g in c("HW","EUR")) {
  cat("\n--- ", g, " ---\n")
  d <- subtype_df %>% filter(GWAS == g) %>% arrange(desc(mean_neglogP))
  print(head(d, 10), row.names = FALSE)
}

cat("\n\nWrote:\n",
    file.path(FIG, "main_HW_EUR_major.pdf"), "\n",
    file.path(FIG, "main_HW_EUR_major.png"), "\n",
    file.path(FIG, "supp_HW_EUR_subtype.pdf"), "\n",
    file.path(FIG, "supp_HW_EUR_subtype.png"), "\n", sep = "")
