# =================================================================
# 4-method integrated SUBTYPE figure (Supp): scPagwas + MAGMA + LDSC + scDRS
# 3 GWAS (HW, EUR, DIV) - 20 subtypes
# Color by combination of significant methods (raw P < 0.05)
# =================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

# Subtypes only: Excitatory + Inhibitory subtypes (14 total)
# Other cell types (Purkinje, Astro, Oligo, etc.) shown in MAJOR figure
CT_ORDER <- c("Ex-mix","Ex-L2.3","Ex-L2.4","Ex-L4.6","Ex-L5","Ex-L5.6","Ex-L6","Ex-NRGN",
              "In-CALM1","In-LAMP5","In-PVALB","In-SHANK2","In-SST","In-VIP")

# Robust cell-type normalizer:
#  spaces, underscores, hyphens, slashes -> periods (THEN squash multi-period to one)
norm_ct <- function(x) {
  y <- gsub("[ /]", ".", x)              # space, slash -> period
  y <- gsub("_", ".", y)                  # underscore -> period (handles In_PVALB)
  # Keep hyphens (Ex-L2.3 stays Ex-L2.3); but harmonize "Ex_L2/3" -> "Ex.L2.3"
  # For subtype we want "Ex-L2.4" style; need to re-hyphenate Ex / In prefix
  y <- gsub("^Ex\\.([LN])", "Ex-\\1", y)  # "Ex.L..." -> "Ex-L..."
  y <- gsub("^Ex\\.mix$", "Ex-mix", y)
  y <- gsub("^In\\.([CLPSV])", "In-\\1", y)  # "In.PVALB" -> "In-PVALB"
  y
}

read_magma <- function(path, gwas) {
  d <- read.table(path, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
  data.frame(GWAS = gwas, cell_type = norm_ct(d$VARIABLE), MAGMA_P = d$P)
}
read_ldsc <- function(path, gwas) {
  if (!file.exists(path)) return(NULL)
  d <- read.csv(path, stringsAsFactors = FALSE)
  data.frame(GWAS = gwas, cell_type = norm_ct(d$cell_type), LDSC_P = d$Coef_p)
}
read_scpagwas <- function(path, gwas) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  data.frame(GWAS = gwas, cell_type = norm_ct(d$celltype), scPagwas_P = d$pvalue)
}
read_scdrs <- function(path, gwas) {
  if (!file.exists(path)) return(NULL)
  d <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data.frame(GWAS = gwas, cell_type = norm_ct(d$group), scDRS_P = d$assoc_mcp)
}

mag <- bind_rows(
  read_magma("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/HW_final_anno_covar_v2.gsa.out",  "HW"),
  read_magma("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/EUR_final_anno_covar_v2.gsa.out", "EUR"),
  read_magma("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/DIV_final_anno_covar_v2.gsa.out", "DIV")
)
ldsc <- bind_rows(
  read_ldsc(file.path(DAT, "ldsc_v2_HW_subtype.csv"),  "HW"),
  read_ldsc(file.path(DAT, "ldsc_v2_EUR_subtype.csv"), "EUR")
)
scp <- bind_rows(
  read_scpagwas("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_HW_paperparams_20260615/cell_type_pvalues_paperparams_subtype.csv","HW"),
  read_scpagwas("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/cell_type_pvalues_v2_subtype.csv",          "EUR"),
  read_scpagwas("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/cell_type_pvalues_paperparams_subtype.csv","DIV")
)
scdrs <- bind_rows(
  read_scdrs("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/downstream_1k/HW.scdrs_group.final_anno",  "HW"),
  read_scdrs("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/downstream_1k/EUR.scdrs_group.final_anno", "EUR"),
  read_scdrs("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/downstream_1k/DIV.scdrs_group.final_anno", "DIV")
)

# Diagnostic: name mismatches
cat("MAGMA cell types (HW):", paste(sort(unique(subset(mag,GWAS=='HW')$cell_type)), collapse=", "), "\n\n")
cat("LDSC cell types (HW):",  paste(sort(unique(subset(ldsc,GWAS=='HW')$cell_type)), collapse=", "), "\n\n")
cat("scP cell types (HW):",   paste(sort(unique(subset(scp,GWAS=='HW')$cell_type)), collapse=", "), "\n\n")
cat("scDRS cell types (HW):", paste(sort(unique(subset(scdrs,GWAS=='HW')$cell_type)), collapse=", "), "\n\n")

df <- mag %>%
  full_join(ldsc,  by = c("GWAS","cell_type")) %>%
  full_join(scp,   by = c("GWAS","cell_type")) %>%
  full_join(scdrs, by = c("GWAS","cell_type")) %>%
  filter(cell_type %in% CT_ORDER)

classify_combo <- function(m, l, s, d) {
  flags <- c(MAGMA = !is.na(m) & m < 0.05,
             LDSC  = !is.na(l) & l < 0.05,
             scPagwas = !is.na(s) & s < 0.05,
             scDRS = !is.na(d) & d < 0.05)
  sig_names <- names(flags)[flags]
  if (length(sig_names) == 0)      "None"
  else if (length(sig_names) == 4) "All four"
  else if (length(sig_names) == 1) paste(sig_names, "only")
  else                             paste(sig_names, collapse = "+")
}

df <- df %>%
  rowwise() %>%
  mutate(group = classify_combo(MAGMA_P, LDSC_P, scPagwas_P, scDRS_P),
         n_sig = sum(c(MAGMA_P, LDSC_P, scPagwas_P, scDRS_P) < 0.05, na.rm = TRUE),
         n_avail = sum(!is.na(c(MAGMA_P, LDSC_P, scPagwas_P, scDRS_P))),
         mean_neglogP = mean(-log10(c(MAGMA_P, LDSC_P, scPagwas_P, scDRS_P)),
                              na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = rev(CT_ORDER)),
         GWAS = factor(GWAS, levels = c("HW","EUR","DIV"))) %>%
  filter(n_avail > 0)

write.csv(df, file.path(DAT, "four_method_subtype.csv"), row.names = FALSE)

# Palette
GROUP_COLORS <- c(
  "All four"                        = "#7a0a0a",
  "MAGMA+LDSC+scPagwas"             = "#cc3a36",
  "MAGMA+LDSC+scDRS"                = "#cc3a36",
  "MAGMA+scPagwas+scDRS"            = "#e85a55",
  "LDSC+scPagwas+scDRS"             = "#e85a55",
  "MAGMA+LDSC"                      = "#fa8072",
  "MAGMA+scDRS"                     = "#fa8072",
  "LDSC+scDRS"                      = "#fa8072",
  "MAGMA+scPagwas"                  = "#ee7e60",
  "LDSC+scPagwas"                   = "#ee7e60",
  "scPagwas+scDRS"                  = "#f4a261",
  "MAGMA only"                      = "#6cba73",
  "LDSC only"                       = "#4c79b6",
  "scDRS only"                      = "#e88f18",
  "scPagwas only"                   = "#b698c5",
  "None"                            = "#b7b7b7"
)

p <- ggplot(df, aes(x = mean_neglogP, y = cell_type, fill = group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", mean_neglogP)),
            hjust = -0.2, size = 2.5) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "gray50") +
  facet_wrap(~ GWAS, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = GROUP_COLORS,
                    name = "Significant\nmethod(s)\n(raw P < 0.05)",
                    drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = expression(Mean~-log[10](P)~"across available methods"),
       y = NULL,
       title = "Supplementary: 4-method subtype-level cell-type enrichment",
       subtitle = "scPagwas + MAGMA + LDSC + scDRS | 20 brain cell subtypes | DIV LDSC=NA") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "right",
        legend.text = element_text(size = 8))

ggsave(file.path(FIG, "four_method_HW_EUR_DIV_subtype.pdf"), p,
       width = 13, height = 7)
ggsave(file.path(FIG, "four_method_HW_EUR_DIV_subtype.png"), p,
       width = 13, height = 7, dpi = 200, bg = "white")
cat("\n========== SUMMARY (top 10 per GWAS) ==========\n")
for (g in c("HW","EUR","DIV")) {
  cat("\n--- ", g, " (top 10) ---\n")
  d <- df %>% filter(GWAS == g) %>% arrange(desc(mean_neglogP)) %>% head(10) %>%
       select(cell_type, MAGMA_P, LDSC_P, scPagwas_P, scDRS_P, group, mean_neglogP) %>%
       mutate(across(ends_with("_P"), ~ signif(.x, 3)),
              mean_neglogP = round(mean_neglogP, 2))
  print(d, row.names = FALSE)
}

cat("\nWrote:\n",
    file.path(FIG, "four_method_HW_EUR_DIV_subtype.pdf"), "\n",
    file.path(FIG, "four_method_HW_EUR_DIV_subtype.png"), "\n", sep = "")
