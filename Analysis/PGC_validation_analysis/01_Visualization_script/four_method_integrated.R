# =================================================================
# 4-method integrated figure: scPagwas + MAGMA + LDSC + scDRS
# 3 GWAS (HW, EUR, DIV) - 8 major cell types
# Color by number/combination of significant methods (raw P < 0.05)
# Robust to missing data (NA cells stay gray)
# =================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

CT_ORDER <- c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons",
              "Endothelial.cells","Oligodendrocytes","OPCs","Microglia","Astrocytes")

norm_ct <- function(x) gsub("[ _]", ".", x)

# ---------- MAGMA gene-covar ----------
read_magma <- function(path, gwas) {
  d <- read.table(path, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
  data.frame(GWAS = gwas, cell_type = norm_ct(d$VARIABLE), MAGMA_P = d$P)
}

# ---------- LDSC v2 ----------
read_ldsc <- function(path, gwas) {
  if (!file.exists(path)) return(NULL)
  d <- read.csv(path, stringsAsFactors = FALSE)
  data.frame(GWAS = gwas, cell_type = norm_ct(d$cell_type), LDSC_P = d$Coef_p)
}

# ---------- scPagwas ----------
read_scpagwas <- function(path, gwas) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  data.frame(GWAS = gwas, cell_type = norm_ct(d$celltype), scPagwas_P = d$pvalue)
}

# ---------- scDRS ----------
read_scdrs <- function(path, gwas) {
  if (!file.exists(path)) return(NULL)
  d <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data.frame(GWAS = gwas, cell_type = norm_ct(d$group), scDRS_P = d$assoc_mcp)
}

# Load all available
mag <- bind_rows(
  read_magma("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/HW_anno_covar_v2.gsa.out",  "HW"),
  read_magma("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/EUR_anno_covar_v2.gsa.out", "EUR"),
  read_magma("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/DIV_anno_covar_v2.gsa.out", "DIV")
)
ldsc <- bind_rows(
  read_ldsc(file.path(DAT, "ldsc_v2_HW_major.csv"),  "HW"),
  read_ldsc(file.path(DAT, "ldsc_v2_EUR_major.csv"), "EUR")
  # DIV LDSC not run (trans-anc not appropriate)
)
scp <- bind_rows(
  read_scpagwas("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_HW_paperparams_20260615/cell_type_pvalues_paperparams.csv","HW"),
  read_scpagwas("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/cell_type_pvalues_v2.csv",          "EUR"),
  read_scpagwas("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/cell_type_pvalues_paperparams.csv","DIV")
)
scdrs <- bind_rows(
  read_scdrs("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/downstream_1k/HW.scdrs_group.anno",  "HW"),
  read_scdrs("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/downstream_1k/EUR.scdrs_group.anno", "EUR"),
  read_scdrs("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/downstream_1k/DIV.scdrs_group.anno", "DIV")
)

# Merge all
df <- mag %>%
  full_join(ldsc,  by = c("GWAS","cell_type")) %>%
  full_join(scp,   by = c("GWAS","cell_type")) %>%
  full_join(scdrs, by = c("GWAS","cell_type")) %>%
  filter(cell_type %in% CT_ORDER)

# -------- HW OVERRIDE: paper-consistent values (user screenshot) --------
# scPagwas + LDSC use the user's published HW Figure 1a values;
# MAGMA + scDRS retained from our v2 / 1000-ctrl pipeline (identical to read above).
hw_override <- data.frame(
  GWAS = "HW",
  cell_type = c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons",
                "Endothelial.cells","Microglia","OPCs","Astrocytes","Oligodendrocytes"),
  MAGMA_P    = c(7.33e-01, 7.48e-04, 1.23e-03, 1.30e-01, 1.74e-02, 5.62e-01, 4.97e-01, 9.88e-01),
  LDSC_P     = c(7.75e-02, 9.96e-07, 3.49e-03, 9.99e-01, 9.94e-01, 5.46e-02, 4.44e-01, 7.43e-01),
  scPagwas_P = c(3.35e-09, 3.14e-06, 5.98e-04, 3.41e-01, 3.64e-01, 3.92e-01, 3.94e-01, 4.50e-01),
  scDRS_P    = c(2.00e-03, 9.99e-04, 9.99e-04, 9.82e-01, 9.78e-01, 2.43e-01, 4.63e-01, 1.00e+00),
  stringsAsFactors = FALSE
)
df <- df %>% filter(GWAS != "HW") %>% bind_rows(hw_override)

# Per row: count significant + assemble label + mean -log10P
classify_n <- function(m, l, s, d) {
  flags <- c(MAGMA = !is.na(m) & m < 0.05,
             LDSC  = !is.na(l) & l < 0.05,
             scPagwas = !is.na(s) & s < 0.05,
             scDRS = !is.na(d) & d < 0.05)
  sig_names <- names(flags)[flags]
  n_sig <- length(sig_names)
  n_avail <- sum(!is.na(c(m, l, s, d)))
  list(n_sig = n_sig, n_avail = n_avail, sigs = paste(sig_names, collapse = ", "))
}

# classify: combination of significant methods
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
  mutate(
    info = list(classify_n(MAGMA_P, LDSC_P, scPagwas_P, scDRS_P)),
    n_sig = info$n_sig,
    n_avail = info$n_avail,
    sig_methods = info$sigs,
    group = classify_combo(MAGMA_P, LDSC_P, scPagwas_P, scDRS_P),
    mean_neglogP = mean(-log10(c(MAGMA_P, LDSC_P, scPagwas_P, scDRS_P)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(-info) %>%
  mutate(
    cell_type = factor(cell_type, levels = rev(CT_ORDER)),
    GWAS = factor(GWAS, levels = c("HW","EUR","DIV"))
  )

# Drop rows where everything is NA (cell_type doesn't appear in any method)
df <- df %>% filter(n_avail > 0)

write.csv(df, file.path(DAT, "four_method_integrated.csv"), row.names = FALSE)
cat("\n========= 4-method integrated table =========\n")
print(df %>% arrange(GWAS, desc(mean_neglogP)) %>%
      select(GWAS, cell_type, MAGMA_P, LDSC_P, scPagwas_P, scDRS_P,
             n_sig, n_avail, sig_methods, mean_neglogP) %>%
      mutate(across(ends_with("_P"), ~ signif(.x, 3))),
      row.names = FALSE)

# Palette (PGC Cell-2025 Figure 3 style, extended for 4 methods)
GROUP_COLORS <- c(
  # ⭐ Strongest: all 4 methods
  "All four"                        = "#7a0a0a",  # dark red

  # 3 methods
  "MAGMA+LDSC+scPagwas"             = "#cc3a36",  # red
  "MAGMA+LDSC+scDRS"                = "#cc3a36",
  "MAGMA+scPagwas+scDRS"            = "#e85a55",
  "LDSC+scPagwas+scDRS"             = "#e85a55",

  # 2 methods
  "MAGMA+LDSC"                      = "#fa8072",  # salmon
  "MAGMA+scDRS"                     = "#fa8072",
  "LDSC+scDRS"                      = "#fa8072",
  "MAGMA+scPagwas"                  = "#ee7e60",
  "LDSC+scPagwas"                  = "#ee7e60",
  "scPagwas+scDRS"                  = "#f4a261",

  # 1 method
  "MAGMA only"                      = "#6cba73",  # green
  "LDSC only"                       = "#4c79b6",  # blue
  "scDRS only"                      = "#e88f18",  # orange
  "scPagwas only"                   = "#b698c5",  # purple

  # None
  "None"                            = "#b7b7b7"   # gray
)

# Plot
p <- ggplot(df, aes(x = mean_neglogP, y = cell_type, fill = group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", mean_neglogP)),
            hjust = -0.2, size = 3.0) +
  geom_text(aes(x = -0.05, label = sprintf("(%d/%d)", n_sig, n_avail)),
            hjust = 1, size = 2.5, color = "gray40") +
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "gray50") +
  facet_wrap(~ GWAS, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = GROUP_COLORS,
                    name = "Significant\nmethod(s)\n(raw P < 0.05)",
                    drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.20))) +
  labs(x = expression(Mean~-log[10](P)~"across available methods"),
       y = NULL,
       title = "4-method integrated cell-type enrichment",
       subtitle = "scPagwas + MAGMA + LDSC + scDRS | raw P < 0.05 | (#sig / #available methods)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom")

ggsave(file.path(FIG, "four_method_HW_EUR_DIV_major.pdf"), p, width = 12, height = 5)
ggsave(file.path(FIG, "four_method_HW_EUR_DIV_major.png"), p, width = 12, height = 5,
       dpi = 200, bg = "white")
cat("\nWrote:\n",
    file.path(FIG, "four_method_HW_EUR_DIV_major.pdf"), "\n",
    file.path(FIG, "four_method_HW_EUR_DIV_major.png"), "\n", sep = "")
