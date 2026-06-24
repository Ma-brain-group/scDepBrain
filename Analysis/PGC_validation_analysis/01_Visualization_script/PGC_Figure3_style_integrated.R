# =================================================================
# PGC Figure 3 style integrated plot
# scPagwas + MAGMA + LDSC across EUR + DIV (+ HW for completeness)
# Color = which methods significant at FDR < 0.05
# Bar  = mean(-log10(P)) across AVAILABLE methods
# DIV  has no LDSC (trans-ancestry / EUR LD ref mismatch)
# =================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617"
FIG <- file.path(OUT, "figures")
DAT <- file.path(OUT, "data")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# ------------------- 8 major cell types canonical order ---------------------
CT_ORDER <- c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons",
              "Endothelial.cells","Oligodendrocytes","OPCs","Microglia","Astrocytes")

# =============================================================
# 1. Load MAGMA gene-covar (Bryois CellTyping standard) — 3 GWAS
# =============================================================
parse_magma_covar <- function(path, gwas) {
  d <- read.table(path, header = TRUE, comment.char = "#",
                  stringsAsFactors = FALSE)
  data.frame(GWAS = gwas,
             cell_type = gsub("_", ".", d$VARIABLE),
             MAGMA_P = d$P,
             MAGMA_z = sign(d$BETA) * qnorm(1 - d$P / 2))
}
M_HW  <- parse_magma_covar("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/HW_anno_covar_v2.gsa.out",  "HW")
M_EUR <- parse_magma_covar("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/EUR_anno_covar_v2.gsa.out", "EUR")
M_DIV <- parse_magma_covar("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/DIV_anno_covar_v2.gsa.out", "DIV")
magma_df <- rbind(M_HW, M_EUR, M_DIV)

# =============================================================
# 2. Load LDSC v2 (with all_gene control) — HW + EUR only
# =============================================================
parse_ldsc <- function(path, gwas) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  data.frame(GWAS = gwas,
             cell_type = d$cell_type,
             LDSC_P = d$Coef_p,
             LDSC_z = d$Coef_z)
}
L_HW  <- parse_ldsc(file.path(DAT, "ldsc_v2_HW_major.csv"),  "HW")
L_EUR <- parse_ldsc(file.path(DAT, "ldsc_v2_EUR_major.csv"), "EUR")
ldsc_df <- rbind(L_HW, L_EUR)

# =============================================================
# 3. Load scPagwas Merge_celltype_p — all 3 GWAS (paper-params)
# =============================================================
parse_scp <- function(path, gwas) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  ct <- gsub(" ", ".", d$celltype)
  data.frame(GWAS = gwas, cell_type = ct, scPagwas_P = d$pvalue)
}
S_HW  <- parse_scp("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_HW_paperparams_20260615/cell_type_pvalues_paperparams.csv", "HW")
S_EUR <- parse_scp("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/cell_type_pvalues_v2.csv",         "EUR")
S_DIV <- parse_scp("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/cell_type_pvalues_paperparams.csv","DIV")
scp_df <- rbind(S_HW, S_EUR, S_DIV)

# =============================================================
# 4. Merge: GWAS × cell_type × {MAGMA, LDSC, scPagwas}
# =============================================================
df <- magma_df %>%
  left_join(ldsc_df, by = c("GWAS","cell_type")) %>%
  left_join(scp_df,  by = c("GWAS","cell_type"))

# Compute FDR per GWAS per method (BH over 8 cell types)
df <- df %>%
  group_by(GWAS) %>%
  mutate(MAGMA_FDR    = p.adjust(MAGMA_P,   "BH"),
         LDSC_FDR     = p.adjust(LDSC_P,    "BH"),
         scPagwas_FDR = p.adjust(scPagwas_P,"BH")) %>%
  ungroup() %>%
  filter(cell_type %in% CT_ORDER)

# -------------------- classification --------------------
classify <- function(magma_sig, ldsc_sig, scp_sig) {
  sigs <- c(MAGMA=magma_sig, LDSC=ldsc_sig, scPagwas=scp_sig)
  sigs <- sigs[!is.na(sigs) & sigs]   # significant methods
  if (length(sigs) == 0)        "None"
  else if (length(sigs) >= 3)   "All three"
  else if (length(sigs) == 2)   paste(names(sigs), collapse = "+")
  else                          paste(names(sigs)[1], "only")
}

df <- df %>%
  rowwise() %>%
  mutate(
    MAGMA_sig    = ifelse(!is.na(MAGMA_P)    & MAGMA_P    < 0.05, TRUE, FALSE),
    LDSC_sig     = ifelse(!is.na(LDSC_P)     & LDSC_P     < 0.05, TRUE, FALSE),
    scPagwas_sig = ifelse(!is.na(scPagwas_P) & scPagwas_P < 0.05, TRUE, FALSE),
    group = classify(MAGMA_sig, LDSC_sig, scPagwas_sig),
    mean_neglogP = mean(c(-log10(MAGMA_P), -log10(LDSC_P), -log10(scPagwas_P)), na.rm = TRUE),
    n_methods    = sum(!is.na(c(MAGMA_P, LDSC_P, scPagwas_P)))
  ) %>%
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = rev(CT_ORDER)),
         GWAS = factor(GWAS, levels = c("HW","EUR","DIV")))

# Save merged table
write.csv(df, file.path(DAT, "integrated_3method_3GWAS.csv"), row.names = FALSE)
cat("\n=========== Integrated Table ===========\n")
print(df %>% select(GWAS, cell_type, MAGMA_P, LDSC_P, scPagwas_P,
                    MAGMA_FDR, LDSC_FDR, scPagwas_FDR,
                    group, mean_neglogP, n_methods))

# -------------------- PGC palette --------------------
GROUP_COLORS <- c(
  "All three"          = "#cc3a36",   # red — strongest evidence
  "MAGMA+LDSC"         = "#fa8072",   # salmon — PGC's "Both"
  "MAGMA+scPagwas"     = "#ee7e60",
  "LDSC+scPagwas"      = "#f4a261",
  "MAGMA only"         = "#6cba73",   # green — MAGMA only
  "LDSC only"          = "#4c79b6",   # blue  — LDSC only
  "scPagwas only"      = "#b698c5",   # purple
  "None"               = "#b7b7b7"    # gray  — none
)

# =============================================================
# 5. Plot — PGC Figure 3 style
# =============================================================
p_main <- ggplot(df, aes(x = mean_neglogP, y = cell_type, fill = group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", mean_neglogP)),
            hjust = -0.15, size = 3) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "gray50") +
  facet_wrap(~ GWAS, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = GROUP_COLORS, name = "Significance\n(raw P < 0.05)") +
  labs(x = expression(Mean~-log[10](P)~"across methods"),
       y = NULL,
       title = "Integrated cell-type enrichment: scPagwas + MAGMA + LDSC",
       subtitle = "PGC-Cell-2025 Figure 3 style. Raw P < 0.05. DIV: LDSC unavailable (EUR LD ref mismatch).") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 9))

ggsave(file.path(FIG, "PGCFig3_integrated_major.pdf"),
       p_main, width = 11, height = 4.5)
ggsave(file.path(FIG, "PGCFig3_integrated_major.png"),
       p_main, width = 11, height = 4.5, dpi = 200, bg = "white")
cat("\nWrote:", file.path(FIG, "PGCFig3_integrated_major.pdf"), "\n")
cat("Wrote:", file.path(FIG, "PGCFig3_integrated_major.png"), "\n")

cat("\nDone.\n")
