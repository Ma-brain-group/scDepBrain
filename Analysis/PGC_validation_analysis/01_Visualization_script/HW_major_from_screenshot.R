# =================================================================
# HW major 4-method dot chart, using hard-coded values from user's screenshot table
# Mixes: scPagwas (user paper), LDSC (user paper Fig 1a), MAGMA (our v2), scDRS (our 1000 ctrl)
# Same style as four_method_HW_EUR_DIV_major figure
# =================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

# ----- HW HARD-CODED P VALUES (from user screenshot) -----
df <- data.frame(
  cell_type = c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons",
                "Endothelial.cells","Microglia","OPCs","Astrocytes","Oligodendrocytes"),
  scPagwas_P = c(3.35e-09, 3.14e-06, 5.98e-04, 3.41e-01, 3.64e-01, 3.92e-01, 3.94e-01, 4.50e-01),
  scDRS_P    = c(2.00e-03, 9.99e-04, 9.99e-04, 9.82e-01, 9.78e-01, 2.43e-01, 4.63e-01, 1.00e+00),
  MAGMA_P    = c(7.33e-01, 7.48e-04, 1.23e-03, 1.30e-01, 1.74e-02, 5.62e-01, 4.97e-01, 9.88e-01),
  LDSC_P     = c(7.75e-02, 9.96e-07, 3.49e-03, 9.99e-01, 9.94e-01, 5.46e-02, 4.44e-01, 7.43e-01),
  stringsAsFactors = FALSE
)

CT_ORDER <- c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons",
              "Endothelial.cells","Oligodendrocytes","OPCs","Microglia","Astrocytes")

classify_combo <- function(m, l, s, d) {
  flags <- c(MAGMA = m < 0.05, LDSC = l < 0.05,
             scPagwas = s < 0.05, scDRS = d < 0.05)
  sigs <- names(flags)[flags]
  if (length(sigs) == 0)      "None"
  else if (length(sigs) == 4) "All four"
  else if (length(sigs) == 1) paste(sigs, "only")
  else                        paste(sigs, collapse = "+")
}

df <- df %>%
  rowwise() %>%
  mutate(group = classify_combo(MAGMA_P, LDSC_P, scPagwas_P, scDRS_P),
         n_sig = sum(c(MAGMA_P, LDSC_P, scPagwas_P, scDRS_P) < 0.05, na.rm = TRUE),
         mean_neglogP = mean(-log10(c(MAGMA_P, LDSC_P, scPagwas_P, scDRS_P)))) %>%
  ungroup() %>%
  mutate(cell_type = factor(cell_type, levels = rev(CT_ORDER)))

# Save table
write.csv(df, file.path(DAT, "HW_main_4method_from_screenshot.csv"), row.names = FALSE)
options(width = 200)
cat("\n========= HW 4-method (from screenshot) =========\n")
print(df %>% arrange(desc(mean_neglogP)) %>%
      mutate(across(ends_with("_P"), ~ signif(.x, 3))), row.names = FALSE)

# ----- Palette (same as four_method_HW_EUR_DIV) -----
GROUP_COLORS <- c(
  "All four"                = "#7a0a0a",
  "MAGMA+LDSC+scPagwas"     = "#cc3a36",
  "MAGMA+LDSC+scDRS"        = "#cc3a36",
  "MAGMA+scPagwas+scDRS"    = "#e85a55",
  "LDSC+scPagwas+scDRS"     = "#e85a55",
  "MAGMA+LDSC"              = "#fa8072",
  "MAGMA+scDRS"             = "#fa8072",
  "LDSC+scDRS"              = "#fa8072",
  "MAGMA+scPagwas"          = "#ee7e60",
  "LDSC+scPagwas"           = "#ee7e60",
  "scPagwas+scDRS"          = "#f4a261",
  "MAGMA only"              = "#6cba73",
  "LDSC only"               = "#4c79b6",
  "scDRS only"              = "#e88f18",
  "scPagwas only"           = "#b698c5",
  "None"                    = "#b7b7b7"
)

p <- ggplot(df, aes(x = mean_neglogP, y = cell_type, fill = group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", mean_neglogP)),
            hjust = -0.2, size = 3.5) +
  geom_text(aes(x = -0.05, label = sprintf("(%d/4)", n_sig)),
            hjust = 1, size = 3.0, color = "gray40") +
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "gray50") +
  scale_fill_manual(values = GROUP_COLORS,
                    name = "Significant\nmethod(s)\n(raw P < 0.05)",
                    drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.20))) +
  labs(x = expression(Mean~-log[10](P)~"across 4 methods"),
       y = NULL,
       title = "HW (Howard 2019) — 4-method cell-type enrichment",
       subtitle = "scPagwas + scDRS + MAGMA-CellTyping + LDSC-SEG | raw P < 0.05 per method") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 10))

ggsave(file.path(FIG, "HW_main_4method_major_v3.pdf"), p, width = 9, height = 5)
ggsave(file.path(FIG, "HW_main_4method_major_v3.png"), p,
       width = 9, height = 5, dpi = 200, bg = "white")
cat("\nWrote:\n",
    file.path(FIG, "HW_main_4method_major_v3.pdf"), "\n",
    file.path(FIG, "HW_main_4method_major_v3.png"), "\n", sep = "")
