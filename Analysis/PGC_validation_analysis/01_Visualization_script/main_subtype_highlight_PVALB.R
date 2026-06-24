# =================================================================
# MAIN FIGURE — subtype-resolved 4-method enrichment, In on top + Ex on bottom
# In-PVALB and Ex-L2/4 highlighted (bold + dark red labels, thick border)
# Inhibitory subtypes ordered first (top); within each, sorted by mean -log10(P)
# across HW/EUR/DIV (descending)
# =================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

# ---------- Load existing subtype data ----------
df <- read.csv(file.path(DAT, "four_method_subtype.csv"), stringsAsFactors = FALSE)
df <- df %>% mutate(class = ifelse(grepl("^Ex", cell_type), "Excitatory", "Inhibitory"))
df$GWAS <- factor(df$GWAS, levels = c("HW","EUR","DIV"))
df$class <- factor(df$class, levels = c("Inhibitory","Excitatory"))   # Inhibitory on top

# Rank within each class by mean -log10P across 3 GWAS, descending
rank_df <- df %>%
  group_by(cell_type, class) %>%
  summarise(across_gwas_mean = mean(mean_neglogP, na.rm = TRUE), .groups = "drop") %>%
  arrange(class, desc(across_gwas_mean))

# rev() so top rank appears at top of y-axis within each facet
ct_levels <- rev(rank_df$cell_type)
df$cell_type <- factor(df$cell_type, levels = ct_levels)

# Highlighted subtypes
HIGHLIGHT <- c("In-PVALB", "Ex-L2.4")

# ---------- Color palette (PGC Cell-2025 Figure 3 inspired,跟 current 主图配色一致) ----
GROUP_COLORS <- c(
  "All four"             = "#7a0a0a",
  "MAGMA+LDSC+scPagwas"  = "#cc3a36",
  "MAGMA+LDSC+scDRS"     = "#cc3a36",
  "MAGMA+scPagwas+scDRS" = "#e85a55",
  "LDSC+scPagwas+scDRS"  = "#e85a55",
  "MAGMA+LDSC"           = "#fa8072",
  "MAGMA+scDRS"          = "#fa8072",
  "LDSC+scDRS"           = "#fa8072",
  "MAGMA+scPagwas"       = "#ee7e60",
  "LDSC+scPagwas"        = "#ee7e60",
  "scPagwas+scDRS"       = "#f4a261",
  "MAGMA only"           = "#6cba73",
  "LDSC only"            = "#4c79b6",
  "scDRS only"           = "#e88f18",
  "scPagwas only"        = "#b698c5",
  "None"                 = "#b7b7b7"
)

# y-axis text styling: bold + red for highlighted, plain otherwise
y_face   <- ifelse(ct_levels %in% HIGHLIGHT, "bold", "plain")
y_color  <- ifelse(ct_levels %in% HIGHLIGHT, "#7a0a0a", "black")
y_size   <- ifelse(ct_levels %in% HIGHLIGHT, 11, 9)

# Border width per row (highlight In-PVALB + Ex-L2/4 with thicker)
df$is_highlight <- df$cell_type %in% HIGHLIGHT

# ---------- Plot ----------
p <- ggplot(df, aes(x = mean_neglogP, y = cell_type, fill = group)) +
  geom_col(aes(linewidth = is_highlight), color = "black", width = 0.78) +
  scale_linewidth_manual(values = c("FALSE" = 0.25, "TRUE" = 0.95), guide = "none") +
  geom_text(aes(label = sprintf("%.2f", mean_neglogP)),
            hjust = -0.15, size = 3.0, color = "black") +
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "gray45", linewidth = 0.4) +
  facet_grid(class ~ GWAS, scales = "free_y", space = "free_y",
             labeller = labeller(GWAS = c(
               "HW"  = "MDD 2019 EUR\n(Howard et al.)",
               "EUR" = "MDD 2025 EUR",
               "DIV" = "MDD 2025 trans-ancestry"))) +
  scale_fill_manual(values = GROUP_COLORS,
                    name = "Significant\nmethod(s)\n(raw P < 0.05)",
                    drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = expression(Mean~-log[10](P)~"across available methods"),
       y = NULL,
       title = "Subtype-resolved cell-type enrichment of MDD genetic risk",
       subtitle = "PVALB+ inhibitory and L2/4 cortical excitatory neurons converge across three GWAS (highlighted in dark red)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text.x = element_text(face = "bold", size = 11),
        strip.text.y = element_text(face = "bold", size = 12, angle = -90),
        strip.background.x = element_rect(fill = "gray95"),
        strip.background.y = element_rect(fill = "gray88"),
        axis.text.y = element_text(face = y_face, color = y_color, size = y_size),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 10, color = "gray35"),
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key.size = unit(0.35, "cm"))

ggsave(file.path(FIG, "Fig4_subtype_4method_PVALB_highlight.pdf"), p,
       width = 13, height = 7)
ggsave(file.path(FIG, "Fig4_subtype_4method_PVALB_highlight.png"), p,
       width = 13, height = 7, dpi = 300, bg = "white")

cat("\n========== Subtype ranking within class (mean across 3 GWAS) ==========\n")
options(width = 200)
print(rank_df, row.names = FALSE)

cat("\nWrote:\n",
    file.path(FIG, "Fig4_subtype_4method_PVALB_highlight.pdf"), "\n",
    file.path(FIG, "Fig4_subtype_4method_PVALB_highlight.png"), "\n", sep = "")
