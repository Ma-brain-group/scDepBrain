# =================================================================
# HW (Howard 2019) — 4-method dotchart per subtype
# Mimics screenshot style: 4 horizontal panels (scPagwas / scDRS /
#   MAGMA_CellTyping / LDSC-SEG), cream background, lollipop dots
# Two figures: Inhibitory (6 subtypes) + Excitatory (8 subtypes)
# Subtype order: by scPagwas -log10(P) descending (matching paper style)
# =================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

# ---------- Load 4-method subtype data ----------
df <- read.csv(file.path(DAT, "four_method_subtype.csv"), stringsAsFactors = FALSE)
# Use HW only
df_hw <- df %>% filter(GWAS == "HW")
cat("HW subtype rows:", nrow(df_hw), "\n")
print(df_hw[, c("cell_type","MAGMA_P","LDSC_P","scPagwas_P","scDRS_P")])

# ---------- Reshape to long format ----------
reshape_long <- function(d) {
  d %>%
    select(cell_type, MAGMA_P, LDSC_P, scPagwas_P, scDRS_P) %>%
    pivot_longer(cols = -cell_type, names_to = "method", values_to = "P") %>%
    mutate(
      method = recode(method,
        "MAGMA_P"   = "MAGMA_CellTyping",
        "LDSC_P"    = "LDSC-SEG",
        "scPagwas_P"= "scPagwas",
        "scDRS_P"   = "scDRS"),
      neglogP = ifelse(is.na(P) | P <= 0, 0, -log10(P))
    )
}

# Excitatory + Inhibitory subtypes
EX_SUBTYPES <- c("Ex-L2.4","Ex-L2.3","Ex-L4.6","Ex-L5","Ex-L5.6","Ex-L6","Ex-NRGN","Ex-mix")
IN_SUBTYPES <- c("In-PVALB","In-LAMP5","In-SST","In-VIP","In-SHANK2","In-CALM1")

# Color palette (consistent across panels)
EX_COLORS <- c(
  "Ex-L2.4"  = "#4F9D4A",     # green (highlight)
  "Ex-L2.3"  = "#D58A6C",     # salmon
  "Ex-L4.6"  = "#7DBE6F",     # light green
  "Ex-L5"    = "#A69333",     # olive
  "Ex-L5.6"  = "#3E8E59",     # darker green
  "Ex-L6"    = "#E54B45",     # red
  "Ex-NRGN"  = "#E66B58",     # coral
  "Ex-mix"   = "#9DA9B1"      # gray-blue
)
IN_COLORS <- c(
  "In-PVALB"  = "#8e4c99",    # purple (highlight, matches "Inhibitory neurons" in main)
  "In-LAMP5"  = "#E68C45",    # orange
  "In-SST"    = "#5BAEAE",    # teal
  "In-VIP"    = "#A69333",    # olive
  "In-SHANK2" = "#E54B45",    # red
  "In-CALM1"  = "#3F8C8C"     # blue-teal
)

# ---------- Plot function ----------
mk_dotchart <- function(d_long, subtypes, palette, sort_method = "scPagwas",
                        ymax = 5.5, panel_title_prefix = "") {
  d_long <- d_long %>% filter(cell_type %in% subtypes)
  # Order subtypes by scPagwas -log10(P) descending
  ord <- d_long %>% filter(method == sort_method) %>%
    arrange(desc(neglogP)) %>% pull(cell_type)
  d_long$cell_type <- factor(d_long$cell_type, levels = rev(ord))
  d_long$method <- factor(d_long$method,
                           levels = c("scPagwas","scDRS","MAGMA_CellTyping","LDSC-SEG"))

  ggplot(d_long, aes(x = neglogP, y = cell_type, color = cell_type)) +
    geom_segment(aes(x = 0, xend = neglogP, yend = cell_type),
                 linewidth = 0.5, color = "gray60") +
    geom_point(size = 7) +
    geom_text(aes(label = sprintf("%.3g", neglogP)),
              hjust = -0.3, size = 3.0, color = "black") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "#3a8fb7", linewidth = 0.5) +
    facet_wrap(~ method, nrow = 1, scales = "free_x") +
    scale_color_manual(values = palette) +
    scale_x_continuous(limits = c(0, ymax),
                       expand = expansion(mult = c(0.02, 0.20))) +
    labs(x = expression(-log[10](P)), y = NULL,
         title = paste0(panel_title_prefix,
                        " (Howard et al. 2019, MDD 2019 EUR)")) +
    theme_minimal(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "#fffef5", color = NA),
      panel.background = element_rect(fill = "#fffef5", color = NA),
      strip.background = element_rect(fill = "#f5e4d4", color = "gray50",
                                       linewidth = 0.4),
      strip.text = element_text(face = "bold", size = 13),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray85"),
      axis.text.y = element_text(face = "bold", size = 11, color = "black"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      panel.spacing.x = unit(0.6, "lines")
    )
}

df_long <- reshape_long(df_hw)

# Highlight In-PVALB and Ex-L2.4 with bold red y-axis text — apply after
modify_highlight <- function(p, highlight) {
  ct_levels <- levels(p$data$cell_type)
  y_face   <- ifelse(ct_levels %in% highlight, "bold.italic", "bold")
  y_color  <- ifelse(ct_levels %in% highlight, "#7a0a0a", "black")
  y_size   <- ifelse(ct_levels %in% highlight, 12, 10)
  p + theme(axis.text.y = element_text(face = y_face, color = y_color, size = y_size))
}

# Compute max neglogP for axis scaling
ex_max <- df_long %>% filter(cell_type %in% EX_SUBTYPES) %>% pull(neglogP) %>% max(na.rm=TRUE)
in_max <- df_long %>% filter(cell_type %in% IN_SUBTYPES) %>% pull(neglogP) %>% max(na.rm=TRUE)
cat("\nEx max -log10(P):", round(ex_max, 2), "\nIn max -log10(P):", round(in_max, 2), "\n")

# ---------- INHIBITORY (6 subtypes, place on top) ----------
p_in <- mk_dotchart(df_long, IN_SUBTYPES, IN_COLORS,
                     sort_method = "scPagwas",
                     ymax = max(in_max, 5) + 1,
                     panel_title_prefix = "Inhibitory neuron subtypes")
p_in <- modify_highlight(p_in, "In-PVALB")

ggsave(file.path(FIG, "Fig_HW_subtype_inhibitory_4method.pdf"), p_in,
       width = 14, height = 4)
ggsave(file.path(FIG, "Fig_HW_subtype_inhibitory_4method.png"), p_in,
       width = 14, height = 4, dpi = 300, bg = "white")

# ---------- EXCITATORY (8 subtypes, place on bottom) ----------
p_ex <- mk_dotchart(df_long, EX_SUBTYPES, EX_COLORS,
                     sort_method = "scPagwas",
                     ymax = max(ex_max, 5) + 1,
                     panel_title_prefix = "Excitatory neuron subtypes")
p_ex <- modify_highlight(p_ex, "Ex-L2.4")

ggsave(file.path(FIG, "Fig_HW_subtype_excitatory_4method.pdf"), p_ex,
       width = 14, height = 5)
ggsave(file.path(FIG, "Fig_HW_subtype_excitatory_4method.png"), p_ex,
       width = 14, height = 5, dpi = 300, bg = "white")

# ---------- COMBINED (Inhibitory top + Excitatory bottom in one figure) ----------
p_combined <- (p_in + theme(plot.title = element_blank())) /
              (p_ex + theme(plot.title = element_blank())) +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(
    title = "Subtype-resolved cell-type enrichment of MDD genetic risk — Howard et al. 2019 (MDD 2019 EUR)",
    subtitle = "Inhibitory subtypes (top) and excitatory subtypes (bottom); In-PVALB and Ex-L2/4 highlighted (bold dark-red labels)",
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10, color = "gray35"))
  )

ggsave(file.path(FIG, "Fig_HW_subtype_combined_4method.pdf"), p_combined,
       width = 14, height = 8.5)
ggsave(file.path(FIG, "Fig_HW_subtype_combined_4method.png"), p_combined,
       width = 14, height = 8.5, dpi = 300, bg = "white")

cat("\n========== HW Subtype Summary (sorted by scPagwas -log10(P)) ==========\n")
options(width = 200)
print(df_long %>%
  filter(cell_type %in% c(EX_SUBTYPES, IN_SUBTYPES)) %>%
  mutate(class = ifelse(grepl("^Ex", cell_type), "Excitatory", "Inhibitory"),
         P_value = signif(P, 3),
         neglogP = round(neglogP, 3)) %>%
  select(class, cell_type, method, P_value, neglogP) %>%
  arrange(class, desc(neglogP)),
  row.names = FALSE)

cat("\nWrote:\n",
    file.path(FIG, "Fig_HW_subtype_inhibitory_4method.pdf"), "\n",
    file.path(FIG, "Fig_HW_subtype_excitatory_4method.pdf"), "\n",
    file.path(FIG, "Fig_HW_subtype_combined_4method.pdf"), "\n", sep = "")
