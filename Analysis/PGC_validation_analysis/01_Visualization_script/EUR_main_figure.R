# =================================================================
# EUR main figure — PGC MDD 2025 EUR (largest, latest)
# 3 methods: scPagwas + MAGMA-CellTyping + LDSC-SEG
# 8 major cell types
# Color: None / 1 method / >= 2 methods (raw P < 0.05)
# =================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

CT_ORDER <- c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons",
              "Endothelial.cells","Oligodendrocytes","OPCs","Microglia","Astrocytes")

m <- read.table("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/EUR_anno_covar_v2.gsa.out",
                header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
M <- data.frame(cell_type = gsub("_", ".", m$VARIABLE), MAGMA_P = m$P)

L <- read.csv(file.path(DAT, "ldsc_v2_EUR_major.csv"), stringsAsFactors = FALSE)
L <- L[, c("cell_type","Coef_p")]; names(L)[2] <- "LDSC_P"

s <- read.csv("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/cell_type_pvalues_v2.csv",
              stringsAsFactors = FALSE)
S <- data.frame(cell_type = gsub(" ", ".", s$celltype), scPagwas_P = s$pvalue)

df <- M %>% left_join(L, by = "cell_type") %>% left_join(S, by = "cell_type") %>%
  filter(cell_type %in% CT_ORDER) %>%
  rowwise() %>%
  mutate(n_sig = sum(c(MAGMA_P, LDSC_P, scPagwas_P) < 0.05, na.rm = TRUE),
         mean_neglogP = mean(-log10(c(MAGMA_P, LDSC_P, scPagwas_P)), na.rm = TRUE),
         sig_methods = paste(c("MAGMA","LDSC","scPagwas")[
           c(MAGMA_P, LDSC_P, scPagwas_P) < 0.05], collapse = ", ")) %>%
  ungroup() %>%
  mutate(group = factor(case_when(n_sig >= 2 ~ "≥2 methods",
                                   n_sig == 1 ~ "1 method",
                                   TRUE       ~ "None"),
                        levels = c("≥2 methods","1 method","None")),
         cell_type = factor(cell_type, levels = rev(CT_ORDER)))

out_tbl <- df %>% arrange(desc(mean_neglogP))
write.csv(out_tbl, file.path(DAT, "EUR_main_integrated.csv"), row.names = FALSE)
cat("\n========= EUR — final integrated (raw P < 0.05) =========\n")
print(out_tbl, row.names = FALSE)

GROUP_COLORS <- c("≥2 methods" = "#cc3a36",   # red — strongest
                  "1 method"   = "#f4a261",   # orange — partial
                  "None"       = "#b7b7b7")   # gray

pA <- ggplot(df, aes(x = mean_neglogP, y = cell_type, fill = group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", mean_neglogP)),
            hjust = -0.2, size = 3.5) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "gray50") +
  scale_fill_manual(values = GROUP_COLORS, name = "Methods\nsignificant\n(raw P < 0.05)",
                    drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = expression(Mean~-log[10](P)~"across scPagwas + MAGMA + LDSC"),
       y = NULL,
       title = "(A) Integrated cell-type enrichment — PGC MDD 2025 EUR",
       subtitle = "8 major brain cell types. Color: number of methods with raw P < 0.05.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "right")

long <- df %>%
  pivot_longer(c(MAGMA_P, LDSC_P, scPagwas_P),
               names_to = "method", values_to = "P") %>%
  mutate(method = factor(sub("_P$","", method),
                         levels = c("scPagwas","MAGMA","LDSC")),
         sig = ifelse(P < 0.05, "Significant", "NS"),
         neglogP = -log10(P))

pB <- ggplot(long, aes(x = neglogP, y = cell_type, color = sig)) +
  geom_segment(aes(x = 0, xend = neglogP, yend = cell_type),
               color = "gray80", linewidth = 1) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "gray50") +
  geom_text(aes(label = sprintf("%.2f", neglogP)),
            hjust = -0.3, size = 2.8, color = "black") +
  facet_wrap(~ method, nrow = 1) +
  scale_color_manual(values = c("NS" = "gray60", "Significant" = "#cc3a36"),
                     name = "raw P < 0.05") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = expression(-log[10](P)),
       y = NULL,
       title = "(B) Per-method cell-type enrichment") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

combined <- pA / pB + plot_layout(heights = c(1.1, 1))
ggsave(file.path(FIG, "EUR_main_figure.pdf"), combined,
       width = 10.5, height = 8)
ggsave(file.path(FIG, "EUR_main_figure.png"), combined,
       width = 10.5, height = 8, dpi = 200, bg = "white")
cat("\nWrote:\n",
    file.path(FIG, "EUR_main_figure.pdf"), "\n",
    file.path(FIG, "EUR_main_figure.png"), "\n", sep = "")
