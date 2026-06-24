# =================================================================
# AddModuleScore boxplot: ONLY excitatory + inhibitory neuron subtypes
# Excludes Purkinje, glia, vascular (those are major-level categories)
# =================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(ggpubr); library(data.table)
})

OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

# Read per-cell scores (gzipped csv ~300 MB)
md <- fread(file.path(DAT, "PGC_MDD_AddModuleScore_per_cell.csv.gz"))
cat("Cells:", nrow(md), "\n")
cat("Subtypes:", paste(sort(unique(md$final_anno)), collapse = ", "), "\n")

# Excitatory + Inhibitory only
NEURONAL <- c(
  "Ex-L2/3", "Ex-L2/4", "Ex-L4/6", "Ex-L5", "Ex-L5/6", "Ex-L6", "Ex-NRGN", "Ex_mix",
  "In_CALM1", "In_LAMP5", "In_PVALB", "In_SHANK2", "In_SST", "In_VIP"
)
nd <- md[final_anno %in% NEURONAL]
cat("Neuronal cells:", nrow(nd), "\n")
cat("Subtypes kept:", paste(sort(unique(nd$final_anno)), collapse = ", "), "\n")

# Custom ordering: Ex first (layered), then In (by marker), Ex_mix / Ex-NRGN end of Ex group
ORDER <- c("Ex-L2/3", "Ex-L2/4", "Ex-L4/6", "Ex-L5", "Ex-L5/6", "Ex-L6", "Ex-NRGN", "Ex_mix",
           "In_PVALB", "In_SST", "In_VIP", "In_LAMP5", "In_CALM1", "In_SHANK2")
nd[, final_anno := factor(final_anno, levels = ORDER)]

# Color by cell class
nd[, class := ifelse(grepl("^Ex", final_anno), "Excitatory", "Inhibitory")]

# ----- Boxplot v1: ordered by canonical (layer / marker) -----
p_canon <- ggplot(nd, aes(x = final_anno, y = PGC_MDD_score, fill = class)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  scale_fill_manual(values = c("Excitatory" = "#4ea64a", "Inhibitory" = "#8e4c99")) +
  labs(x = NULL,
       y = "PGC-MDD-295 module score",
       title = "PGC 295 MDD genes — expression by neuronal subtype",
       subtitle = sprintf("Excitatory (n=%d) + Inhibitory (n=%d) only | %d cells | n_genes=278",
                          sum(nd$class == "Excitatory"),
                          sum(nd$class == "Inhibitory"),
                          nrow(nd))) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        legend.position = "top")

ggsave(file.path(FIG, "PGC_MDD_AddModuleScore_neurons_canonical.pdf"), p_canon,
       width = 9, height = 5)
ggsave(file.path(FIG, "PGC_MDD_AddModuleScore_neurons_canonical.png"), p_canon,
       width = 9, height = 5, dpi = 200, bg = "white")

# ----- Boxplot v2: ordered by median score (top-expressing first) -----
p_sorted <- ggplot(nd, aes(x = reorder(final_anno, PGC_MDD_score, FUN = median),
                            y = PGC_MDD_score, fill = class)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  scale_fill_manual(values = c("Excitatory" = "#4ea64a", "Inhibitory" = "#8e4c99")) +
  labs(x = NULL,
       y = "PGC-MDD-295 module score",
       title = "PGC 295 MDD genes — by neuronal subtype (sorted by median score)",
       subtitle = sprintf("Excitatory (n=%d) + Inhibitory (n=%d) | %d cells | n_genes=278",
                          sum(nd$class == "Excitatory"),
                          sum(nd$class == "Inhibitory"),
                          nrow(nd))) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        legend.position = "top")

ggsave(file.path(FIG, "PGC_MDD_AddModuleScore_neurons_sorted.pdf"), p_sorted,
       width = 9, height = 5)
ggsave(file.path(FIG, "PGC_MDD_AddModuleScore_neurons_sorted.png"), p_sorted,
       width = 9, height = 5, dpi = 200, bg = "white")

# Summary table
summ <- nd[, .(n = .N,
               mean_score = round(mean(PGC_MDD_score), 4),
               median_score = round(median(PGC_MDD_score), 4),
               sd_score = round(sd(PGC_MDD_score), 4)),
           by = final_anno] %>% .[order(-mean_score)]
print(summ)
write.csv(summ, file.path(DAT, "PGC_MDD_AddModuleScore_neurons_summary.csv"),
          row.names = FALSE)

cat("\nWrote:\n",
    file.path(FIG, "PGC_MDD_AddModuleScore_neurons_canonical.pdf"), "\n",
    file.path(FIG, "PGC_MDD_AddModuleScore_neurons_sorted.pdf"), "\n", sep = "")
