# =================================================================
# scDRS PGC295 dot charts:
#   (1) MAJOR (8 cell types) × 3 GWAS facets
#   (2) SUBTYPE (14 Ex+In) × 3 GWAS facets
# Uses MAGMA-z-weighted PGC 295 gene sets per GWAS (HW/EUR/DIV)
# =================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(ggpubr)
})

OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

norm_ct <- function(x) {
  y <- gsub("[ /]", ".", x)
  y <- gsub("_", ".", y)
  y <- gsub("^Ex\\.([LN])", "Ex-\\1", y)
  y <- gsub("^Ex\\.mix$", "Ex-mix", y)
  y <- gsub("^In\\.([CLPSV])", "In-\\1", y)
  y
}

read_scdrs <- function(path, gwas) {
  d <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data.frame(GWAS = gwas,
             cell_type = norm_ct(d$group),
             P = d$assoc_mcp,
             z = d$assoc_mcz,
             n_cell = d$n_cell)
}

# ----- MAJOR (8) -----
CT_MAJOR <- c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons",
              "Endothelial.cells","Oligodendrocytes","OPCs",
              "Microglia","Astrocytes")
PAL_MAJOR <- c("Excitatory.neurons" = "#4ea64a",
               "Inhibitory.neurons" = "#8e4c99",
               "Purkinje.neurons"   = "#a05528",
               "Endothelial.cells"  = "#3777ac",
               "Oligodendrocytes"   = "#e47faf",
               "OPCs"               = "#b698c5",
               "Microglia"          = "#e88f18",
               "Astrocytes"         = "#d5231d")

major <- bind_rows(
  read_scdrs(file.path(OUT, "downstream/PGC295_HW.scdrs_group.anno"),  "HW"),
  read_scdrs(file.path(OUT, "downstream/PGC295_EUR.scdrs_group.anno"), "EUR"),
  read_scdrs(file.path(OUT, "downstream/PGC295_DIV.scdrs_group.anno"), "DIV")
) %>%
  mutate(neglogP = -log10(P),
         GWAS = factor(GWAS, levels = c("HW","EUR","DIV")),
         cell_type = factor(cell_type, levels = rev(CT_MAJOR)))

write.csv(major, file.path(DAT, "scdrs_PGC295_major.csv"), row.names = FALSE)

p_major <- ggdotchart(
  major, x = "cell_type", y = "neglogP",
  color = "cell_type", palette = PAL_MAJOR[CT_MAJOR],
  sorting = "descending",
  add = "segments", add.params = list(color = "lightgray", size = 2),
  group = "cell_type", dot.size = 6,
  label = round(major$neglogP, 2),
  font.label = list(color = "black", size = 7, vjust = 0.5),
  facet.by = "GWAS",
  ggtheme = theme_pubr()) +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "lightgray") +
  theme(axis.text.x = element_text(size = 9, color = 'black'),
        axis.text.y = element_text(size = 9, color = 'black'),
        legend.position = 'none',
        strip.text = element_text(face = "bold", size = 11)) +
  labs(y = expression(-log[10](assoc_mcp)),
       title = "scDRS — PGC 295 high-confidence MDD genes (MAGMA z-weighted)",
       subtitle = "Major cell types | n_ctrl = 500 (min P = 1/501 ≈ 0.002)") +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1))

ggsave(file.path(FIG, "scDRS_PGC295_major_dotchart.pdf"), p_major,
       width = 12, height = 4)
ggsave(file.path(FIG, "scDRS_PGC295_major_dotchart.png"), p_major,
       width = 12, height = 4, dpi = 200, bg = "white")

# ----- SUBTYPE (14 Ex+In only) -----
CT_NEU <- c("Ex-L2.3","Ex-L2.4","Ex-L4.6","Ex-L5","Ex-L5.6","Ex-L6","Ex-NRGN","Ex-mix",
            "In-PVALB","In-SST","In-VIP","In-LAMP5","In-CALM1","In-SHANK2")
PAL_NEU <- c(
  rep("#4ea64a", 8),    # Ex green
  rep("#8e4c99", 6)     # In purple
)
names(PAL_NEU) <- CT_NEU

sub <- bind_rows(
  read_scdrs(file.path(OUT, "downstream/PGC295_HW.scdrs_group.final_anno"),  "HW"),
  read_scdrs(file.path(OUT, "downstream/PGC295_EUR.scdrs_group.final_anno"), "EUR"),
  read_scdrs(file.path(OUT, "downstream/PGC295_DIV.scdrs_group.final_anno"), "DIV")
) %>%
  filter(cell_type %in% CT_NEU) %>%
  mutate(neglogP = -log10(P),
         GWAS = factor(GWAS, levels = c("HW","EUR","DIV")),
         cell_type = factor(cell_type, levels = rev(CT_NEU)))

write.csv(sub, file.path(DAT, "scdrs_PGC295_subtype_neurons.csv"), row.names = FALSE)

p_sub <- ggdotchart(
  sub, x = "cell_type", y = "neglogP",
  color = "cell_type", palette = PAL_NEU,
  sorting = "descending",
  add = "segments", add.params = list(color = "lightgray", size = 2),
  group = "cell_type", dot.size = 5,
  label = round(sub$neglogP, 2),
  font.label = list(color = "black", size = 6, vjust = 0.5),
  facet.by = "GWAS",
  ggtheme = theme_pubr()) +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "lightgray") +
  theme(axis.text.x = element_text(size = 9, color = 'black'),
        axis.text.y = element_text(size = 9, color = 'black'),
        legend.position = 'none',
        strip.text = element_text(face = "bold", size = 11)) +
  labs(y = expression(-log[10](assoc_mcp)),
       title = "scDRS — PGC 295 high-confidence MDD genes by neuronal subtype",
       subtitle = "Excitatory (green) + Inhibitory (purple) subtypes | n_ctrl = 500") +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1))

ggsave(file.path(FIG, "scDRS_PGC295_subtype_dotchart.pdf"), p_sub,
       width = 12, height = 6)
ggsave(file.path(FIG, "scDRS_PGC295_subtype_dotchart.png"), p_sub,
       width = 12, height = 6, dpi = 200, bg = "white")

# Print summary
options(width = 200)
cat("\n========== scDRS PGC295 MAJOR ==========\n")
print(major %>% select(GWAS, cell_type, P, z, neglogP) %>%
      arrange(GWAS, desc(neglogP)) %>%
      mutate(P = signif(P, 3), z = round(z, 2), neglogP = round(neglogP, 2)),
      row.names = FALSE)
cat("\n========== scDRS PGC295 SUBTYPE (14 Ex+In) ==========\n")
print(sub %>% select(GWAS, cell_type, P, z, neglogP) %>%
      arrange(GWAS, desc(neglogP)) %>%
      mutate(P = signif(P, 3), z = round(z, 2), neglogP = round(neglogP, 2)),
      row.names = FALSE)

cat("\nWrote:\n",
    file.path(FIG, "scDRS_PGC295_major_dotchart.pdf"), "\n",
    file.path(FIG, "scDRS_PGC295_subtype_dotchart.pdf"), "\n", sep = "")
