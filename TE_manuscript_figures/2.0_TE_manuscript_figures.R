#' Make updated figures for TE manuscript as requested
#' by CM & JS.

# req'd pkgs
x <- c("Seurat", "dplyr", "data.table", "ggplot2", "extrafont", 
       "extrafontdb", "tidyr", "umap", "Cairo")
sapply(x, library, character.only = TRUE)

# load upd seurat object
seur <- readRDS("./upd_tenx.seurat.RDS")

# add annotation to plot for day of differentiation
add <- c(unique(seur@meta.data$day), "diff")
add <- c(paste(add[1], "(hPSC)", sep = " "), add[2], paste(add[3], add[4], sep = "/"), add[5])

loadfonts()

# update figure 4b: Re-label clusters so they
# are 0, 1, 2, 3, 4, 5 & print in sans serif
p <- Seurat::DimPlot(seur, reduction = "umap", label = TRUE, label.size = 10, pt.size = 1.2) + 
  theme(legend.position = "top", legend.justification = "center", 
        text = element_text(size = 36, family = "NotoSans-Condensed"),
        axis.text = element_text(size = 24),
        legend.text = element_text(size=28),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.key.width=unit(0.3,"in"),
        legend.spacing.x=unit(0.01,"in"), 
        legend.spacing.y=unit(0.2,"in")
        ) +
  guides(color = guide_legend(label.position = "top", override.aes = list(size = 5), nrow = 1, byrow = T, 
                              label.theme = element_text(family = "NotoSans-Condensed", size = 24)))
p <- p + annotate("text", x = c(8, 0, 2, 0), y = c(-5.5, -9, 10, 5), label = add, size = 14, family = "NotoSans-Condensed")
ggsave("Figure_4B_updated.png", plot = p, width=8.2, height=8.6)
ggsave("Figure_4B_updated.pdf", plot = p, device=cairo_pdf, width=8.2, height=8.6)


# update figure 4g: filter to remove pseudogenes,
# MT, AC genes; use sky-blue:biege:firebrick red
# schema; add legend; use top 5-6 top DGE genes
# by cluster
gen_to_filt <- fread("./35387_genes_filter.txt", header = F)

# keep top 50 DGE genes
de <- fread("./dge_clus_1_3_merge.csv")
adj <- de[abs(de$avg_logFC) > 0.25 & de$p_val_adj <= 0.001, ]
adj <- split(adj, adj$cluster)
adj <- unlist(lapply(adj, function(x) {
  y <- x %>%
    arrange(-abs(avg_logFC)) %>%
    filter(!gene %in% gen_to_filt[,1]) %>%
    filter(!grepl("^AC|^MT-", gene)) %>%
    pull(gene)
  y <- unique(y)
  y <- y[c(1:10)]
  }))
adj <- as.vector(adj)
adj <- unique(adj)
# re-adjust ordering of genes for better viewing
# in heatmap based on high exp in clusters
adj <- adj[c(10, 12, 13, 14, 16, 18, 8, 15, 3, 19,
             6, 9, 7, 30, 11, 2, 4, 29, 27, 28, 1, 
             5, 22, 17, 20, 24, 23, 25, 21, 26, 2)]

seur$day <- as.factor(seur$day)
# set identity to updated cluster names
Idents(seur) <- 'upd_res'

# create heatmap for specified top 30 DGE genes 
# by cluster w/specified color palette, legend
CairoPDF(file = "Figure_4G_updated.pdf", width = 32, height = 16)
DoHeatmap(seur, features = adj, raster = F, size = 10) +
  scale_fill_gradientn("Scaled\nExpression", 
                       colors = colorRampPalette(c('sky blue', 'beige', 'firebrick4'))(100)) +  
  guides(color = FALSE) +
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size = 26, color = "black"),
        text = element_text(size = 26, family = "NotoSans-Condensed"),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 26),
        legend.key.width = unit(0.1,"in"),
        legend.key.height = unit(0.4, "in")
        )
dev.off()

# DoMultiBarHeatmap(seur, features = adj,
#                   group.by = 'upd_res',
#                   additional.group.by = 'day') +
#   scale_fill_gradientn(colors = colorRampPalette(c('sky blue', 'beige', 'firebrick4'))(100)) +  
#   guides(color=FALSE)

# 4e/f avg exp or module score umaps for list of
# trophectoderm or syncytiotrophoblast markers
te_list <- c("GATA3", "TFAP2A", "CDX2", "IGFBP3",
             "TEAD3", "TEAD4", "PAPPA2", "KRT7",
             "KRT18", "MSX2", "TAC3", "XAGE2", 
             "VGLL1", "PERP", "GABRP", "TPM1")
syn_list <- c("CGA", "CGB3", "DAB2", "PGF", "GDF15",
              "INSL4", "PHLDA2", "S100P", "TFPI",
              "LGALS16", "ERVW-1", "ERVFRD-1", 
              "CYP11A1", "TCL1B", "HOPX", "KRT23")

# compute average gene expression for targeted
# goi
te_avg <- colMeans(x = as.matrix(seur@assays$RNA@data[te_list, ]), na.rm = TRUE)
if (all(names(te_avg) == rownames(x = seur@meta.data))) {
seur@meta.data$te_avg <- te_avg
}
syn_avg <- colMeans(x = as.matrix(seur@assays$RNA@data[syn_list, ]), na.rm = TRUE)
if (all(names(syn_avg) == rownames(x = seur@meta.data))) {
  seur@meta.data$syn_avg <- syn_avg
}

# plot avg exp umap for targeted goi
CairoPDF(file = "Figure_4E_avg_exp.pdf", width = 24, height = 24)
FeaturePlot(seur, features = "te_avg", reduction = "umap", pt.size = 1.2, min.cutoff="q9", order = TRUE) + 
  ggtitle("Trophectoderm/trophoblast average\nexpression of selected markers") +
  theme(text = element_text(family = "NotoSans-Condensed"),
        plot.title=element_text(size=32, hjust = 0.5),
        axis.title = element_blank(),
        legend.key.size=unit(0.3,"in"), 
        legend.text=element_text(size = 18),
        legend.title = element_text(size = 18)
        ) + labs(color = 'Average expression') 
dev.off()

CairoPDF(file = "Figure_4F_avg_exp.pdf", width = 24, height = 24)
FeaturePlot(seur, features = "syn_avg", reduction = "umap", pt.size = 1.2, min.cutoff="q9", order = TRUE) + 
  ggtitle("Syncytiotrophoblast average\nexpression of selected markers") +
  theme(text = element_text(family = "NotoSans-Condensed"),
        plot.title=element_text(size=32, hjust = 0.5),
        axis.title = element_blank(),
        legend.key.size=unit(0.3,"in"), 
        legend.text=element_text(size = 18),
        legend.title = element_text(size = 18)
  ) + labs(color = 'Average expression') 
dev.off()

# compute module score for targeted goi
te_list <- list(c("GATA3", "TFAP2A", "CDX2", "IGFBP3",
             "TEAD3", "TEAD4", "PAPPA2", "KRT7",
             "KRT18", "MSX2", "TAC3", "XAGE2", 
             "VGLL1", "PERP", "GABRP", "TPM1"))
syn_list <- list(c("CGA", "CGB3", "DAB2", "PGF", "GDF15",
              "INSL4", "PHLDA2", "S100P", "TFPI",
              "LGALS16", "ERVW-1", "ERVFRD-1", 
              "CYP11A1", "TCL1B", "HOPX", "KRT23"))
seur <- AddModuleScore(seur, features = te_list, name = "te_mod_score")
seur <- AddModuleScore(seur, features = syn_list, name = "syn_mod_score")

CairoPDF(file = "Figure_4E_mod_score.pdf", width = 24, height = 24)
FeaturePlot(seur, features = "te_mod_score1", reduction = "umap", pt.size = 1.2, min.cutoff="q9", order = TRUE) + 
  ggtitle("Trophectoderm/trophoblast module\nscore of selected markers") +
  theme(text = element_text(family = "NotoSans-Condensed"),
        plot.title=element_text(size=32, hjust = 0.5),
        axis.title = element_blank(),
        legend.key.size=unit(0.3,"in"), 
        legend.text=element_text(size = 18),
        legend.title = element_text(size = 18)
  ) + labs(color = 'Average module score') 
dev.off()

CairoPDF(file = "Figure_4F_mod_score.pdf", width = 24, height = 24)
FeaturePlot(seur, features = "syn_mod_score1", reduction = "umap", pt.size = 1.2, min.cutoff="q9", order = TRUE) + 
  ggtitle("Syncytiotrophoblast module\nscore of selected markers") +
  theme(text = element_text(family = "NotoSans-Condensed"),
        plot.title=element_text(size=32, hjust = 0.5),
        axis.title = element_blank(),
        legend.key.size=unit(0.3,"in"), 
        legend.text=element_text(size = 18),
        legend.title = element_text(size = 18)
  ) + labs(color = 'Average module score') 
dev.off()
