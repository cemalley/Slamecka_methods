## 08-25-2020
## Evaluating clustering at multiple resolutions for IS021 data

setwd("Z:/NGS_related/Chromium/IS021/")
require(Seurat)
require(clustree)
load("2019-6-TE.sc.RData")

# UMAP at default resolution
DimPlot(tenx.seurat, reduction = "umap", label = T)

# Run clustering at range of resolutions
seuratData = tenx.seurat
seuratData = FindClusters(seuratData, 
                          resolution = c(seq(0.01,0.2,0.01),
                                         seq(0.3,0.8,0.1)))
# clustree plot at resolution 0.01-0.08
g = clustree(seuratData@meta.data[, colnames(seuratData@meta.data) %in% paste0("RNA_snn_res.", seq(0.01,0.2,0.01))], 
             prefix = "RNA_snn_res.")
ggsave("IS021_Clustree_0.01-0.2_08262020.pdf", plot = g, width = 7, height = 12)

# clustree plot at resolution 0.2-0.8
g = clustree(seuratData@meta.data[, colnames(seuratData@meta.data) %in% paste0("RNA_snn_res.", seq(0.2,0.8,0.1))], 
             prefix = "RNA_snn_res.")
ggsave("IS021_Clustree_0.2-0.8_08262020.pdf", plot = g, width = 7, height = 10)

# UMAP at all resolutions
for(i in c(seq(0.01,0.2,0.01), seq(0.3,0.8,0.1))) {
  resColumn = paste0("RNA_snn_res.", i)
  g = DimPlot(seuratData, reduction = "umap", group.by = resColumn,
              label = T) + 
    labs(title = paste0("Resolution: ", i)) + 
    theme(
      # text = element_text(size = 24),
      plot.title = element_text(hjust = 0.5)
      # axis.text = element_text(size = 24)
    )
  fileName = paste0("UMAP_Multiple_Resolutions_08262020/IS021_UMAP_Res_", i, ".pdf")
  ggsave(fileName, plot = g, width = 7, height = 7)
}

# FindAllMarkers at range of resolutions
resVec = paste0("RNA_snn_res.", seq(0.01,0.07,0.01))
markerList = lapply(resVec, function(resColumn) {
  Idents(seuratData) = resColumn
  FindAllMarkers(seuratData)
})
names(markerList) = resVec

# Keep unique positive DE genes for each cluster with logFC > 0.25 and adjusted p-value <= 0.001
markerListUnique = markerList
markerListUnique = lapply(markerListUnique, function(res) {
  res[res$avg_logFC > 0.25 & res$p_val_adj <= 0.001, ]
})
markerListUnique = lapply(markerListUnique, function(res) {
  res[!(res$gene %in% unique(res$gene[duplicated(res$gene)])), ]
})

# Barplots showing # DE genes per cluster and total DE genes for range of resolutions
pdf("IS021_Unique-Positive-DE-Per-Cluster_0.01-0.07_08262020.pdf", width = 7, height = 10)
par(mfrow = c(4, 2))
for(res in names(markerListUnique)) {
  resVal = sub("RNA_snn_res.(.+)$", replacement = "\\1", x = res)
  uniqueDEByCluster = table(markerListUnique[[res]]$cluster)
  xx = barplot(uniqueDEByCluster, 
               main = paste0("Resolution: ", resVal, "\nTotal: ", sum(uniqueDEByCluster)), 
               xlab = "Cluster", ylab = "# DE genes",
               ylim = c(0, max(uniqueDEByCluster)*1.2))
  text(x = xx, y = uniqueDEByCluster, label = uniqueDEByCluster, pos = 3, cex = 1)
}
dev.off()

# Fig 4B: UMAP with resolution = 0.07
g = DimPlot(seuratData, reduction = "umap", group.by = "RNA_snn_res.0.07",
            label = T) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = "top", legend.direction = "horizontal", legend.justification = "center",
        legend.text = element_text(angle = 45)) +
  guides(color = guide_legend(label.position = "top", override.aes = list(size = 2), 
                              nrow = 1, byrow = T))
ggsave("Fig4B_IS021_UMAP_Res0.07.pdf", plot = g, width = 4, height = 4)

# Heatmap with top 10 positive DE genes per cluster at resolution = 0.07
genes = c()
for(clust in unique(markerListUnique$RNA_snn_res.0.07$cluster)) {
  print(clust)
  df = markerListUnique$RNA_snn_res.0.07[markerListUnique$RNA_snn_res.0.07$cluster == clust, ]
  df = df[order(df$avg_logFC, decreasing = T), ]
  genes = c(genes, df$gene[1:10])
}
genes = unique(genes)
seuratData_scaleAll = ScaleData(seuratData, features = rownames(seuratData))
seuratData_scaleAll$RNA_snn_res.0.07 = factor(seuratData_scaleAll$RNA_snn_res.0.07,
                                              levels = c(0,5,2,4,1,3,6))
ramp = colorRamp(c("skyblue4", "beige", "firebrick4"))
g = DoHeatmap(seuratData_scaleAll, features = genes, group.by = "RNA_snn_res.0.07", raster = F) +
  NoLegend() +
  scale_y_discrete(position = "right") +
  scale_fill_gradientn(colors = rgb(ramp(seq(0,1,length=256)),max=256))
  # scale_fill_gradientn(colors = colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100))

ggsave("Fig4G_IS021_Heatmap_ClusterMarkers_08262020.pdf", plot = g, width = 4, height = 8)

# Module scores
te_tb_markers = c("GATA3","TFAP2A","CDX2","IGFBP3","TEAD3","TEAD4","PAPPA2","KRT7",
                  "KRT18","MSX2","TAC3","XAGE2","VGLL1","PERP","GABRP","TPM1")
all(te_tb_markers %in% rownames(seuratData))

sct_markers = c("CGA","CGB3","DAB2","PGF","GDF15","INSL4","PHLDA2","S100P","TFPI","LGALS16",
                "ERVW-1","ERVFRD-1","CYP11A1","TCL1B","HOPX","KRT23")
all(sct_markers %in% rownames(seuratData))

markerList = list(te_tb = te_tb_markers, sct = sct_markers)
seuratData_withModuleScore = AddModuleScore(seuratData, features = markerList)
colnames(seuratData_withModuleScore@meta.data)[colnames(seuratData_withModuleScore@meta.data) == "Cluster1"] = "te_tb"
colnames(seuratData_withModuleScore@meta.data)[colnames(seuratData_withModuleScore@meta.data) == "Cluster2"] = "sct"
Idents(seuratData_withModuleScore) = "RNA_snn_res.0.07"
FeaturePlot(seuratData_withModuleScore, features = "te_tb", label = T)

# Fig 4 E: Module scores for trophectoderm/trophoblast markers
g = FeaturePlot(seuratData_withModuleScore, features = "te_tb", reduction = "umap", label = T) +
  ggtitle("Trophectoderm/trophoblast markers") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) 
ggsave("Fig4E_IS021_ModuleScore_UMAP_TE_TB_Res0.07.pdf", plot = g, width = 5, height = 4)

# Fig 4 F: Module scores for syncytiotrophoblast markers
g = FeaturePlot(seuratData_withModuleScore, features = "sct", reduction = "umap", label = T) +
  ggtitle("Syncytiotrophoblast markers") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) 
ggsave("Fig4F_IS021_ModuleScore_UMAP_SCT_Res0.07.pdf", plot = g, width = 5, height = 4)

# Violin plots for module scores
g = VlnPlot(seuratData_withModuleScore, features = "te_tb", 
            group.by = "RNA_snn_res.0.07", pt.size = 0.1) +
  ggtitle("Trophectoderm/trophoblast markers") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  NoLegend() +
  xlab("Cluster") +
  ylab("Module score")
ggsave("Fig4_IS021_ModuleScore_VlnPlot_TE_TB_Points_Res0.07_08262020.pdf", plot = g, width = 9, height = 6)

g = VlnPlot(seuratData_withModuleScore, features = "te_tb", 
            group.by = "RNA_snn_res.0.07", pt.size = 0) +
  ggtitle("Trophectoderm/trophoblast markers") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  NoLegend() +
  xlab("Cluster") +
  ylab("Module score")
ggsave("Fig4_IS021_ModuleScore_VlnPlot_TE_TB_NoPoints_Res0.07_08262020.pdf", plot = g, width = 9, height = 6)


g = VlnPlot(seuratData_withModuleScore, features = "sct", 
            group.by = "RNA_snn_res.0.07", pt.size = 0.1) +
  ggtitle("Syncytiotrophoblast markers") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), 
        plot.title = element_text(hjust = 0.5)) +
  NoLegend() +
  xlab("Cluster") +
  ylab("Module score")
ggsave("Fig4_IS021_ModuleScore_VlnPlot_SCT_Points_Res0.07_08262020.pdf", plot = g, width = 9, height = 6)

g = VlnPlot(seuratData_withModuleScore, features = "sct", 
            group.by = "RNA_snn_res.0.07", pt.size = 0) +
  ggtitle("Syncytiotrophoblast markers") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  NoLegend() +
  xlab("Cluster") +
  ylab("Module score")
ggsave("Fig4_IS021_ModuleScore_VlnPlot_SCT_NoPoints_Res0.07_08262020.pdf", plot = g, width = 9, height = 6)
