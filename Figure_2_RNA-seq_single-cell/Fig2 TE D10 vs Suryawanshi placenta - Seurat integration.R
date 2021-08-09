
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# written by Joseph Ernestt, modified by Jaroslav Slamecka below
# scRNA-seq integration
# integrated with dataset from: Suryawanshi et al., 2018: A single-cell survey of the human first-trimester placenta and decidua

# Use integration approach with TE vs. placenta data ----
load("TE_vs_placenta_09102019/placenta.reorder.Rda")
load("TE_vs_placenta_09102019/IS021.seurat.RData")
plSeurat = subset(x = placenta.reorder, subset = orig.ident %in% c("P17","P23"))
plSeurat$cellLabels = Idents(plSeurat)
plSeuratProcessed = SeuratProcess(plSeurat)
Idents(plSeuratProcessed) = "cellLabels"
saveRDS(plSeuratProcessed, file = "TE_vs_placenta_09102019/Placenta_Seurat_Processed_11012019.rds")

teSeurat = IS021.seurat
teSeurat$cellLabels = teSeurat@meta.data$orig.ident
teSeuratProcessed = SeuratProcess(teSeurat)

seuratList = list(te = teSeuratProcessed, pl = plSeuratProcessed)
anchors = FindIntegrationAnchors(seuratList, dims = 1:30)
seuratIntegrated = IntegrateData(anchors, dims = 1:30)
DefaultAssay(seuratIntegrated) = "integrated"

seuratIntegratedProcessed = ScaleData(seuratIntegrated)
seuratIntegratedProcessed = RunPCA(seuratIntegratedProcessed)
seuratIntegratedProcessed = RunUMAP(seuratIntegratedProcessed, dims = 1:30)
Idents(seuratIntegratedProcessed) = "cellLabels"
saveRDS(seuratIntegratedProcessed, file = "TE_vs_placenta_09102019/TE-Placenta_Seurat_Processed_Integrated_10312019.rds")

DimPlot(seuratIntegratedProcessed, reduction = "umap", label = T)

png("TE_vs_placenta_09102019/umap_TE-PL-integrated_10312019.png", width = 960, height = 960)
DimPlot(seuratIntegratedProcessed, reduction = "umap", label = T, label.size = 8, pt.size = 1.5) +  
  labs(title = "Placenta + TE, Seurat integrated") + 
  theme(text = element_text(size = 24),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 24))
dev.off()

# Trying just processing them separately and merging (actually processing them separately and running through
# SeuratProcess but skipping normalizing)
mergedSeurat = merge(teSeuratProcessed, plSeuratProcessed)
mergedSeuratProcessed = SeuratProcess(mergedSeurat, prenormalized = T)
Idents(mergedSeuratProcessed) = "cellLabels"
DimPlot(mergedSeuratProcessed, reduction = "umap", label = T)
png("TE_vs_placenta_09102019/umap_TE-PL_normalized-separately_10312019.png", width = 960, height = 960)
DimPlot(mergedSeuratProcessed, reduction = "umap", label = T, label.size = 8, pt.size = 1.5) +  
  labs(title = "Placenta + TE, normalized separately") + 
  theme(text = element_text(size = 24),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 24))
dev.off()
# Gives same result as normalizing together after merging
mergedSeuratProcessedNormTogether = SeuratProcess(mergedSeurat)
Idents(mergedSeuratProcessedNormTogether) = "cellLabels"
DimPlot(mergedSeuratProcessedNormTogether, reduction = "umap", label = T)
# still same


# Jaroslav Slamecka - additional work ----

library(Seurat)
library(ggplot2)
library(extrafont)
library(Cairo)
library(viridis)
library(RColorBrewer)
library(ggsci)

seuratIntegratedProcessed = readRDS(file="TE_PL_Int_AllAnchs_CCA_11182019.rds")

Idents(seuratIntegratedProcessed)

levels(seuratIntegratedProcessed) = c("H9TE_D0", "H9TE_D3_i527", "H9TE_D6_i524", "H9TE_D9_i521",
                                      "VCT", "SCT", "EVT", "VEC", "HC", "EB", "FB1", "FB2", "FB3")

seuratIntegratedProcessed = RenameIdents(seuratIntegratedProcessed,
                                         "H9TE_D0"="D0",
                                         "H9TE_D3_i527"="D3",
                                         "H9TE_D6_i524"="D6",
                                         "H9TE_D9_i521"="D9")


# viridis colors
plot.colors = c(viridis_pal(begin=0.2,end=0.7, option="D")(4),
                viridis_pal(begin=0.2,end=0.8, direction=-1, option="C")(9))

# update_geom_defaults - makes it possible to manipulate the labels
# because otherwise they are not affected by theme changes
DimPlot(seuratIntegratedProcessed, reduction="umap", label=T, label.size=10, pt.size=1.2) +
  # scale_color_viridis(discrete = TRUE) +
  # scale_color_manual(values=c("#481C6EFF", "#3A548CFF", "#297A8EFF", "#9BD93CFF", "cyan4", "tomato2", "lightpink3", "chartreuse4")) +
  scale_color_manual(values=plot.colors) +
  # scale_color_discrete(c=100, l=60) +
  update_geom_defaults("text", list(family="Noto Sans Cond", size=36)) + 
  theme(text=element_text(size = 36, family="Noto Sans Cond"),
        axis.text=element_text(size=24),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.text=element_text(size=28, angle=50),
        # legend.key=element_rect(size=2, color="black"), # allows to see what happens when adjusting legend parameters
        legend.key.width=unit(0.1,"in"),
        legend.spacing.x=unit(0.01,"in"), legend.spacing.y=unit(0.2,"in"),
        legend.direction="horizontal", legend.position="top", legend.justification="left") +
  # theme_update(legend.key=element_rect(size=2, color="black"), legend.key.size=unit(0.2,"in")) +
  guides(color=guide_legend(label.position="top", nrow=1, override.aes=list(size=5)))
ggsave(filename="TE vs Suryawanshi placenta sc survey - Seurat integrated.pdf",
       device=cairo_pdf,
       width=8.2, height=8.6)


