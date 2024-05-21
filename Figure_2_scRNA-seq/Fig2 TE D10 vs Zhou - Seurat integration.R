
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# written by Joseph Ernestt, modified by Jaroslav Slamecka below
# scRNA-seq integration
# integrated with dataset from: Zhou et al., 2019: Reconstituting the transcriptome and DNA methylome landscapes of human implantation



# Load packages
require(Seurat)

# Source custom R functions
source("../rscripts_09102019/r_functions_09092019.R")



# Process Zhou data ----

# Zhou metadata (table already created from supplementary info)
metaData = read.table("metaData_10302019.txt", header = T, as.is = T)
rownames(metaData) = metaData$Sample

# Download Zhou TPM data 
# At command line:
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109555/suppl/GSE109555_All_Embryo_TPM.txt.gz
# gunzip GSE109555_All_Embryo_TPM.txt.gz

# Read in TPM data
embCounts = as.matrix(read.table("GSE109555_All_Embryo_TPM.txt"))

# Create Seurat object
zhouSeurat = CreateSeuratObject(embCounts, meta.data = metaData)

# Process Zhou using standard Seurat workflow. Use isTPM = T to skip Normalize step and use log(TMP + 1) instead
zhouSeuratProcessed = SeuratProcess(zhouSeurat, isTPM = T)

# Create cellLabels metadata column (same as Lineage) to facilitate combining with TE data
zhouSeuratProcessed$cellLabels = zhouSeuratProcessed@meta.data$Lineage
Idents(zhouSeuratProcessed) = "cellLabels"

# Save Seurat object if necessary
saveRDS(zhouSeuratProcessed, file = "Zhou_Seurat_Processed_11012019.rds")


# Seurat integration procedure ----

# Read in TE data (already processed and with cellLabels column) and set Idents to cellLabels
teSeuratProcessed = readRDS("../TE_Seurat_Processed_11012019.rds")
Idents(teSeuratProcessed) = "cellLabels"

# Read in Zhou Seurat data if previous steps weren't performed
# zhouSeuratProcessed = readRDS("Zhou_Seurat_Processed_11012019.rds")
# Idents(zhouSeuratProcessed) = "cellLabels"

# Combine Seurat objects into R list
seuratList = list(te = teSeuratProcessed, zhou = zhouSeuratProcessed)

# Find anchors
anchors = FindIntegrationAnchors(seuratList, dims = 1:30)

# Integrate data. This creates a new element called "integrated" in Assays (use this instead of RNA)
seuratIntegrated = IntegrateData(anchors, dims = 1:30)

DefaultAssay(seuratIntegrated) = "integrated"

# Run regular Seurat processing on integrated data
seuratIntegratedProcessed = ScaleData(seuratIntegrated)
seuratIntegratedProcessed = RunPCA(seuratIntegratedProcessed)
seuratIntegratedProcessed = RunUMAP(seuratIntegratedProcessed, dims = 1:30)

# Create UMAP plot of integrated and processed data
png("umap_zhou-TE-integrated_lineage_10312019.png", width = 960, height = 960)
DimPlot(seuratIntegratedProcessed, reduction = "umap", label = T, label.size = 8, pt.size = 1.5) +  
  labs(title = "Zhou 2019 + TE, Seurat integrated") + 
  theme(text = element_text(size = 24),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 24))
dev.off()

# Save integrated and processed Seurat object
saveRDS(seuratIntegratedProcessed, file = "TE_Zhou_Integrated_Processed_11042019.rds")



# Jaroslav Slamecka - additional work ----

library(Seurat)
library(ggplot2)
library(extrafont)
library(Cairo)
library(viridis)
library(RColorBrewer)
library(ggsci)

seuratIntegratedProcessed = readRDS(file="TE_Zhou_Integrated_AllAncFeats_CCA_Abbrev_11182019.rds")

Idents(seuratIntegratedProcessed)

# viridis colors
plot.colors = c(viridis_pal(begin=0.2,end=0.7, option="D")(4),
                viridis_pal(begin=0.3,end=0.7, direction=-1, option="C")(4))

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
        legend.text=element_text(size=28, angle=45),
        # legend.key=element_rect(size=2, color="black"), # allows to see what happens when adjusting legend parameters
        legend.key.width=unit(0.1,"in"),
        legend.spacing.x=unit(0.01,"in"), legend.spacing.y=unit(0.2,"in"),
        legend.direction="horizontal", legend.position="top", legend.justification="center") + # legend.position=c(0.42,0.96) for an inset legend
  # theme_update(legend.key=element_rect(size=2, color="black"), legend.key.size=unit(0.2,"in")) +
  guides(color=guide_legend(label.position="top", nrow=1, override.aes=list(size=5)))
ggsave(filename="TE vs Zhou Nature 2019 - Seurat integrated.pdf",
       device=cairo_pdf,
       width=8.2, height=8.6)


