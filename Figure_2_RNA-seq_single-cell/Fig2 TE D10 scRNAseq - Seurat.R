
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# load libraries ----

library(dplyr)
library(Seurat)
library(ggplot2)
library(extrafont)
library(viridis)
library(Cairo)



# load 10X datasets ----

# tenx.data = list()
samples = c("H9TE_D0", "H9TE_D3_i527", "H9TE_D6_i524", "H9TE_D9_i521")

# individual samples
tenx.data.D0 = Read10X(data.dir=paste0("../COUNTS/",samples[1],"/outs/filtered_gene_bc_matrices/GRCh38/"))
tenx.data.D3 = Read10X(data.dir=paste0("../COUNTS/",samples[2],"/outs/filtered_gene_bc_matrices/GRCh38/"))
tenx.data.D6 = Read10X(data.dir=paste0("../COUNTS/",samples[3],"/outs/filtered_gene_bc_matrices/GRCh38/"))
tenx.data.D9 = Read10X(data.dir=paste0("../COUNTS/",samples[4],"/outs/filtered_gene_bc_matrices/GRCh38/"))

# initialize the Seurat object with the raw (non-normalized data) ----
# give each sample a unique Project ID
tenx.seurat.D0 = CreateSeuratObject(counts=tenx.data.D0, project=samples[1], min.cells=3, min.features=200)
tenx.seurat.D3 = CreateSeuratObject(counts=tenx.data.D3, project=samples[2], min.cells=3, min.features=200)
tenx.seurat.D6 = CreateSeuratObject(counts=tenx.data.D6, project=samples[3], min.cells=3, min.features=200)
tenx.seurat.D9 = CreateSeuratObject(counts=tenx.data.D9, project=samples[4], min.cells=3, min.features=200)

# clean up
rm(tenx.data.D0, tenx.data.D3, tenx.data.D6, tenx.data.D9)

# merge Seurat objects ----
tenx.seurat = merge(x = tenx.seurat.D0, y = c(tenx.seurat.D3, tenx.seurat.D6, tenx.seurat.D9), project="H9TE", add.cell.ids=samples)

# check number of cells from each original Seurat object
table(tenx.seurat$orig.ident)

# clean up
rm(tenx.seurat.D0, tenx.seurat.D3, tenx.seurat.D6, tenx.seurat.D9)



# PRE-PROCESS ----

tenx.seurat = NormalizeData(tenx.seurat)
tenx.seurat = FindVariableFeatures(tenx.seurat)
tenx.seurat = ScaleData(tenx.seurat)
tenx.seurat = RunPCA(tenx.seurat)
tenx.seurat = RunUMAP(tenx.seurat, dim="pca", dims=1:20)
tenx.seurat = FindNeighbors(tenx.seurat)
tenx.seurat = FindClusters(tenx.seurat, resolution=0.07)
tenx.seurat = RunTSNE(tenx.seurat)



# UMAP plot ----

UMAPPlot(tenx.seurat, pt.size=1.2, label=TRUE, label.size=10) +
  update_geom_defaults("text", list(family="Noto Sans Cond", size=36)) +
  scale_color_discrete(c=100, l=58) +
  theme(text=element_text(size = 36, family="Noto Sans Cond"),
        axis.text=element_text(size=24),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.text=element_text(size=28, angle=45),
        # legend.key=element_rect(size=2, color="black"), # allows to see what happens when adjusting legend parameters
        legend.key.width=unit(0.1,"in"),
        legend.spacing.x=unit(0.01,"in"), legend.spacing.y=unit(0.2,"in"),
        legend.direction="horizontal", legend.position="top", legend.justification="center") + # legend.position=c(0.42,0.96) for an inset legend
  guides(color=guide_legend(label.position="top", nrow=1, override.aes=list(size=5)))
ggsave(filename="umap.clusters.(6).pdf", device=cairo_pdf, width=8.2, height=8.6)

# color by sample

plot.colors = c(viridis_pal(begin=0.2,end=0.7, option="D")(4))
UMAPPlot(tenx.seurat, group.by="orig.ident", pt.size=1.2) +
  scale_color_discrete(c=100, l=58, labels=c("D0","D3","D6","D9")) +
  # scale_color_manual(values=plot.colors, labels=c("D0","D3","D6","D9")) +
  update_geom_defaults("text", list(family="Noto Sans Cond", size=36)) +
  theme(text=element_text(size = 36, family="Noto Sans Cond"),
        axis.text=element_text(size=24),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.text=element_text(size=28, angle=45),
        # legend.key=element_rect(size=2, color="black"), # allows to see what happens when adjusting legend parameters
        legend.key.width=unit(0.1,"in"),
        legend.spacing.x=unit(0.01,"in"), legend.spacing.y=unit(0.2,"in"),
        legend.direction="horizontal", legend.position="top", legend.justification="center") + # legend.position=c(0.42,0.96) for an inset legend
  guides(color=guide_legend(label.position="top", nrow=1, override.aes=list(size=5)))
ggsave(filename="umap.color.by.sample.pdf", device=cairo_pdf, width=8.2, height=8.6) # add "viridis"
rm(plot.colors)



# FEATURE PLOTS ----

library(extrafont)
library(Cairo)
library(patchwork)

gene.group = list(general=c("GATA3","TFAP2A","CDX2","IGFBP3","KRT7",
                            "TEAD4","XAGE2","VGLL1","MSX2","NR2F2",
                            "EPAS1","PERP","PAPPA2","GABRP","TPM1"),
                  synTB=c("CGA","DAB2","PGF","INSL4","PHLDA2",
                          "S100P","LGALS16","TCL1B","ERVW-1","ERVFRD-1"))

# general TE markers
CairoPDF(file=paste0("feature.plot.final-general.pdf"), width=30, height=18) # 24x24 for 4x4 markers
feature.plots = FeaturePlot(tenx.seurat, features=gene.group$general, pt.size=1.2, min.cutoff="q9", order=TRUE, combine=FALSE)
feature.plots = lapply(X=feature.plots,
                       FUN=function(x) {x + theme(text=element_text(family="Noto Sans Cond"),
                                                  plot.title=element_text(size=32),
                                                  axis.title.x=element_blank(), axis.title.y=element_blank(),
                                                  legend.key.size=unit(0.3,"in"), legend.text=element_text(size=18)
                       )}
)
wrap_plots(feature.plots, ncol=5)
dev.off()

# synTB markers
CairoPDF(file=paste0("feature.plot.final-synTB.pdf"), width=30, height=12)
feature.plots = FeaturePlot(tenx.seurat, features=gene.group$synTB, pt.size=1.2, min.cutoff="q9", order=TRUE, combine=FALSE)
feature.plots = lapply(X=feature.plots,
                       FUN=function(x) {x + theme(text=element_text(family="Noto Sans Cond"),
                                                  plot.title=element_text(size=32),
                                                  axis.title.x=element_blank(), axis.title.y=element_blank(),
                                                  legend.key.size=unit(0.3,"in"), legend.text=element_text(size=18)
                       )}
)
wrap_plots(feature.plots, ncol=5)
dev.off()

rm(gene.group, feature.plots)

# pluripotency markers
CairoPDF(file=paste0("feature.plot.pluripotency.markers.pdf"), width=30, height=6) # 24x24 for 4x4 markers
feature.plots = FeaturePlot(tenx.seurat, features=c("POU5F1","SOX2","NANOG","TDGF1","UTF1"), pt.size=1.2, min.cutoff="q9", order=TRUE, combine=FALSE)
feature.plots = lapply(X=feature.plots,
                       FUN=function(x) {x + theme(text=element_text(family="Noto Sans Cond"),
                                                  plot.title=element_text(size=32),
                                                  axis.title.x=element_blank(), axis.title.y=element_blank(),
                                                  legend.key.size=unit(0.3,"in"), legend.text=element_text(size=18)
                       )}
)
wrap_plots(feature.plots, ncol=5)
dev.off()



# MARKER analysis ----

# find markers for every cluster compared to all remaining cells, report only the positive ones
te.markers = FindAllMarkers(tenx.seurat, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# display top 2 for each cluster
te.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
# write into a file
library(data.table)
fwrite(all.markers, file="all.markers.txt", sep="\t")



# HEATMAP ----

# find markers for every cluster compared to all remaining cells, report only the positive ones
te.markers = FindAllMarkers(tenx.seurat, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# display top 3 for each cluster
te.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

top.genes = te.markers %>%
  group_by(cluster) %>%
  top_n(n=20, wt=avg_log2FC) %>%
  filter(!grepl(gene, pattern=c("^MT-"))) %>% # filter out mitochondrial and AC genes sequentially
  filter(!grepl(gene, pattern=c("^AC")))
top.genes %>% print(n=Inf) # to display the markers, check if they are correct
library(data.table)
fwrite(top.genes, file="heatmap.genes.txt", sep="\t")

CairoPDF(file="heatmap.pdf", width=26, height=28)
DoHeatmap(tenx.seurat, features=top.genes$gene, raster=FALSE, size=10) + 
  # NoLegend() +
  scale_fill_gradientn(colors=heatcols) +
  scale_y_discrete(position="right") +
  theme(text=element_text(size=26))
dev.off()





