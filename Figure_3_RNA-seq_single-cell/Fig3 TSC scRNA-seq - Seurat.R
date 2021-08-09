
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# load libraries ----

library(dplyr)
library(Seurat)
library(ggplot2)
library(extrafont)
library(viridis)
library(Cairo)
library(patchwork)



# load 10X datasets ----

basedir="../CELLRANGER-AGGR/"
samples = c("H9TE_p24")

tenx.data = Read10X(data.dir=paste0(basedir,samples[1],"/outs/count/filtered_feature_bc_matrix/"))

# initialize the Seurat object with the raw (non-normalized data)
tenx.seurat = CreateSeuratObject(counts=tenx.data, project=samples[1], min.cells=3, min.features=200)

# clean up
rm(tenx.data)

# merge Seurat objects
tenx.seurat = merge(x = tenx.seurat.r1, y = c(tenx.seurat.r2), project="H9TE.exp", add.cell.ids=samples)

# check number of cells from each original Seurat object
table(tenx.seurat$orig.ident)

# clean up
rm(tenx.seurat.r1, tenx.seurat.r2)



# QC -----

tenx.seurat = PercentageFeatureSet(tenx.seurat, "^MT-", col.name = "percent_mito")
tenx.seurat = PercentageFeatureSet(tenx.seurat, "^RP[SL]", col.name = "percent_ribo")
# or
tenx.seurat[["percent_mito"]] = PercentageFeatureSet(tenx.seurat, pattern = "^MT-")
tenx.seurat[["percent_ribo"]] = PercentageFeatureSet(tenx.seurat, pattern = "^RP[SL]")

feats = c("nFeature_RNA","nCount_RNA","percent_mito","percent_ribo")
CairoPDF(file=paste0("QC.pdf"), width=16, height=8)
VlnPlot(tenx.seurat, group.by="orig.ident", features=feats, pt.size=0.1,ncol=4) + NoLegend()
dev.off()
rm(feats)

# feature scatter
CairoPDF(file=paste0("feature.scatter.pdf"), width=16, height=8)
plot1 = FeatureScatter(tenx.seurat, feature1 = "nCount_RNA", feature2 = "percent_mito")
plot2 = FeatureScatter(tenx.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

# subset based on QC
tenx.seurat = subset(tenx.seurat, subset=nFeature_RNA>500 & nFeature_RNA<8000 & percent_mito<25 & percent_ribo<35)



# PRE-PROCESS ----

tenx.seurat = NormalizeData(tenx.seurat)
tenx.seurat = FindVariableFeatures(tenx.seurat)
tenx.seurat = ScaleData(tenx.seurat)
tenx.seurat = RunPCA(tenx.seurat)
tenx.seurat = RunUMAP(tenx.seurat, dim="pca", dims=1:20)
tenx.seurat = FindNeighbors(tenx.seurat)
tenx.seurat = FindClusters(tenx.seurat, resolution=0.2)
tenx.seurat = RunTSNE(tenx.seurat)



# CLUSTREE ----

tenx.seurat = FindClusters(tenx.seurat, resolution=seq(from=0.1, to=1, by=0.1))

library(clustree)
CairoPDF(file=paste0("clustree.pdf"), width=18, height=12)
clustree(tenx.seurat, prefix="RNA_snn_res.")
dev.off()



# Dimensionality reduction plots ----

# UMAP
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
ggsave(filename="umap.clusters.res0.2.pdf", device=cairo_pdf, width=6.8, height=7)

# UMAP - color by sample
plot.colors = c(viridis_pal(begin=0.2,end=0.7, option="D")(4))
UMAPPlot(tenx.seurat, group.by="orig.ident", pt.size=1.2) +
  scale_color_discrete(c=100, l=58, labels=c("r1","r2")) +
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



# JackStraw and elbow plot ----

tenx.seurat = JackStraw(object = tenx.seurat, reduction = "pca")
tenx.seurat = ScoreJackStraw(tenx.seurat, dims = 1:20)
JackStrawPlot(tenx.seurat, dims = 1:15)
ElbowPlot(object=tenx.seurat, ndims=10, reduction="pca")



# RENAME or MERGE clusters ----

levels(tenx.seurat)
new.cluster.ids = c("0", "0", "1", "0", "2", "3")
names(new.cluster.ids) = levels(tenx.seurat)
tenx.seurat = RenameIdents(tenx.seurat, new.cluster.ids)



# FEATURE PLOTS ----

sel.features = c("PAGE4","PEG10","PARP1","NR2F2","LRP2","VGLL1","VIM","COL1A2","CGA","INSL4") 

CairoPDF(file=paste("feature.plot",paste(sel.features[1:3],collapse="."),"pdf",sep="."), width=25, height=10)
feature.plots = FeaturePlot(tenx.seurat, features=sel.features, pt.size=1.2, min.cutoff="q9", order=TRUE, combine=FALSE)
feature.plots = lapply(X=feature.plots,
                       FUN=function(x) {x + theme(text=element_text(family="Noto Sans Cond"),
                                                  plot.title=element_text(size=32),
                                                  axis.title.x=element_blank(), axis.title.y=element_blank(),
                                                  legend.key.size=unit(0.3,"in"), legend.text=element_text(size=18)
                       )}
)
wrap_plots(feature.plots, ncol=5)
dev.off()

rm(sel.features)



# MARKER analysis ----

# find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers = FindAllMarkers(tenx.seurat, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# display top n for each cluster
all.markers %>% group_by(cluster) %>% top_n(n=50, wt=avg_log2FC) %>% write.table(file="cluster.markers.txt", sep="\t", row.names=FALSE)
all.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC) %>% print(n=Inf)



# HEATMAP ----

# define custom color palette
ramp = colorRamp(c("skyblue4", "beige", "firebrick4"))
heatcols = rgb(ramp(seq(0,1,length=256)),max=256)

# to re-order clusters on the heatmap so that they are grouped by original sample
# re-level the idents of the tenx.seurat object
levels(tenx.seurat) = unlist(cluster.groups)
# authors: the development version of Seurat contains a new argument in DoHeatmap group.order
# as of now in stable version 3.1.3, it is not implemented

top.genes = all.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC) %>%
  # arrange (factor(cluster, levels=unlist(cluster.groups)))
top.genes %>% print(n=Inf)
library(data.table)
fwrite(top.genes, file="heatmap.genes.txt", sep="\t")

# downsampling of the Seurat object will make the columns for each cluster on the heatmap the same size
CairoPDF(file="heatmap.custom.colors.pdf", width=16, height=24)
# CairoPNG(file="heatmap.custom.colors.png", width=8400, height=4200) # this can get around colors being muted upon import of PDF into Photoshop, use size=40 in DoHeatmap and theme text size 80, however, the top margin seems too wide
DoHeatmap(tenx.seurat, features=top.genes$gene, raster=FALSE, size=10) + # subset(tenx.seurat, downsample=100)
  NoLegend() +
  scale_fill_gradientn(colors=heatcols) +
  scale_y_discrete(position="right") +
  theme(text=element_text(size=26))
dev.off()




