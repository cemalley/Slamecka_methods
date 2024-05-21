
# SLINGSHOT ----

library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(SingleCellExperiment)
library(scater)
library(RColorBrewer)


## prepare data for slingshot analysis ----

library(scry)
sce = as.SingleCellExperiment(tenx.seurat)



## explore the SCE object

dim(reducedDim(sce, "UMAP"))
plotReducedDim(sce, dimred="UMAP", colour_by="orig.ident")
plotReducedDim(sce, dimred="UMAP", colour_by="seurat_clusters")



## run slignshot ----

# the output is a SingleCellExperiment object with slingshot results incorporated
# all of the results are stored in a PseudotimeOrdering object, which is added to the colData of the original object
# and can be accessed via colData(sce)$slingshot
sce = slingshot(sce, reducedDim = "UMAP",
                clusterLabels = colData(sce)$seurat_clusters,
                stretch=0, extend="n", maxit=30)

library(TrajectoryUtils)



## UMAP plots with pseudotimes and trajectory lines ----

lapply(grep(colnames(colData(sce)), pattern="slingPseudotime_", value=TRUE), function(X) {
  require(grDevices)
  require(slingshot)
  require(BUSpaRse)
  require(tidyverse)
  require(tidymodels)
  require(Seurat)
  require(scales)
  require(viridis)
  require(Matrix)
  require(SingleCellExperiment)
  require(scater)
  require(RColorBrewer)
  CairoPDF(file=paste0("slingshot.trajectories.",X,".pdf"), width=10, height=10)
  colors = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  plotcol = colors[cut(sce[[X]], breaks=100)]
  plotcol[is.na(plotcol)] = "grey" # assign grey color to NA values - values that do not belong to the lineage
  plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2, col='black')
  dev.off()
})

## save the original colors
# 11 is the maximum number of colors"
CairoPDF(file="pseudotime.colors.pdf", width=10, height=3)
display.brewer.pal(11, "Spectral")
dev.off()


