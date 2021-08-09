
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# LIMMA ----

# no subsetting/partitioning
tar.tmp = targets
DGE.f.subset = DGE.f[ ,tar.tmp$SampleID]



# create a file name base ----

fn = list()
fn$base = "counts"



# TMM normalization ----
DGE.fn = calcNormFactors(DGE.f.subset)
rm(DGE.f.subset)
plotMDS(DGE.fn)



# MDS plot ----

# set order of cell types and cell lines to appear on the legend
cols.df = unique(tar.tmp[ ,c("CellType", "PlotColor")])
cols.df = cols.df[c(1,2,3), ]
points.df = unique(tar.tmp[ ,c("Line", "PlotPCH")])

library(extrafont)
library(Cairo) # will embed fonts automatically, no need to use embed_fonts() from "extrafont" package
CairoPDF(file=paste0("MDS.",fn$base,".symbols.pdf"), width=4, height=4.4) # family="Noto Sans Cond" - as an argument doesn't seem to work here, setting par in the next line does
par(family="Noto Sans Cond") # cex=1.2, lwd=2 - to control text and point size; and line width
plotMDS(DGE.fn, pch=tar.tmp$PlotPCH, bg=tar.tmp$PlotColor, col="white", cex=2, lwd=0.8)
legend(x="topright", fill=cols.df$PlotColor, legend=cols.df$CellType,
       border="white", inset=0.01, pt.cex=1.2, cex=1.2, text.width=1, y.intersp = 0.8, x.intersp=0.2#,
       # bty="n" # comment for box  around the legend - default value is "o"
       )
dev.off()

rm(cols.df, points.df)



# design matrix ----

factor = DGE.fn$samples$group
design = model.matrix(~0+factor)
colnames(design) = levels(factor)

# voom transformation ----

pdf(file=paste0("voom.",fn$base,".pdf"), width=5, height=5)
DGE.fnv = voom(DGE.fn, design, plot=TRUE)
dev.off()
rm(DGE.fn)

# lmFit ----

fit = lmFit(DGE.fnv, design)

# choose contrasts ----

fn$contrast = paste0(fn$base,".TSC-TE.p0.TE.p0-hPSC.TSC-hPSC")
contrast.matrix = makeContrasts(TSC.JHU198i-TE.p0.JHU198i,TE.p0.JHU198i-hPSC.JHU198i,TSC.JHU198i-hPSC.JHU198i,
                                levels=design)

# fit and eBayes smoothing ----

fit = contrasts.fit(fit, contrast.matrix)
fit = eBayes(fit)



# FILTERING, ORDERING based on P-value ----

# re-order all genes based on evidence of differential expression
top = topTable(fit, adjust='fdr', number=Inf)
# write results to file
write.table(top,file=paste0("DGE.",fn$contrast,".txt"), sep="\t", quote=FALSE, col.names=NA)



# GENES OF INTEREST ----

C19MC = readLines("C19MC.source2.txt")
C19MC = paste0("hsa-",C19MC)
table(C19MC %in% top$Mature_ID)



# HEATMAP ----

heatmap = list()
# choose overall font size for heatmap labels
heatmap$fontsize = 22
heatmap$fontfamily = "Noto Sans Cond"
heatmap$fontfamily = "" # default font



# choose features for heatmap ----

heatmap$features = C19MC



# SUBSET for heatmap ----

heatmap$DGE = DGE.fnv[DGE.fnv$genes$Mature_ID %in% heatmap$features, ]



# attach libraries ----

library(ComplexHeatmap)
library(circlize)
library(extrafont)
library(Cairo)



# choose color scale ----

ramp = colorRamp(c("skyblue4", "beige", "firebrick4"))
heatmap$colors = rgb(ramp(seq(0,1,length=256)),max=256)
rm(ramp)



# create the heatmap expression matrix ----

# include gene symbols
heatmap$matrix = heatmap$DGE$E + abs(min(DGE.fnv$E)) # adding the absolute value of the (global) minimum to avoid negative values
rownames(heatmap$matrix) = heatmap$DGE$genes$Mature_ID
colnames(heatmap$matrix) = colnames(heatmap$DGE)
# keep unique rows
heatmap$matrix = unique(heatmap$matrix)



# heatmap ANNOTATION - side colors ----
# can be column and row side colors
# get levels of the group factor and manually assign colors
levels(heatmap$DGE$targets$group)
# construct from colors assigned earlier in targets
# first, create a data.frame from factor and assigned colors
# use unique to remove duplicate rows
col.df = unique(data.frame(paste(as.character(tar.tmp$CellType),as.character(tar.tmp$Line),sep="."),as.character(tar.tmp$PlotColor)))
col.list = list(CellType=col.df[,2]) # list named "CellType" containing the colors
names(col.list$CellType) = col.df[,1] # name the elements of the list with the factor elements
rm(col.df)
heatmap$col.annotation = HeatmapAnnotation(df = data.frame(CellType = heatmap$DGE$targets$group),
                                           col = col.list,
                                           annotation_legend_param=list(labels_gp=gpar(fontsize=heatmap$fontsize, fontfamily=heatmap$fontfamily), # modify legend label font size
                                                                        title_gp=gpar(fontsize=heatmap$fontsize, fontfamily=heatmap$fontfamily), # legend title font size
                                                                        grid_height=unit(0.8,"cm"), grid_width=unit(0.8,"cm") # legend tile dimensions
                                                                        )
                                           )
rm(col.list)



# DENDROGRAM ----

# get a global dendrogram to reuse if plotting smaller groups of genes
# get dendrogram from a simple DGE test
get.col.dendro = function (dge=DGE.fnv, # DGEList - filtered, normalized and voom-transformed
                           top.df=top, # top table data.frame
                           col.to.use, # column with contrast to use
                           co.logFC, # logFC cut-off
                           co.P) { # P value cut-off
  heatmap = list()
  heatmap$features = top.df[top.df[col.to.use] > co.logFC & top.df$adj.P.Val < co.P, ]$Mature_ID
  heatmap$DGE = dge[dge$genes$Mature_ID %in% heatmap$features, ]
  heatmap$matrix = heatmap$DGE$E
  rownames(heatmap$matrix) = heatmap$DGE$genes$Mature_ID
  colnames(heatmap$matrix) = colnames(heatmap$DGE)
  col.dendro = column_dend(Heatmap(heatmap$matrix))
  return(col.dendro)
}
heatmap$col.dendro = get.col.dendro(col.to.use="TSC.JHU198i...TE.p0.JHU198i", co.logFC=2, co.P=1e-03) # CC samples inlcuded
# re-order samples on the dendrogram - no JHU191i
labels(rotate(heatmap$col.dendro, order=c(5:7,8:10,1:4)))
heatmap$col.dendro = rotate(heatmap$col.dendro, order=c(5:7,8:10,1:4))
# save and load
saveRDS(heatmap$col.dendro, file="col.dendro.rds")
col.dendro = readRDS("col.dendro.rds")
heatmap$col.dendro = col.dendro
rm(col.dendro)



# DRAW heatmap and save PDF ----

# color key legend parameters are modified within the Heatmap function
# annotation legend parameters are modified above in the HeatmapAnnotation function
heatmap$pdf.w = 12

CairoPDF(file=paste0("heatmap.",fn$base,".C19MC.",nrow(heatmap$matrix),"features.pdf"), width=heatmap$pdf.w, height=18)
heatmap$heatmap = Heatmap(heatmap$matrix,
                          col=heatmap$colors,
                          top_annotation=heatmap$col.annotation,
                          column_dend_height=unit(3,"cm"), row_dend_width=unit(3,"cm"),
                          row_names_gp=gpar(fontsize=heatmap$fontsize, fontfamily=heatmap$fontfamily),
                          column_names_gp=gpar(fontsize=12, fontfamily=heatmap$fontfamily),
                          row_names_max_width = unit(8,"cm"),
                          cluster_columns=heatmap$col.dendro,
                          heatmap_legend_param=list(labels_gp=gpar(fontsize=heatmap$fontsize, fontfamily=heatmap$fontfamily), # font size of legend labels
                                                    title=NULL, # remove title
                                                    #legend_height=unit(4,"cm"),# overall height of the legend, if it needs to be wider, use the parameters below
                                                    grid_height=unit(1,"cm"), grid_width=unit(1,"cm"))
                          )
draw(heatmap$heatmap, heatmap_legend_side="right", show_annotation_legend=FALSE, merge_legends=FALSE)
dev.off()


