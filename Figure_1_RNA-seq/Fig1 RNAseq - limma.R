
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# LIMMA ----


# create a file name base ----

fn = list()
fn$base = "star-featureCounts"



# TMM normalization ----

DGE.fn = calcNormFactors(DGE.f)
rm(DGE.f)
plotMDS(DGE.fn)
# SUBSET the DGE now if not subset before
# DGE.fn = DGE.fn[ ,as.character(targets$SampleName)]



# MDS plot ----

library(extrafont)
library(Cairo)
CairoPDF(file=paste0("MDS.",fn$base,".symbols.pdf"), width=4, height=4.4)
par(family="Noto Sans Cond") # cex=1.2, lwd=2 - to control text and point size; and line width
plotMDS(DGE.fn, pch=targets$PlotPCH, bg=as.character(targets$PlotColor), col="white", cex=2, lwd=0.8)
# legend(x="left", pch=unique(targets$PlotPCH), legend=unique(targets$Line),
#        inset=0.01, pt.cex=1.2, cex=1.2, text.width=1.4, y.intersp = 0.8#,
#        # bty="n" # comment for box  around the legend - default value is "o"
#        )
legend(x="bottomleft", fill=unique(targets$PlotColor)[c(2,1,3,4)], legend=unique(targets$Factor)[c(2,1,3,4)],
      border="white", inset=0.01, pt.cex=1.2, cex=1.2, text.width=1.4, y.intersp = 0.8#,
       # bty="n" # comment for box  around the legend - default value is "o"
       )
dev.off()



# design matrix ----

factor = DGE.fn$samples$group
design = model.matrix(~0+factor)
colnames(design) = levels(factor)



# voom transformation ----

pdf(file=paste0("voom.",fn$base,".pdf"), width=5, height=5)
DGE.fnv = voom(DGE.fn, design, plot=TRUE)
dev.off()
# rm(DGE.fn)



# lmFit ----
fit = lmFit(DGE.fnv, design)

# choose contrasts
fn$contrast = paste0(fn$base,".time-course")
contrast.matrix = makeContrasts(D10-D0, D7-D0, D3-D0,
                                levels=design)
# fit and eBayes smoothing
fit = contrasts.fit(fit, contrast.matrix)
fit = eBayes(fit)



# FILTERING, ORDERING based on P-value ----

# re-order all genes based on evidence of differential expression
top = topTable(fit, adjust='fdr', number=Inf)
# write results to file
write.table(top,file=paste0("DGE.",fn$contrast,".txt"), sep="\t", quote=FALSE, col.names=NA)



# choose CUT-OFFs for MARKER list ----

CO = list("logFC"=list("D10"=8.00,
                       "D7"=8.00,
                       "D3"=8.00,
                       "D0"=-7.50),
          "adjPval"=list("global"=1e-4))
# filter based on FDR-adjusted P value first
top.pval = top[top$adj.P.Val < CO$adjPval$global, ]
# make a named list of markers for each group
markers = list("D10"=top.pval[top.pval$D10...D0 > CO$logFC$D10 & !(top.pval$hgnc_symbol == "") , ]$hgnc_symbol,
               "D7"=top.pval[top.pval$D7...D0 > CO$logFC$D7 & !(top.pval$hgnc_symbol == "") , ]$hgnc_symbol,
               "D3"=top.pval[top.pval$D3...D0 > CO$logFC$D3 & !(top.pval$hgnc_symbol == "") , ]$hgnc_symbol,
               "D0"=top.pval[top.pval$D10...D0 < CO$logFC$D0 & !(top.pval$hgnc_symbol == "") , ]$hgnc_symbol
               )
# get number of genes in each marker group
sapply(markers, length)
rm(top.pval)

# write cutoffs into file
writeLines(capture.output(CO), con=paste0("DGE.",fn$base,".CO.file.txt"))

# add filter of the top table to the "fn" filename list ----
fn$selected = paste0(fn$base,".log",CO$logFC$global,".p",CO$adjPval$global)
fn$selected = paste0(fn$base,".D10.v.D0.log",CO$logFC$D10,".p",CO$adjPval$global)
fn$selected = paste0(fn$base,".D0.v.D10.log",CO$logFC$D0,".p",CO$adjPval$global)
fn$selected = paste0(fn$base,".cutoffs.in.CO.file")



# HEATMAP ----

heatmap = list()
# choose overall font size for heatmap labels
heatmap$fontsize = 22
heatmap$fontfamily = "" # default font, size 22 works well



# choose features for heatmap ----

# top.sel genes
heatmap$features = top.sel$ensembl_gene_id
# top DE genes in each marker group
heatmap$features = markers$D10
# further restrict these genes manually for final panel in Figure 1
markers$D10
markers$manual = list()
markers$manual$trophoblast = c("IGFBP3","H19","CDX2","GATA2","GATA3",
                               "TFAP2A","TFAP2B","HAND1","MSX2","PAPPA2",
                               "EPAS1","STS","GABRP","VGLL1","XAGE2",
                               "PAGE4","SIGLEC6","PPARG","TBX3","XAGE2",
                               "GABRP","EPAS1","MSX2","IGF2","STS",
                               "MFAP5","PLAC4","DIO3", "SP6", "CRH",
                               "TFAP2C","NR2F2","LRP2","KRT7","KRT18") # this row of genes requested by Ilyas
markers$manual$syncytiotrophoblast = c("ERVFRD-1","HOPX","ERVV-2","ERVW-1","CGA",
                                       "KRT23","ERVV-1","GCM1","ERVE-1","DLX3",
                                       "HSD3B1","CGB3","CGB8","CYP19A1","MUC15",
                                       "CSF3R","S100P","INSL4","LGALS16",
                                       "TEAD3","SDC1") # this row of genes requested by Ilyas
heatmap$features = c(markers$manual$trophoblast, markers$manual$syncytiotrophoblast)
# all markers

# hPSC D0 markers
heatmap$features = markers$D0


  
# SUBSET for heatmap ----

# subsetting based on Ensembl gene IDs and symbols
heatmap$DGE = DGE.fnv[DGE.fnv$genes$hgnc_symbol %in% heatmap$features | DGE.fnv$genes$ensembl_gene_id %in% heatmap$features, ]



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
heatmap$matrix = heatmap$DGE$E + abs(min(heatmap$DGE$E)) # adding the absolute value of the (global) minimum to avoid negative values
rownames(heatmap$matrix) = heatmap$DGE$genes$external_gene_name
colnames(heatmap$matrix) = colnames(heatmap$DGE)



# heatmap ANNOTATION - top colors ----

col.df = unique(data.frame(targets$Factor,as.character(targets$PlotColor)))
col.list = list(Stage=col.df[,2]) # list named "CellType" containing the colors
names(col.list$Stage) = col.df[,1] # name the elements of the list with the factor elements
rm(col.df)
heatmap$col.annotation = HeatmapAnnotation(df = data.frame(Stage = heatmap$DGE$targets$group),
                                           col = col.list,
                                           annotation_legend_param=list(labels_gp=gpar(fontsize=heatmap$fontsize, fontfamily=heatmap$fontfamily), # modify legend label font size
                                                                        title_gp=gpar(fontsize=heatmap$fontsize, fontfamily=heatmap$fontfamily), # legend title font size
                                                                        grid_height=unit(0.8,"cm"), grid_width=unit(0.8,"cm") # legend tile dimensions
                                                                        )
                                           )
rm(col.list)




# SLICE heatmap based on feature groups ----
make.heatmap.feature.groups = function(genes.df, # data.frame with genes
                                       col.to.use = 2, # which column of the data.frame contains the gene categories
                                       heatmap.dge=heatmap$DGE, # subset of the DGE.fnv object created for heatmap
                                       group.order=NULL # order of feature groups to display on heatmap, from top to bottom
                                       ) {
  if (class(genes.df) == "list") {
    library(reshape2)
    heatmap.feature.groups = reshape2::melt(genes.df)
  } else if (class(genes.df) == "data.frame") {
    heatmap.feature.groups = genes.df
  } else print("input must be a named list or a data.frame")
  heatmap.feature.groups = heatmap.feature.groups[match(heatmap.dge$genes$external_gene_name, heatmap.feature.groups[ ,1]), ] # assuming gene symbol is in the first column
  if (is.null(group.order)) {return(heatmap.feature.groups[ ,col.to.use])
    } else if (all(group.order %in% unique(heatmap.feature.groups[ ,col.to.use]))) { # all evaluates into a single value
      heatmap.feature.groups[ ,col.to.use] = factor(heatmap.feature.groups[ ,col.to.use], levels=group.order)
      return(heatmap.feature.groups[ ,col.to.use])
    } else {
      print("group name mismatch!")
      return(NULL)
        }
}

# construct the feature groups
# for final panel in Figure 1
heatmap$feature.groups = make.heatmap.feature.groups(markers$manual,
                                                     2,
                                                     heatmap$DGE,
                                                     group.order=c("trophoblast","syncytiotrophoblast"))



# DENDROGRAM ----

# get a global dendrogram to reuse if plotting smaller groups of genes
# get dendrogram from a simple DGE test
get.col.dendro = function (dge=DGE.fnv, # DGEList - filtered, normalized and voom-transformed
                           top.df=top, # top table data.frame
                           col.to.use, # column with contrast to use
                           co.logFC, # logFC cut-off
                           co.P) { # P value cut-off
  heatmap = list()
  heatmap$features = top.df[top.df[col.to.use] > co.logFC & top.df$adj.P.Val < co.P, ]$ensembl_gene_id
  heatmap$DGE = dge[dge$genes$ensembl_gene_id %in% heatmap$features, ]
  heatmap$matrix = heatmap$DGE$E
  rownames(heatmap$matrix) = heatmap$DGE$genes$external_gene_name
  colnames(heatmap$matrix) = colnames(heatmap$DGE)
  col.dendro = column_dend(Heatmap(heatmap$matrix))
  return(col.dendro)
}
heatmap$col.dendro = get.col.dendro(col.to.use="D10...D0", co.logFC=2, co.P=1e-04)
# re-order samples on the dendrogram
library(dendextend)
labels(heatmap$col.dendro) # get order of labels of the dendrogram
heatmap$col.dendro = rotate(heatmap$col.dendro, order=c(10:12,7:9,4:6,1:3)) # reorder dendrogram and re-write the slot in heatmap list
# save
saveRDS(heatmap$col.dendro, file="col.dendro.rds")



# DRAW heatmaps and save PDFs ----

# color key legend parameters are modified within the Heatmap function
# annotation legend parameters are modified above in the HeatmapAnnotation function
CairoPDF(file=paste0("heatmap.",fn$contrast,".",nrow(heatmap$matrix),"features.pdf"), width=10, height=17) # for final panel in Figure 1
CairoPDF(file=paste0("heatmap.",fn$selected,".",nrow(heatmap$matrix),".hPSC.markers.pdf"), width=12, height=32)
heatmap$heatmap = Heatmap(heatmap$matrix,
                          col=heatmap$colors,
                          top_annotation=heatmap$col.annotation,
                          column_dend_height=unit(3,"cm"), row_dend_width=unit(3,"cm"),
                          row_names_gp=gpar(fontsize=heatmap$fontsize, fontfamily=heatmap$fontfamily),
                          column_names_gp=gpar(fontsize=12, fontfamily=heatmap$fontfamily),
                          row_names_max_width = unit(8,"cm"),
                          cluster_columns=heatmap$col.dendro,
                          row_split=heatmap$feature.groups, cluster_row_slices=FALSE, row_title_gp=gpar(fontsize=heatmap$fontsize), row_gap=unit(3,"mm"), # turn off row slice clustering to manually order slices based on previously specified order of groups using function make.heatmap.feature.groups
                          heatmap_legend_param=list(labels_gp=gpar(fontsize=heatmap$fontsize, fontfamily=heatmap$fontfamily), # font size of legend labels
                                                    title=NULL, # remove title
                                                    #legend_height=unit(4,"cm"),# overall height of the legend, if it needs to be wider, use the parameters below
                                                    grid_height=unit(1,"cm"), grid_width=unit(1,"cm"))
                          )
draw(heatmap$heatmap, heatmap_legend_side="right", show_annotation_legend=FALSE, merge_legends=FALSE)
dev.off()

