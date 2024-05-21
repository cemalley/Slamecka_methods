
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# LIMMA ----



# subset  ----

# subset dataset based on criteria for the purpose of creating an MDS plot without choriocarcinoma cell lines
tar.tmp = targets[!(targets$CellType == "CC"), ] # exclude choriocarcinoma lines
DGE.f.subset = DGE.f[ ,tar.tmp$SampleName.LIMS]

# no subsetting
tar.tmp = targets
DGE.f.subset = DGE.f



# create a file name base ----

fn = list()
fn$base = "bbmap-featureCounts"



# TMM normalization ----

DGE.fn = calcNormFactors(DGE.f.subset)
rm(DGE.f.subset)
plotMDS(DGE.fn)



# MDS plot ----

library(extrafont)
library(Cairo) # will embed fonts automatically, no need to use embed_fonts() from "extrafont" package
CairoPDF(file=paste0("MDS.",fn$base,".symbols.pdf"), width=4, height=4.4)
par(family="Noto Sans Cond") # cex=1.2, lwd=2 - to control text and point size; and line width
plotMDS(DGE.fn, pch=tar.tmp$PlotPCH, bg=as.character(tar.tmp$PlotColor), col="white", cex=2, lwd=0.8)
legend(x="left", pch=unique(tar.tmp$PlotPCH), legend=unique(tar.tmp$Line),
       inset=0.01, pt.cex=1.2, cex=1.2, text.width=1.4, y.intersp = 0.8#,
       # bty="n" # comment for box  around the legend - default value is "o"
       )
legend(x="topright", fill=unique(tar.tmp$PlotColor), legend=unique(tar.tmp$CellType),
       border="white", inset=0.01, pt.cex=1.2, cex=1.2, text.width=1.4, y.intersp = 0.8#,
       # bty="n" # comment for box  around the legend - default value is "o"
       )
dev.off()

# noCC (legend not needed)
CairoPDF(file=paste0("MDS.",fn$base,".symbols.pdf"), width=4, height=4.4)
par(family="Noto Sans Cond")
plotMDS(DGE.fn, pch=tar.tmp$PlotPCH, bg=as.character(tar.tmp$PlotColor), col="white", cex=2, lwd=0.8)
# legend(x="top", pch=unique(tar.tmp$PlotPCH), legend=unique(tar.tmp$Line),
#        inset=0.01, pt.cex=1.2, cex=1.2, text.width=1.4, y.intersp = 0.8#,
#        # bty="n" # comment for box  around the legend - default value is "o"
#        )
# legend(x="bottomright", fill=unique(tar.tmp$PlotColor), legend=unique(tar.tmp$CellType),
#        border="white", inset=0.01, pt.cex=1.2, cex=1.2, text.width=1.4, y.intersp = 0.8#,
#        # bty="n" # comment for box  around the legend - default value is "o"
#        )
dev.off()



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

fn$contrast = paste0(fn$base,".TE-PSC.CC-all.Ecto-PSC,Meso-PSC,Endo-PSC.PSC-all") # CC included
contrast.matrix = makeContrasts(TE.WA09-hPSC.WA09, TE.WA14-hPSC.WA14, TE.WA17-hPSC.WA17,
                                Ecto.WA09-hPSC.WA09, Meso.WA09-hPSC.WA09, Endo.WA09-hPSC.WA09,
                                TE.vs.hPSC = (TE.WA09+TE.WA14+TE.WA17)/3-(hPSC.WA09+hPSC.WA14+hPSC.WA17)/3,
                                hPSC.vs.all = (hPSC.WA09+hPSC.WA14+hPSC.WA17)/3-(TE.WA09+TE.WA14+TE.WA17+Ecto.WA09+Meso.WA09+Endo.WA09+CC.BeWo+CC.JEG3)/8, # CC included as well
                                CC.vs.all = (CC.BeWo+CC.JEG3)/2-(TE.WA09+TE.WA14+TE.WA17+Ecto.WA09+Meso.WA09+Endo.WA09+hPSC.WA09+hPSC.WA14+hPSC.WA17)/9, # comment this line if CC are excluded
                                levels=design)

# fit and eBayes smoothing ----

fit = contrasts.fit(fit, contrast.matrix)
fit = eBayes(fit)



# FILTERING, ORDERING based on P-value ----

# re-order all genes based on evidence of differential expression
top = topTable(fit, adjust='fdr', number=Inf)
# write results to file
write.table(top,file=paste0("DGE.",fn$contrast,".txt"), sep="\t", quote=FALSE, col.names=NA)



# choose CUT-OFFs for MARKER list ----

CO = list("logFC"=list("global"=8.00,
                       "hPSC"=8.00,
                       "TE"=8.00,
                       "Ecto"=7.00,
                       "Meso"=9.00,
                       "Endo"=9.00,
                       "CC"=8.00),
          "adjPval"=list("global"=1e-4))

# filter based on FDR-adjusted P value first
top.pval = top[top$adj.P.Val < CO$adjPval$global, ]
# make a named list of markers for each group
markers = list("hPSC"=top.pval[top.pval$hPSC.vs.all > CO$logFC$hPSC & !(top.pval$hgnc_symbol == ""), ]$hgnc_symbol, # external_gene_name
               "TE"=top.pval[top.pval$TE.vs.hPSC > CO$logFC$TE & !(top.pval$hgnc_symbol == ""), ]$hgnc_symbol, # external_gene_name
               "Ecto"=top.pval[top.pval$Ecto.WA09...hPSC.WA09 > CO$logFC$Ecto & !(top.pval$hgnc_symbol == ""), ]$hgnc_symbol, # external_gene_name
               "Meso"=top.pval[top.pval$Meso.WA09...hPSC.WA09 > CO$logFC$Meso & !(top.pval$hgnc_symbol == ""), ]$hgnc_symbol, # external_gene_name
               "Endo"=top.pval[top.pval$Endo.WA09...hPSC.WA09 > CO$logFC$Endo & !(top.pval$hgnc_symbol == ""), ]$hgnc_symbol, # external_gene_name
               "CC"=top.pval[top.pval$CC.vs.all > CO$logFC$CC & !(top.pval$hgnc_symbol == ""), ]$hgnc_symbol) # external_gene_name
# get number of genes in each marker group
sapply(markers, length)
rm(top.pval)

# write cutoffs into file
writeLines(capture.output(CO), con=paste0("DGE.",fn$base,".CO.file.txt"))



# add filter of the top table to the "fn" filename list ----
fn$selected = paste0(fn$base,".log",CO$logFC$global,".p",CO$adjPval$global)
fn$selected = paste0(fn$base,".TE.log",CO$logFC$TE,".p",CO$adjPval$global)
fn$selected = paste0(fn$base,".cutoffs.in.CO.file")




# GENES OF INTEREST ----

basepath = "H:/"

library(tidyverse)
library(readxl)
library(data.table)

GOI = list()
# PA = Protein Atlas - enriched/specific genes
GOI$placenta.PA.enriched = read.delim(paste0(basepath,"DATA ANALYSIS/NGS/GOI/[2020-4-6] Protein Atlas - placenta-specific genes (91).tsv"),header=TRUE,sep="\t",comment.char="#")
# [2019, Nature] Reconstituting the transcriptome and DNA methylome landscapes of human implantation
GOI$scTrMet = read_xlsx(paste0(basepath,"ARTICLES/PSC Differentiation/Trophectoderm/[2019, Nature] Reconstituting the transcriptome and DNA methylome landscapes of human implantation - SUPP TABLES/Supplementary Table 3 Lineage_Markers.xlsx"))
GOI$scTrMet = GOI$scTrMet[ ,c(1,7,2:6)] # reorder columns, function make.feature.groups expects symbols in column 1 and group in column 2
GOI$scTrMet = as.data.frame(GOI$scTrMet)



# HEATMAP ----

heatmap = list()
# choose overall font size for heatmap labels
heatmap$fontsize = 24
heatmap$fontfamily = "" # default font



# choose features for heatmap ----

# top DE genes - TE markers and PSC, 3L, CC markers separately
heatmap$features = markers$TE
heatmap$features = c(markers$hPSC,
                     markers$Ecto,
                     markers$Meso,
                     markers$Endo,
                     markers$CC)

# GOI - genes of interest

heatmap$features = as.character(GOI$placenta.PA.enriched$Gene)
heatmap$features = as.character(GOI$scTrMet[GOI$scTrMet$cluster == "TE", ]$gene)



# SUBSET for heatmap ----

# subsetting based on either Ensembl gene IDs and symbols
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
heatmap$matrix = heatmap$DGE$E + abs(min(DGE.fnv$E)) # adding the absolute value of the (global) minimum to avoid negative values
rownames(heatmap$matrix) = heatmap$DGE$genes$external_gene_name
colnames(heatmap$matrix) = colnames(heatmap$DGE)



# heatmap ANNOTATION - side colors ----

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



# heatmap ANNOTATION - mark annotation ----

# highlight selected genes on a crowded heatmap
# anno_mark() needs indices of the positions and the corresponding labels
# choose genes to highlight
writeLines(markers$TE, con="selected.genes.txt") # write into txt and manually select
selected.genes = readLines("selected.genes.txt") # read back in
# selected.genes = c("") # manual
selected.genes = unique(selected.genes)
selected.genes = heatmap$DGE$genes$external_gene_name[heatmap$DGE$genes$external_gene_name %in% selected.genes]
heatmap$mark.annotation = HeatmapAnnotation(placenta.genes=anno_mark(at=match(selected.genes, heatmap$DGE$genes$external_gene_name),
                                                                     labels=selected.genes,
                                                                     labels_gp=gpar(fontsize=heatmap$fontsize, fontfamily=heatmap$fontfamily),
                                                                     link_gp=gpar(lwd=2.4),
                                                                     link_width=unit(20,"mm")),
                                            which="row"
                                            )
rm(selected.genes)



# SLICE heatmap based on feature groups ----

make.heatmap.feature.groups = function(genes.df, # data.frame with genes
                                       col.to.use = 2, # which column of the data.frame contains the gene categories
                                       heatmap.dge=heatmap$DGE, # subset of the DGE.fnv object created for heatmap
                                       group.order=NULL # order of feature groups to display on heatmap, from top to bottom
                                       ) {
  # input - named list or data.frame containing gene groups
  # if input is data.frame, it should have the gene symbols in column 1 and categories/groups in column 2
  # if input is a list, it will first be converted to data.frame with genes in first column and gene groups in second
  # heatmap.dge$genes$external_gene_name - gene symbols from the DGE object that was subset to create the heatmap
  # the function will reorder the data.frame so that the genes match the gene symbols of the subset DGE for heatmap
  # finally, it will return a character vector with gene groups (second column of the data.frame)
  # this will be fed as categorical variable to Heatmap function to slice the heatmap with
  # using is.list() instead of class(input)=="list" won't work because a data.frame is also a list
  # using inherits() can fix this but it's more cumbersome
  if (class(genes.df) == "list") {
    library(reshape2)
    heatmap.feature.groups = reshape2::melt(genes.df)
  } else if (class(genes.df) == "data.frame") {
    heatmap.feature.groups = genes.df
  } else print("input must be a named list or a data.frame")
  heatmap.feature.groups = heatmap.feature.groups[match(heatmap.dge$genes$external_gene_name, heatmap.feature.groups[ ,1]), ] # assuming gene symbol is in the first column
  if (is.null(group.order)) {return(heatmap.feature.groups[ ,col.to.use])
    } else if (all(group.order %in% unique(heatmap.feature.groups[ ,col.to.use]))) { # all evaluates into a single value
      # if two conditions used, && operator must be used instead of & for proper conditional evaluation
      heatmap.feature.groups[ ,col.to.use] = factor(heatmap.feature.groups[ ,col.to.use], levels=group.order)
      return(heatmap.feature.groups[ ,col.to.use])
    } else {
      print("group name mismatch!")
      return(NULL)
        }
}

# construct the feature groups
# specify the group.order for the heatmap if needed
heatmap$feature.groups = make.heatmap.feature.groups(GOI$scTrMet, "cluster", heatmap$DGE)
# delete the TE slot from the markers list in order to exclude those genes (and plot them separately)
markers$TE = NULL
heatmap$feature.groups = make.heatmap.feature.groups(markers, 2, heatmap$DGE,
                                                     group.order=c("hPSC","Ecto","Endo","Meso","CC"))



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
heatmap$col.dendro = get.col.dendro(col.to.use="hPSC.vs.all", co.logFC=2, co.P=1e-03) # CC samples inlcuded
# save and load
saveRDS(heatmap$col.dendro, file="col.dendro.rds")
heatmap$col.dendro = readRDS("col.dendro.rds")



# DRAW heatmaps and save PDFs ----

# color key legend parameters are modified within the Heatmap function
# annotation legend parameters are modified above in the HeatmapAnnotation function
CairoPDF(file=paste0("heatmap.",fn$selected,".pdf"), width=24, height=48) # no TE feature group
CairoPDF(file=paste0("heatmap.",fn$selected,".",nrow(heatmap$matrix),"features.pdf"), width=22, height=42)
CairoPDF(file=paste0("heatmap.",fn$selected,".",nrow(heatmap$matrix),"features.mark-a.pdf"), width=22, height=20) # for top DEG - TE markers - with MARK ANNOTATION
CairoPDF(file=paste0("heatmap.",fn$base,".placental.genes.PA.enriched.91.pdf"), width=24, height=28)
CairoPDF(file=paste0("heatmap.",fn$base,".placental.genes.scTrMet.TE.pdf"), width=24, height=34)

heatmap$heatmap = Heatmap(heatmap$matrix,
                          col=heatmap$colors,
                          top_annotation=heatmap$col.annotation,
                          right_annotation=heatmap$mark.annotation, show_row_names=FALSE, # set to TRUE to check if the labels were matched properly
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



# EXPRESSION plots ----

expr.plot = function(genes.to.plot,
                     plot.h.f=6, # plot height factor
                     plot.w.f=4.6, # plot width factor
                     facet.nrow=2, # number of facet rows
                     facet.ncol=5, # number of facet columns
                     counts=FALSE, log.tr=FALSE, # use normalized counts, either intact or log2-transformed, instead of transformed values
                     group.order=NULL # order of groups on plot
                     ) {
  # subset DGE based on genes to plot
  # use DGE.fnv for transformed values
  # use DGE.fn for normalized counts - with argument counts=TRUE
  if (counts) {
    expr.dge = DGE.fn[DGE.fn$genes$external_gene_name %in% genes.to.plot | DGE.fn$genes$ensembl_gene_id %in% genes.to.plot, ]
    if (log.tr) {
      expr.matrix = log2(expr.dge$counts + 1)
      plot.scales="fixed"
      } else {
        expr.matrix = expr.dge$counts
        plot.scales="free"
        }
    } else {
    expr.dge = DGE.fnv[DGE.fnv$genes$external_gene_name %in% genes.to.plot | DGE.fnv$genes$ensembl_gene_id %in% genes.to.plot, ]
    expr.matrix = expr.dge$E + abs(min(expr.dge$E))
    plot.scales="fixed"
    }
  # name rows with genes
  rownames(expr.matrix) = expr.dge$genes$external_gene_name
  # transpose the matrix
  expr.matrix = t(expr.matrix)
  # re-order columns based on genes.to.plot so that the order of plotted facets is as the original order of genes supplied
  # first subset genes.to.plot to keep only those that were found in the DGE object
  # otherwise, the subsetting will be out of bounds
  genes.to.plot = genes.to.plot[genes.to.plot %in% colnames(expr.matrix)]
  expr.matrix = expr.matrix[ ,genes.to.plot]
  # build expr data.frame
  expr.df = data.frame("SampleID"=rownames(expr.matrix),
                       "CellType"=tar.tmp$CellType,
                       "PlotColor"=tar.tmp$PlotColor,
                       expr.matrix)
  # create a named vector containing color assignemnts
  # this will be supplied to ggplot to scale colors manually with
  plot.colors = as.character(expr.df$PlotColor)
  names(plot.colors) = as.character(expr.df$CellType)
  # clean up
  rm(expr.dge,expr.matrix)
  rownames(expr.df)=NULL
  # melt the data.frame into long format
  library(reshape2)
  expr.df = melt(expr.df, id.vars=c("SampleID","CellType","PlotColor")) # measure.vars argument left undeclared, all other than id.vars will be used 
  # rename columns
  colnames(expr.df)[colnames(expr.df)=="variable"] = "Gene"
  colnames(expr.df)[colnames(expr.df)=="value"] = "Expr"
  # if the supplied order is empty, the default ggplot2 order will be used
  if (!is.null(group.order)) {expr.df$CellType = factor(expr.df$CellType, levels=group.order)}
  # plot
  library(ggplot2)
  library(extrafont)
  ggplot(expr.df, aes(x=CellType, y=Expr, colour=CellType, group=CellType)) +
    geom_violin(size=1.6, show.legend=FALSE, trim=FALSE, scale="width") + # trim=TRUE violins trimmed to min-max
    geom_point(size=4, shape=19, alpha=0.4, fill="white", show.legend=FALSE, position=position_dodge2(width=0.6, padding=0.1)) +
    scale_color_manual(values=plot.colors) +
    xlab(label=element_blank()) +
    ylab(label=element_blank()) +
    theme_minimal(base_size=44, base_family="Noto Sans Cond") + # base_size=24
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    facet_wrap(~ Gene, nrow=facet.nrow, ncol=facet.ncol, scales=plot.scales)
  plot.name = paste0("expr.plot.",
                     paste(genes.to.plot[1:3], collapse="."),
                     ".pdf")
  ggsave(plot.name,
         device=cairo_pdf,
         width=plot.w.f*facet.ncol, # width=(trunc(length(genes.to.plot)/facet.nrow))*plot.w.f*facet.nrow
         height=plot.h.f*facet.nrow, # (trunc(length(genes.to.plot)/facet.ncol))*plot.h.f*facet.nrow
         limitsize=TRUE)
}

# EGFR, MET only
expr.plot(c("EGFR","MET"),
          plot.h.f=5.8, plot.w.f=5.8, facet.nrow=2, facet.ncol=1,
          group.order=c("hPSC","Ecto","Meso","Endo","TE","CC"))



