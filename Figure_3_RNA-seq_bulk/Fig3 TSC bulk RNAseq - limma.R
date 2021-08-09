
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# LIMMA ----

# if no subsetting is needed
tar.tmp = targets
DGE.f.subset = DGE.f



# create a file name base ----

fn = list()
fn$base = "star-featureCounts"



# TMM normalization ----

DGE.fn = calcNormFactors(DGE.f.subset)
rm(DGE.f.subset)
plotMDS(DGE.fn)



# MDS plot ----

# set order of cell types and cell lines to appear on the legend
cols.df = unique(tar.tmp[ ,c("CellType", "PlotColor")])
cols.df = cols.df[c(1,3,2,4), ]
points.df = unique(tar.tmp[ ,c("Line", "PlotPCH")])

library(extrafont)
library(Cairo) # will embed fonts automatically, no need to use embed_fonts() from "extrafont" package
CairoPDF(file=paste0("MDS.",fn$base,".symbols.pdf"), width=4, height=4.4) # family="Noto Sans Cond" - as an argument doesn't seem to work here, setting par in the next line does
par(family="Noto Sans Cond") # cex=1.2, lwd=2 - to control text and point size; and line width
plotMDS(DGE.fn, pch=tar.tmp$PlotPCH, bg=tar.tmp$PlotColor, col="white", cex=2, lwd=0.8)
legend(x="left", pch=points.df$PlotPCH, legend=points.df$Line,
       inset=0.01, pt.cex=1.2, cex=1.2, text.width=2.6, y.intersp = 0.8, x.intersp=0.8#,
       # bty="n" # comment for box  around the legend - default value is "o"
       )
legend(x=-0.38, y=2.8, fill=cols.df$PlotColor, legend=cols.df$CellType, # x=-0.62, y=3.4 with hypoxia group
       border="white", inset=0.01, pt.cex=1.2, cex=1.2, text.width=3.2, y.intersp = 0.8, x.intersp=0.2#,
       # bty="n" # comment for box  around the legend - default value is "o"
       )
dev.off()

# text version to show all labels
pdf(file=paste0("MDS.",fn$base,".labels.pdf"), width=5, height=5)
plotMDS(DGE.fn, cex=0.2)
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

# lmFit ----

fit = lmFit(DGE.fnv, design)

# choose contrasts ----

fn$contrast = paste0(fn$base,".TSC-TEp0.late-early.TEp0-hPSC.TSC-hPSC")
contrast.matrix = makeContrasts(TSC.vs.TEp0 = (TSC.early.WA07+TSC.early.WA09+TSC.late.WA09+TSC.early.JHU191i+TSC.late.JHU191i+TSC.early.JHU198i+TSC.late.JHU198i+TSC.early.MCW032i)/8 - (TE.p0.WA07+TE.p0.WA09+TE.p0.JHU191i+TE.p0.JHU198i+TE.p0.MCW032i)/5, # TSC early+late (no hypoxia) vs. TE.p0
                                late.vs.early = (TSC.late.WA09+TSC.late.JHU191i+TSC.late.JHU198i)/3 - (TSC.early.WA09+TSC.early.JHU191i+TSC.early.JHU198i)/3, # late vs. early in 3 cell lines
                                TEp0.vs.hPSC = (TE.p0.WA07+TE.p0.WA09+TE.p0.JHU191i+TE.p0.JHU198i+TE.p0.MCW032i)/5 - (hPSC.WA07+hPSC.WA09+hPSC.JHU191i+hPSC.JHU198i+hPSC.MCW032i)/5, # TE.p0 vs. hPSC in all 5 cell lines
                                TSC.vs.hPSC = (TSC.early.WA07+TSC.early.WA09+TSC.late.WA09+TSC.early.JHU191i+TSC.late.JHU191i+TSC.early.JHU198i+TSC.late.JHU198i+TSC.early.MCW032i)/8 - (hPSC.WA07+hPSC.WA09+hPSC.JHU191i+hPSC.JHU198i+hPSC.MCW032i)/5, # TSC vs. hPSC
                                levels=design)

# fit and eBayes smoothing ----

fit = contrasts.fit(fit, contrast.matrix)
fit = eBayes(fit)



# FILTERING, ORDERING based on P-value ----

# re-order all genes based on evidence of differential expression
top = topTable(fit, adjust='fdr', number=Inf)
# write results to file
write.table(top,file=paste0("DGE.",fn$contrast,".txt"), sep="\t", quote=FALSE, col.names=NA)



# MISSING GENE SYMBOLS ----

# fill missing HGNC gene symbols with ensembl gene IDs
# skip if genes without a symbol are not to be included
top[top$hgnc_symbol == "", ]$hgnc_symbol = top[top$hgnc_symbol == "", ]$GeneID
DGE.fnv$genes[DGE.fnv$genes$hgnc_symbol == "", ]$hgnc_symbol = DGE.fnv$genes[DGE.fnv$genes$hgnc_symbol == "", ]$GeneID
# check how many HGNC symbols now have Ensembl gene IDs
length(grep(DGE.fnv$genes$hgnc_symbol, pattern="ENSG", value=TRUE))



# choose CUT-OFFs for MARKER list ----

CO = list("logFC"=list("global"=8.00,
                       "TSC"=3.00,
                       "TE.p0"=6.00,
                       "TSC.late"=2.00,
                       "hypoxia"=3.00,
                       "TSC.vs.hPSC"=8.00),
          "adjPval"=list("global"=1e-4, # only global used at the moment
                         "TSC"=NULL,
                         "TEp0"=NULL,
                         "TSC.late"=NULL,
                         "hypoxia"=NULL,
                         "TSC.vs.hPSC"=NULL))

# filter based on FDR-adjusted P value first
top.pval = top[top$adj.P.Val < CO$adjPval$global, ]
markers = list("TSC"=top.pval[top.pval$TSC.vs.TEp0 >= CO$logFC$TSC & !(top.pval$hgnc_symbol == ""), ]$hgnc_symbol,
               "TE.p0"=top.pval[top.pval$TEp0.vs.hPSC >= CO$logFC$TE.p0 & !(top.pval$hgnc_symbol == ""), ]$hgnc_symbol,
               "TSC.late"=top.pval[top.pval$late.vs.early >= CO$logFC$TSC.late & !(top.pval$hgnc_symbol == ""), ]$hgnc_symbol,
               "hypoxia"=top.pval[top.pval$hypoxia.vs.TSC >= CO$logFC$hypoxia & !(top.pval$hgnc_symbol == ""), ]$hgnc_symbol,
               "TSC.vs.hPSC"=top.pval[top.pval$TSC.vs.hPSC >= CO$logFC$TSC.vs.hPSC & !(top.pval$hgnc_symbol == ""), ]$hgnc_symbol
               )
# get number of genes in each marker group
sapply(markers, length)
rm(top.pval)

# write cutoffs into file
writeLines(capture.output(CO), con=paste0("DGE.",fn$base,".CO.file.txt"))

# get a simple top selected genes subset
top.sel = top.pval[top.pval$TE.WA09...hPSC.WA09 > CO$logFC$global , ]



# add filter of the top table to the "fn" filename list ----

fn$selected = paste0(fn$base,".cutoffs.in.CO.file")



# GENES OF INTEREST ----

basepath = "./"

library(tidyverse)
library(readxl)
library(data.table)

GOI = list()
# look at placental genes
GOI$placenta = read.delim(paste0(basepath,"placental.genes.txt"),header=TRUE,sep="\t",comment.char="#")
# manually modify the placental genes according to the desired way of plotting the heatmap
GOI$placenta = GOI$placenta[GOI$placenta$Include.2019.6.TE.53==TRUE, ]
GOI$placenta$Origin = as.character(GOI$placenta$Origin) # change factor to character, matchup won't work otherwise
# Protein Atlas - placenta enriched/specific genes
GOI$placenta.PA.enriched = read.delim(paste0(basepath,"DATA ANALYSIS/NGS/GOI/[2020-4-6] Protein Atlas - placenta-specific genes (91).tsv"),header=TRUE,sep="\t",comment.char="#")
# Zhou et al., 2019: Reconstituting the transcriptome and DNA methylome landscapes of human implantation
GOI$scTrMet = read_xlsx(paste0(basepath,"ARTICLES/PSC Differentiation/Trophectoderm/[2019, Nature] Reconstituting the transcriptome and DNA methylome landscapes of human implantation - SUPP TABLES/Supplementary Table 3 Lineage_Markers.xlsx"))
GOI$scTrMet = GOI$scTrMet[ ,c(1,7,2:6)] # reorder columns, function make.feature.groups expects symbols in column 1 and group in column 2
GOI$scTrMet = as.data.frame(GOI$scTrMet)
rm(basepath)



# HEATMAP ----

heatmap = list()
# choose overall font size for heatmap labels
heatmap$fontsize = 22
heatmap$fontfamily = "Noto Sans Cond"
heatmap$fontfamily = ""

# choose features for heatmap ----

# select features manually for mark annotation
heatmap$features = c(markers$TSC,
                     markers$TE.p0)
markers$manual = list()
markers$manual$TSC = c("ENPEP","FOLR1","MKI67","SIGLEC6","DUSP6","XAGE3","PSG8","MAGEA8",
                       "PSG7","PSG11","PSG1","MAGEB16","PSG2","PSG5","MAGEA10","PSG9","PSG6",
                       "PSG3","HAVCR1","MAGEB2","PSG4")
markers$manual$TE.p0 = c("STS","MBNL3","VGLL1","EPAS1","TFAP2B","PRLR","PPARG","CRH","XAGE2","H19",
                         "PAPPA2","GATA3","IGFBP3","PAGE4","ELF5","TFAP2A","PSG2","TP63","HAND1","SVEP1",
                         "PSG4","PHLDA2","CDX2","PGF","PLAC4","INSL4",
                         "GCM1","CGA","DIO3")
# GOI - genes of interest
heatmap$features = as.character(GOI$placenta$GeneSymbol)
heatmap$features = as.character(GOI$placenta.PA.enriched$Gene)
heatmap$features = as.character(GOI$scTrMet$gene)
heatmap$features = as.character(GOI$scTrMet[GOI$scTrMet$cluster == "TE", ]$gene)



# SUBSET for heatmap ----

# subsetting based on either Ensembl gene IDs and symbols
# choose hgnc_symbol if no genes with missing symbols are plotted on heatmap (unless missing symbols were replaced by Ensembl Gene IDs above)
heatmap$DGE = DGE.fnv[DGE.fnv$genes$external_gene_name %in% heatmap$features | DGE.fnv$genes$ensembl_gene_id %in% heatmap$features, ]
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
# limma voom values
# include gene symbols
heatmap$matrix = heatmap$DGE$E + abs(min(DGE.fnv$E)) # adding the absolute value of the (global) minimum to avoid negative values
rownames(heatmap$matrix) = heatmap$DGE$genes$external_gene_name
colnames(heatmap$matrix) = colnames(heatmap$DGE)



# heatmap ANNOTATION - side colors ----
# can be column and row side colors
# get levels of the group factor and manually assign colors
levels(heatmap$DGE$targets$group)
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



# heatmap ANNOTATION - MARK ----
# highlight selected genes on a crowded heatmap
# anno_mark() needs indices of the positions and the corresponding labels
# choose genes to highlight
selected.genes = c(markers$manual$TSC, markers$manual$TE.p0)
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
# delete a slot from the markers list in order to exclude those genes (and plot them separately)
heatmap$feature.groups = make.heatmap.feature.groups(markers[1:2], 2, heatmap$DGE,
                                                     group.order=c("TE.p0","TSC"))



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
heatmap$col.dendro = get.col.dendro(col.to.use="TEp0.vs.hPSC", co.logFC=4, co.P=1e-04) # TSC.vs.TEp0"
heatmap$col.dendro = get.col.dendro(col.to.use="TSC.vs.hPSC", co.logFC=4, co.P=1e-04) # TSC.vs.hPSC"
# re-order samples on the dendrogram
# rotate
labels(heatmap$col.dendro)
heatmap$col.dendro = rotate(heatmap$col.dendro, order=c(53:74,1:22,23:52))
heatmap$col.dendro = rotate(heatmap$col.dendro, order=c(53:74,31:52,1:30))
# save and load
saveRDS(heatmap$col.dendro, file="col.dendro.noHypoxia.rds")
col.dendro = readRDS("col.dendro.noHypoxia.rds")
heatmap$col.dendro = col.dendro
rm(col.dendro)
# no hypoxia and no TE.p0 groups



# DRAW heatmaps and save PDFs ----
# color key legend parameters are modified within the Heatmap function
# annotation legend parameters are modified above in the HeatmapAnnotation function
heatmap$pdf.w = 16 
CairoPDF(file=paste0("heatmap.",fn$selected,".",nrow(heatmap$matrix),"features.pdf"), width=heatmap$pdf.w, height=24) # for top DEG - TSC markers
heatmap$heatmap = Heatmap(heatmap$matrix,
                          col=heatmap$colors,
                          top_annotation=heatmap$col.annotation,
                          right_annotation=heatmap$mark.annotation, show_row_names=FALSE,
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

# draw expression plot ----

expr.plot = function(genes.to.plot,
                     plot.h.f=6, # plot height factor
                     plot.w.f=4.6, # plot width factor
                     facet.nrow=2, # number of facet rows
                     facet.ncol=5, # number of facet columns
                     counts=FALSE, log.tr=FALSE, # use normalized counts, either intact or log2-transformed, instead of transformed values
                     group.order=NULL, # order of groups on plot
                     plot.per.line=FALSE # split x axis to plot expression values of individual cell lines
                     ) {
  # subset DGE based on genes to plot
  # use DGE.fnv for transformed values
  # use DGE.fn for normalized counts - with argument counts=TRUE
  if (counts) {
    expr.dge = DGE.fn[DGE.fn$genes$external_gene_name %in% genes.to.plot | DGE.fn$genes$ensembl_gene_id %in% genes.to.plot, tar.tmp$SampleID]
    if (log.tr) {
      expr.matrix = log2(expr.dge$counts + 1)
      plot.scales="fixed"
      } else {
        expr.matrix = expr.dge$counts
        plot.scales="free"
        }
    } else {
    expr.dge = DGE.fnv[DGE.fnv$genes$external_gene_name %in% genes.to.plot | DGE.fnv$genes$ensembl_gene_id %in% genes.to.plot, tar.tmp$SampleID]
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
  if (plot.per.line) {expr.df = data.frame("SampleID"=rownames(expr.matrix),
                                           "CellType"=paste(tar.tmp$CellType, tar.tmp$Line, sep="-"), # combine cell type and line to plot expression values split per line
                                           "PlotColor"=tar.tmp$PlotColor,
                                           expr.matrix)}
  else {expr.df = data.frame("SampleID"=rownames(expr.matrix),
                       "CellType"=tar.tmp$CellType,
                       "PlotColor"=tar.tmp$PlotColor,
                       expr.matrix)}
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
    geom_violin(size=1.6, show.legend=FALSE, trim=FALSE, scale="width") + # trim=TRUE trimmed violins to min-max
    geom_point(size=4, shape=19, alpha=0.4, fill="white", show.legend=FALSE, position=position_dodge2(width=0.6, padding=0.1)) +
    scale_color_manual(values=plot.colors) +
    xlab(label=element_blank()) +
    ylab(label=element_blank()) +
    theme_minimal(base_size=44, base_family="Noto Sans Cond") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    facet_wrap(~ Gene, nrow=facet.nrow, ncol=facet.ncol, scales=plot.scales)
  plot.name = paste0("expr.plot.",
                     paste(genes.to.plot[1:3], collapse="."),
                     ".pdf")
  ggsave(plot.name,
         device=cairo_pdf,
         width=plot.w.f*facet.ncol,
         height=plot.h.f*facet.nrow,
         limitsize=TRUE)
}


# modify targets
tar.tmp[tar.tmp$CellType == "TE.p0", ]$CellType = "TE (D10)"
tar.tmp[tar.tmp$CellType == "TSC.early", ]$CellType = "TSC (p7-10)"
tar.tmp[tar.tmp$CellType == "TSC.late", ]$CellType = "TSC (p16-21)"

expr.plot(c("TP63","CDX2","IGFBP3","GATA3","VGLL1", # ELF5 plotted separately
            "TFAP2A","TFAP2B","TFAP2C",
            "PAGE4","PEG10","PARP1", # villous cytotrophoblast markers (Suryawanshi et al., 2018)
            "NR2F2","LRP2", # +PEG10 hTSC markers (Castel et al., 2020)
            "XAGE3","HAVCR1","SIGLEC6","MAGEA10",
            "PSG2","HLA-A","HLA-B"),
          plot.h.f=5.6, plot.w.f=4.4, facet.nrow=4, facet.ncol=5,
          group.order=c("hPSC","TE (D10)","TSC (p7-10)","TSC (p16-21)"))

# expression of ELF5
expr.plot(c("ELF5","GATA3"),
          plot.h.f=7, plot.w.f=2.4, facet.nrow=1, facet.ncol=4,
          group.order=c("hPSC","TE (D10)","TSC (p7-10)","TSC (p16-21)"))


