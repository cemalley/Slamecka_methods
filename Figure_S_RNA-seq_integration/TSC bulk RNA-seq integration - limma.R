
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# LIMMA ----


# SUBSET ----

# if no subsetting is needed
tar.tmp = targets
DGE.a.subset = DGE.a



# batch correction with ComBat-seq ----

library(sva)
adjusted.counts = ComBat_seq(DGE.a.subset$counts, batch=tar.tmp$Batch, group=NULL)

DGE.a.subset$counts = adjusted.counts
rm(adjusted.counts)



# FILTERING ----

filter.DGE = function(DGE.object) {
require(limma)
require(edgeR)
# design matrix for filtering
factor = DGE.object$samples$group
design = model.matrix(~0+factor)
colnames(design) = levels(factor)
# filtering, new function
keep = filterByExpr(DGE.object, design=design)
print(table(keep))
DGE.object = DGE.object[keep, ]
return(DGE.object)
}

DGE.f = filter.DGE(DGE.a.subset)



# TMM normalization ----

DGE.fn = calcNormFactors(DGE.f)



# create a file name base ----

fn = list()
fn$base = "star-featureCounts.SCTL-Liu-Io"



# MDS plot function ----

plot.MDS.from.DGE = function(DGE.object = DGE.fn,
                             plot.targets = tar.tmp, # must have columns: Line, CellType, PlotPCH, PlotColor
                             line.legend = TRUE,
                             celltype.legend = TRUE,
                             celltype.order = NULL, # indices of vector
                             line.legend.pos = c(-2,2), # position of cell line legend, can be c("topleft",NULL)
                             celltype.legend.pos = c(2,2), # position of cell type legend, can be c("topright",NULL)
                             plot.w = 8,
                             plot.h = 8.2,
                             plot.filename = paste0("MDS.",fn$base,".symbols.pdf")
                             ) {
        # if font change is needed, uncomment par() in line after opening up the device to change the font
        # library(extrafont)
        # Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.52/bin/gswin64c.exe") # part of the initial setup, not needed every time
        
        # set order of cell types and cell lines to appear on the legend
        cols.df = unique(plot.targets[ ,c("CellType", "PlotColor")])
        # reorder cell types if celltype.order vector provided
        if (!is.null(celltype.order)) {cols.df = cols.df[celltype.order, ]}
        # make data.frame of line and point shape association
        points.df = unique(plot.targets[ ,c("Line", "PlotPCH")])
        
        require(extrafont)
        require(Cairo) # will embed fonts automatically, no need to use embed_fonts() from "extrafont" package
        
        CairoPDF(file=plot.filename, width=plot.w, height=plot.h) # family="Noto Sans Cond" - as an argument doesn't seem to work here, setting par in the next line does
        par(xpd=TRUE, family="Noto Sans Cond") # cex=1.2, lwd=2 - to control text and point size; and line width
        plotMDS(DGE.fn, pch=plot.targets$PlotPCH, bg=plot.targets$PlotColor, col="white", cex=2, lwd=0.8)
        if (line.legend) {
                legend(x=line.legend.pos[1], y=line.legend.pos[2], pch=points.df$PlotPCH, legend=points.df$Line,
                       inset=0.01, pt.cex=1.2, cex=1.2, text.width=1.2, y.intersp = 0.8, x.intersp=0.8#,
                       # bty="n" # comment for box  around the legend - default value is "o"
                )}
        if (celltype.legend) {
                legend(x=celltype.legend.pos[1], y=celltype.legend.pos[2], fill=cols.df$PlotColor, legend=cols.df$CellType,
                       border="white", inset=0.01, pt.cex=1.2, cex=1.2, text.width=2, y.intersp = 0.8, x.intersp=0.2#,
                       # bty="n" # comment for box  around the legend - default value is "o"
                )}
        dev.off()
}


# draw the MDS plots ----

# SCTL-Liu-Io integration
plot.MDS.from.DGE(celltype.order = c(1,3,2,4, # SCTL
                                     6,8,9,7,5,10, # Liu
                                     11,23:26,12:20,27,21,22), # Io
                  line.legend.pos = c("topleft",NULL), celltype.legend.pos = c("topright",NULL),
                  plot.h = 10, plot.w = 10,
                  plot.filename = paste0("MDS.",fn$base,".symbols.filtered.+legend.pdf"))
# SCTL-Liu-Io integration, no legend
plot.MDS.from.DGE(celltype.order = c(1,3,2,4,
                                     6,8,9,7,5,10,
                                     11,23:26,12:20,27,21,22),
                  line.legend = FALSE, celltype.legend = FALSE,
                  plot.w = 10, plot.h = 10,
                  plot.filename = paste0("MDS.",fn$base,".symbols.filtered.-legend.pdf"))


