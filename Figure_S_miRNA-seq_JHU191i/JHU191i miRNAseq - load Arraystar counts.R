
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# load targets in ----

targets = read.delim("targets.txt",
                     header=TRUE,sep="\t",comment.char="#",
                     stringsAsFactors=FALSE)



# choose COUNTS file ----

counts.file = "Expression_miRNA.xlsx"

library(readxl)
raw = read_xlsx(path=counts.file, range="A22:AB1069")

writeLines(colnames(raw)[19:37], con="samples.txt") # for building of targets file



# create a DGEList from raw counts ----

library(edgeR)
DGE = DGEList(counts=raw[,19:ncol(raw)],
              genes=raw[,1:18],
              group=paste(targets$CellType,targets$Line,sep="."))



# FILTERING ----

# design matrix for filtering
factor = DGE$samples$group
design = model.matrix(~0+factor)
colnames(design) = levels(factor)
# filtering, new function
keep = filterByExpr(DGE, design=design)
table(keep)
DGE.f = DGE[keep, ]
# clean up
rm(factor,design,keep)
