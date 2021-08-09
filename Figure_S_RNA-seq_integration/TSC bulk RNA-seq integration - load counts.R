
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# TARGETS ----

library(readxl)
library(data.table)
library(dplyr)

# load targets and select columns
targets.sctl = fread("../TSC_RNA-seq/ANALYSIS/targets.txt", skip=13) %>% select("SampleID","Line","CellType","PlotPCH","PlotColor")
targets.Liu = read_xlsx(path="../Liu 2020 iTSC/PRJNA632917/SRA/PRJNA632917 SraRunTable.xlsx", range="A1:Y24") %>% select("Run","Cell_type")
targets.Io = fread("../Io 2021 nhPSC-TSC/PRJNA605646/SRA/SraRunTable.txt") %>% select("Run","JaroAdded") # JaroAdded is a column I added because the table is a mess, the groups are based on these columns: "source_name","Cell_Line","treatment","Cell_type"



# harmonize columns ----

targets.sctl$PlotPCH = 21
targets.sctl$Batch = 1
targets.sctl$Line = "SCTL"
targets.sctl[targets.sctl$CellType == "TE.p0", ]$CellType = "TE (D10)"
targets.sctl[targets.sctl$CellType == "TSC.early", ]$CellType = "TSC (p7-10)"
targets.sctl[targets.sctl$CellType == "TSC.late", ]$CellType = "TSC (p16-21)"

targets.Liu$SampleID = targets.Liu$Run
targets.Liu$Line = "Liu"
targets.Liu$CellType = targets.Liu$Cell_type
targets.Liu$PlotPCH = 22
targets.Liu$PlotColor = "grey"
targets.Liu[targets.Liu$CellType == "EVTs", ]$PlotColor = "blueviolet"
targets.Liu[targets.Liu$CellType == "human fibroblasts", ]$PlotColor = "bisque3"
targets.Liu[targets.Liu$CellType == "iTSCs", ]$PlotColor = "blue3"
targets.Liu[targets.Liu$CellType == "Naive iPSCs", ]$PlotColor = "darkorange3"
targets.Liu[targets.Liu$CellType == "Primed iPSCs", ]$PlotColor = "cadetblue4"
targets.Liu[targets.Liu$CellType == "STs", ]$PlotColor = "aquamarine3"
targets.Liu$Run = NULL
targets.Liu$Cell_type = NULL
targets.Liu$Batch = 2

targets.Io$SampleID = targets.Io$Run
targets.Io$Line = "Io"
targets.Io$CellType = targets.Io$JaroAdded
targets.Io$PlotPCH = 25
targets.Io$PlotColor = "grey"
targets.Io[targets.Io$CellType == "naive H9", ]$PlotColor = "darkorange3"
targets.Io[targets.Io$CellType == "naive H9-TE D1", ]$PlotColor = "tomato"
targets.Io[targets.Io$CellType == "naive H9-TE D2", ]$PlotColor = "tomato1"
targets.Io[targets.Io$CellType == "naive H9-TE D3", ]$PlotColor = "tomato2"
targets.Io[targets.Io$CellType == "naive H9-CT P3", ]$PlotColor = "darkolivegreen3"
targets.Io[targets.Io$CellType == "naive H9-CT P10", ]$PlotColor = "darkgreen"
targets.Io[targets.Io$CellType == "naive H9-CT P15", ]$PlotColor = "darkgreen"
targets.Io[targets.Io$CellType == "cR-H9-CT P15", ]$PlotColor = "darkgreen"
targets.Io[targets.Io$CellType == "409B2-iPSC-CT P15", ]$PlotColor = "darkgreen"
targets.Io[targets.Io$CellType == "Okae-CT", ]$PlotColor = "blue3"
targets.Io[targets.Io$CellType == "naive H9-CT-ST", ]$PlotColor = "aquamarine3"
targets.Io[targets.Io$CellType == "naive H9-CT-EVT", ]$PlotColor = "blueviolet"
targets.Io[targets.Io$CellType == "primed H9", ]$PlotColor = "cadetblue4"
targets.Io[targets.Io$CellType == "primed H9-BAP D1", ]$PlotColor = "pink"
targets.Io[targets.Io$CellType == "primed H9-BAP D2", ]$PlotColor = "pink1"
targets.Io[targets.Io$CellType == "primed H9-BAP D3", ]$PlotColor = "pink2"
targets.Io[targets.Io$CellType == "placental CT", ]$PlotColor = "blue4"
targets.Io$Run = NULL
targets.Io$JaroAdded = NULL
targets.Io$Batch = 3


# combine ----

targets = rbind(targets.sctl, targets.Liu, targets.Io)

# write into a file
fwrite(targets, file="targets.combined.txt", sep="\t")



# list COUNTS file ----

counts.files.sctl = paste0("../2020-8-TE.exp/COUNTS/STAR-featureCounts/",targets.sctl$SampleID,"-star-counts.txt")
counts.files.Liu = paste0("../PRJNA632917 (Liu 2020 iTSC)/COUNTS/STAR-featureCounts/",targets.Liu$SampleID,"-star-counts.txt")
counts.files.Io = paste0("../Io 2021 nhPSC-TSC/PRJNA605646/COUNTS/STAR-featureCounts/",targets.Io$SampleID,"-star-counts.txt")

# combine in the same order as targets
counts.files = c(counts.files.sctl,counts.files.Liu,counts,counts.files.Io)



# load up and cbind the count data ----

library(data.table)
anno = fread(counts.files[1])[ ,1:6]
# include spot check to see if the order of the genes is the same
# they should since they were all pre-processed with the same pipeline
anno2 = fread(counts.files.Liu[1])[ ,1:6]
anno3 = fread(counts.files.Io[1])[ ,1:6]
all(anno$Geneid == anno2$Geneid)
all(anno$Geneid == anno3$Geneid)
rm(anno2, anno3)
# read all files in
list = lapply(counts.files, function (X) {fread(X)[ ,7]})
# cbind the whole list into a data.frame
# is relatively slower than data.table's rbindlist but there is no "cbindlist" in data.table currently
expr = do.call(cbind,list)
# cbind annotation and expression matrix
raw = cbind(anno,expr)

# clean up
rm(anno,list,expr)
rm(counts.files, counts.files.sctl, counts.files.Liu, counts.files.Io)
rm(targets.sctl, targets.Liu, targets.Io)

# modify column names ----
colnames(raw) = c("GeneID","Chr","Start","End","Strand","Length",targets$SampleID)



# create a DGEList from raw counts ----

library(edgeR)
DGE = DGEList(counts=raw[,7:ncol(raw)],
              genes=raw[,1:6],
              group=paste(targets$CellType,targets$Line,sep="."))



# ANNOTATE the DGE object ----

# to look for attributes available through biomaRt
library(biomaRt)
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart.df = listAttributes(mart)
mart.df[grep(mart.df$name, pattern="entrez"), ]
rm(mart.df)



annotate.DGE = function(ANN=DGE) {
  library(biomaRt)
  mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  annotated.genes = getBM(c("ensembl_gene_id",
                            "external_gene_name",
                            "hgnc_symbol",
                            "description",
                            "chromosome_name",
                            "band",
                            "strand",
                            "start_position",
                            "end_position"),
                          "ensembl_gene_id",ANN$genes$GeneID,mart)
  # remove duplicates - look into differences in duplicated lines!
  genes = annotated.genes[!duplicated(annotated.genes$ensembl_gene_id),]
  # subset ANN
  ANN = ANN[as.character(ANN$genes$GeneID) %in% genes$ensembl_gene_id, ]
  # match IDs
  genes = genes[match(as.character(ANN$genes$GeneID),genes$ensembl_gene_id),] # should handle the fact that fewer genes got annotated, if that's the case
  if (all(as.character(ANN$genes$GeneID) == genes$ensembl_gene_id)) {print("ID match")}
  else {print("ID mismatch!")}
  genes = cbind(ANN$genes,genes,stringsAsFactors=TRUE)
  genes = genes[ ,c("GeneID","ensembl_gene_id","external_gene_name","hgnc_symbol","description",
                    "chromosome_name","band","start_position","end_position","Length","strand")]
  print(colnames(genes))
  ANN$genes = genes
  return(ANN)
}

# use the declared function to annotate DGE
DGE.a = annotate.DGE(DGE)


