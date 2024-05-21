
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# TARGETS ----

# define base path ----
basepath = c("../")

# load targets in ----
targets = read.delim(paste0(basepath,"ANALYSIS/targets.txt"),
                     header=TRUE,sep="\t",comment.char="#",
                     stringsAsFactors=FALSE)


# rename cell lines for consistency ----
# "ESC">>"hPSC" | "H9">>"WA09"
targets[targets$CellType == "ESC", ]$CellType = "hPSC"
targets[targets$Line == "H9", ]$Line = "WA09"
targets[targets$Line == "H14", ]$Line = "WA14"



# choose COUNTS file ----

# single-file counts
counts.files = paste0(basepath,"COUNTS/STAR-featureCounts/",targets$SampleID,"-star-counts.txt")

# load up and cbind the count data ----
library(data.table)
anno = fread(counts.files[1])[ ,1:6]
# read all files in
list = lapply(counts.files, function (X) {fread(X)[ ,7]})
# cbind the whole list into a data.frame
expr = do.call(cbind,list)
# cbind annotation and expression matrix
raw = cbind(anno,expr)

# clean up
rm(anno,list,expr,basepath,counts.files)

# modify column names ----
colnames(raw) = c("GeneID","Chr","Start","End","Strand","Length",targets$SampleName.LIMS)

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
  genes = genes[match(as.character(ANN$genes$GeneID),genes$ensembl_gene_id),]
  if (all(as.character(ANN$genes$GeneID) == genes$ensembl_gene_id)) {print("ID match")}
  else {print("ID mismatch!")}
  # [1] TRUE
  genes = cbind(ANN$genes,genes,stringsAsFactors=TRUE)
  # genes = genes[ ,c(1,14,2:13)]
  genes = genes[ ,c("GeneID","ensembl_gene_id","external_gene_name","hgnc_symbol","description",
                    "chromosome_name","band","start_position","end_position","Length","strand")]
  print(colnames(genes))
  ANN$genes = genes
  return(ANN)
  }

# use the declared function to annotate DGE
DGE.a = annotate.DGE(DGE)



# FILTERING ----

# design matrix for filtering
factor = DGE.a$samples$group
design = model.matrix(~0+factor)
colnames(design) = levels(factor)
# filtering, new function
keep = filterByExpr(DGE.a, design=design)
table(keep)
DGE.f = DGE.a[keep, ]
# clean up
rm(factor,design,keep)
