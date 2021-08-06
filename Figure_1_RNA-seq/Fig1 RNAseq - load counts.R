
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# load TARGETS ----

# define base path ----
basepath = c("../")

# load targets in ----
targets = read.delim("targets.txt",
                     header=TRUE,sep="\t",comment.char="#",
                     stringsAsFactors=FALSE)

# check if I assigned ahspes properly
unique(targets[ ,c("Factor", "PlotPCH")])



# load COUNTS files ----

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
colnames(raw) = c("GeneID","Chr","Start","End","Strand","Length",targets$SampleName)

# create a DGEList from raw counts ----
library(edgeR)
DGE = DGEList(counts=raw[,7:ncol(raw)],
              genes=raw[,1:6],
              group=factor(targets$Factor, levels=c("D0","D3","D7","D10")))



# ANNOTATE the DGE object ----

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
  # remove duplicates
  genes = annotated.genes[!duplicated(annotated.genes$ensembl_gene_id),]
  # subset ANN
  ANN = ANN[as.character(ANN$genes$GeneID) %in% genes$ensembl_gene_id, ]
  # match IDs
  genes = genes[match(as.character(ANN$genes$GeneID),genes$ensembl_gene_id),]
  if (all(as.character(ANN$genes$GeneID) == genes$ensembl_gene_id)) {print("ID match")}
  else {print("ID mismatch!")}
  # [1] TRUE
  genes = cbind(ANN$genes,genes,stringsAsFactors=TRUE)
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
# filtering
keep = filterByExpr(DGE.a, design=design)
table(keep)
DGE.f = DGE.a[keep, ]
# clean up
rm(factor,design,keep)
