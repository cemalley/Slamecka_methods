# save Jaro's single cell data for partek for pcas there------
is021 <- mergedSeuratProcessed
rm(mergedSeuratProcessed)

is021@active.ident
is021@meta.data$sample <- is021@meta.data$orig.ident
is021@meta.data$sample <- gsub('H9TE_D3_i527','H9TE_D3',is021@meta.data$sample)
is021@meta.data$sample <- gsub('H9TE_D6_i524','H9TE_D6',is021@meta.data$sample)
is021@meta.data$sample <- gsub('H9TE_D9_i521','H9TE_D9',is021@meta.data$sample)

Idents(is021) <- 'sample'
is021@active.ident

metadata <- as.data.frame(is021@meta.data)
metadata$barcode <- row.names(metadata)
metadata <- as.data.table(metadata)

metadata

unique(metadata$sample) # "H9TE_D0" "H9TE_D3" "H9TE_D6" "H9TE_D9" "P17"     "P23"

rawcounts <- as.data.frame(is021@assays$RNA@counts)


H9TE_D0 <- subset(rawcounts, select=c(metadata[sample=='H9TE_D0',barcode]))
H9TE_D3 <- subset(rawcounts, select=c(metadata[sample=='H9TE_D3',barcode]))
H9TE_D6 <- subset(rawcounts, select=c(metadata[sample=='H9TE_D6',barcode]))
H9TE_D9 <- subset(rawcounts, select=c(metadata[sample=='H9TE_D9',barcode]))
P17 <- subset(rawcounts, select=c(metadata[sample=='P17',barcode]))
P23 <- subset(rawcounts, select=c(metadata[sample=='P23',barcode]))

H9TE_D0[1:5,1:5]

write.csv(H9TE_D0, file='/Volumes/ncatssctl/NGS_related/Chromium/IS021/Partek/H9TE_D0.csv')
write.csv(H9TE_D3, file='/Volumes/ncatssctl/NGS_related/Chromium/IS021/Partek/H9TE_D3.csv')
write.csv(H9TE_D6, file='/Volumes/ncatssctl/NGS_related/Chromium/IS021/Partek/H9TE_D6.csv')
write.csv(H9TE_D9, file='/Volumes/ncatssctl/NGS_related/Chromium/IS021/Partek/H9TE_D9.csv')
write.csv(P17, file='/Volumes/ncatssctl/NGS_related/Chromium/IS021/Partek/P17.csv')
write.csv(P23, file='/Volumes/ncatssctl/NGS_related/Chromium/IS021/Partek/P23.csv')


# reformat for GSEA--------
load("/Volumes/ncatssctl/NGS_related/Chromium/IS021/2019-6-TE.sc.RData")
Idents(tenx.seurat) <- 'orig.ident'

counts <- as.data.frame(as.matrix(tenx.seurat@assays$RNA@counts))
counts[1:3,1:3]
counts$NAME <- row.names(counts)
counts$DESCRIPTION <- 'na'
counts <- as.data.table(counts)
counts <- counts[,c(11055,11056, 1:11054)]

genes.filter.2 <- unique(c(genes.filter$V1,
                         counts$NAME[grep('^AL[0:9]',counts$NAME)],
                         counts$NAME[grep('^LINC[0:9]',counts$NAME)]))

counts <- counts[NAME %nin% c(genes.filter.2),]

fwrite(counts, '/Volumes/ncatssctl/NGS_related/Chromium/IS021/GSEA/IS021_raw_counts.gct', col.names=T, sep='\t')

phenotypes <- as.data.frame(tenx.seurat@meta.data)
unique(phenotypes$orig.ident)
phenotypes$Day <- phenotypes$orig.ident
phenotypes$Day <- gsub('H9TE_D0',0,phenotypes$Day)
phenotypes$Day <- gsub('H9TE_D3_i527',3,phenotypes$Day)
phenotypes$Day <- gsub('H9TE_D6_i524',6,phenotypes$Day)
phenotypes$Day <- gsub('H9TE_D9_i521',9,phenotypes$Day)
unique(phenotypes$Day)

cat(phenotypes$Day, sep='\t')

#
