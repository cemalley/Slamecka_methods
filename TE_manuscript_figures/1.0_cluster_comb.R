#' Test cluster merging for TE analysis for Jaro
#' Compile basic stats
#' 
#' Compile adjusted figures as CM discussed on call

# req'd pkgs
x <- c("Seurat", "dplyr", "data.table", "ggplot2", "gridExtra")
sapply(x, library, character.only = TRUE)

# load seurat obj
load("./2019-6-TE.sc.RData")
seur <- tenx.seurat
# update cluster resolution to 0.07 as decided by
# CM & BE
seur <- Seurat::FindClusters(seur, resolution = 0.07)

# Create DGE plot w/o merging cluster 1 & 3
Idents(seur) <- "seurat_clusters"
# find all DGE genes for cluster x vs rest for all clusters
de_first <- FindAllMarkers(seur)
# filter to retain genes > logfc 0.25 and below adj-pval <= 0.001
adj_first <- de_first[de_first$avg_logFC > 0.25 & de_first$p_val_adj <= 0.001, ]
# keep unique genes per cluster
adj_first <- adj_first[!(adj_first$gene %in% unique(adj_first$gene[duplicated(adj_first$gene)])), ]
# summarize count of unique, positive DGE genes
sum_first <- adj_first %>% group_by(cluster) %>% summarise(n = n())
sum_first$cluster <- factor(sum_first$cluster, levels = c("0", "1", "2", "3", "4", "5", "6"))

# plot count of unique, positive dge genes by cluster
a <- ggplot(sum_first, aes(x = cluster, y = n)) + geom_bar(stat = "identity", fill = "grey", color = "black") + 
  geom_text(size = 3, label = sum_first$n, vjust = -1) +
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     panel.border = element_blank(),
                     axis.line.y = element_line(colour = "black"),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_text(size = 14),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title.y = element_text(size = 14),
                     plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  xlab("Cluster") + ylab("# Unique, positive DE genes") +
  scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + ylim(c(0, 650)) +
  ggtitle(paste("Resolution: 0.07\nTotal: ", sum(sum_first$n), sep = ""))

# merge clusters 1 & 3
clus <- c(1, 3)
# extract cluster data from seurat_clusters column in meta.data
res <- seur@meta.data$seurat_clusters
res <- as.numeric(res)
# subtract 1 from cluster numbers to line up w/seurat
res <- res - 1
# combine clusters 1 & 3 into 100
res <- ifelse(res %in% clus, 100, res)
res <- as.character(res)
# replace 100 w/comb_1_3
res <- ifelse(res == "100", "comb_1_3", res)
res <- as.factor(res) 

# add updated cluster numbers to seurat object
# in clus_comb
seur$clus_comb <- res
Idents(seur) <- "clus_comb"

# repeat steps above for updated cluster #'s
de <- FindAllMarkers(seur)
adj <- de[de$avg_logFC > 0.25 & de$p_val_adj <= 0.001, ]
adj <- adj[!(adj$gene %in% unique(adj$gene[duplicated(adj$gene)])), ]
sum <- adj %>% group_by(cluster) %>% summarise(n = n())
sum$cluster <- factor(sum$cluster, levels = c("0", "comb_1_3", "2", "4", "5", "6"))

b <- ggplot(sum, aes(x = cluster, y = n)) + geom_bar(stat = "identity", fill = "grey", color = "black") + 
  geom_text(size = 3, label = sum$n, vjust = -1) +
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     panel.border = element_blank(),
                     axis.line.y = element_line(colour = "black"),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_text(size = 14),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title.y = element_text(size = 14),
                     plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  xlab("Cluster") + ylab("# Unique, positive DE genes") +
  scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + ylim(c(0, 650)) +
  ggtitle(paste("Resolution: 0.07\nTotal: ", sum(sum$n), sep = ""))

p <- grid.arrange(a, b, nrow = 1)
ggsave("./combine_cluster1_3.png", plot = p)

# re-label cluster numbering based on merge of 
# cluster 1 & 3 so it's contiguous
clus <- as.character(seur@meta.data$clus_comb)
clus <- ifelse(clus == "comb_1_3", "4",
               ifelse(clus == "2", "3",
                  ifelse(clus == "4", "2",
                      ifelse(clus == "5", "1", 
                             ifelse(clus == "6", "5", clus)))))
clus <- as.factor(clus)
seur$upd_res <- clus
Idents(seur) <- seur$upd_res
seur@meta.data$day <- ifelse(seur@meta.data$upd_res %in% c("0", "1"), "D0",
                       ifelse(seur@meta.data$upd_res %in% c("2", "3"), "D3",
                        ifelse(seur@meta.data$upd_res == "4", "D6/D9",
                         ifelse(seur@meta.data$upd_res == "5", "diff", "NA"))))
seur@meta.data$day <- as.factor(seur@meta.data$day)
saveRDS(seur, file = "./upd_tenx.seurat.RDS")
