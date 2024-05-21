# Ben Ernest
# R functions

# Reads in CSV and returns Seurat object
MetaData2Seurat = function(file) {
  #file = "/root/ec2_volume/ernest_09102019/seurat_testing_09102019/IS021_merged_metadata.csv"
  ext = strsplit(x = file, split = "\\.")[[1]][-1]
  if(! ext %in% c("csv", "txt")) {
    stop("Input file must be csv or tab-delimited txt.")
  }
  if(ext == "csv") sep = ","
  if(ext == "txt") sep = "\t"
  df = read.table(file = file, header = T, sep = sep, row.names = 1, 
                  as.is = T, check.names = F)
  return(df)
}

# Performs standard processing of scRNA-seq data using Seurat methods.
# The function checks whether each step has been performed and performs 
# steps that haven't been performed.
SeuratProcess = function(object, uniqueFeatureCountRange = c(200, 2500),
                         mtPercent = 5, allGenes = T, res = 0.8, isTPM = F,
                         prenormalized = F, assay = "RNA", justRed = F,
                         reprocess = F) {
  # load(file = "/root/ec2_volume/ernest_09102019/TE_vs_placenta_09102019/placenta.reorder.Rda")
  # load(file = "/root/ec2_volume/ernest_09102019/TE_vs_placenta_09102019/IS021.seurat.RData")
  # 
  # placentaSeurat = subset(x = placenta.reorder, subset = orig.ident %in% c("P17", "P23"))
  # is021Seurat = subset(x = IS021.seurat, subset = orig.ident == "H9TE_D9_i521")
  # mergedSeurat = merge(x = is021Seurat, y = placentaSeurat)
  # object = mergedSeurat
  # object = readRDS("TE_vs_placenta_09102019/seurat_processed_cell-labels_10082019.rds")
  # object = zhouEmbSeurat
  # object = seuratIntegrated
  # object = mergedSeurat
  # uniqueFeatureCountRange = c(200, 2500)
  # mtPercent = 5
  # allGenes = T
  # res = 0.8
  # isTPM = F
  # prenormalized = T
  # assay = "RNA"
  
  # Check whether percentage of reads mapping to mitochondrial genome has been 
  # calculated
  if(!justRed) {
    if(is.null(object@meta.data$percent.mt) | reprocess) {
      print("Calculating mitochondrial percentages...")
      object[["percent.mt"]] = PercentageFeatureSet(object = object, pattern = "^MT-", assay = assay)
    }
  }
  
  # cellsToKeep = WhichCells(object = object, nFeature_RNA > uniqueFeatureCountRange[1] &
  #                            nFeature_RNA < uniqueFeatureCountRange[2] &
  #                            percent.mt < mtPercent)
  
  # uniqueFeatureCountRange = uniqueFeatureCountRange
  # mtPercent = mtPercent
  # object = subset(x =object, subset = nFeature_RNA > uniqueFeatureCountRange[1] &
  #                   nFeature_RNA < uniqueFeatureCountRange[2] &
  #                   percent.mt < mtPercent)
  # object = subset(object = object, cells = cellsToKeep)
  
  if(!justRed) {
    print("Filtering cells for number of features and mitochondrial percentages...")
    if(isTPM | prenormalized) {
      object = object[, which(object@meta.data$percent.mt < mtPercent)]
    } else {
      object = object[, which(object@meta.data$nFeature_RNA > uniqueFeatureCountRange[1] &
                              object@meta.data$nFeature_RNA < uniqueFeatureCountRange[2] &
                              object@meta.data$percent.mt < mtPercent)]
    }
  }

  if(!justRed) {
    if((is.null(object@commands$NormalizeData.RNA) & !prenormalized) | reprocess) {
      if(isTPM) {
        object@assays$RNA@data = log1p(object@assays$RNA@counts)
      } else {
        object = NormalizeData(object = object)
      }
    }
  }
  
  if(!justRed) {
    if(is.null(object@commands$FindVariableFeatures.RNA) | reprocess) {
      object = FindVariableFeatures(object = object, assay = assay)
    }
  }
  
  geneNames = NULL
  if(allGenes) geneNames = rownames(x = object)
   
  if(is.null(object@commands$ScaleData.RNA) | reprocess) {
    object = ScaleData(object = object, features = geneNames, assay = assay)
  }
  
  
  if(is.null(object@commands$RunPCA.RNA) | reprocess) {
    object = RunPCA(object = object, 
                    features = geneNames, assay = assay)
  }
  
  if(is.null(object@commands$FindNeighbors.RNA.pca) | reprocess) {
    object = FindNeighbors(object = object, features = geneNames, assay = assay)
  }
  
  if(is.null(object@commands$FindClusters) | reprocess) {
    for(i in res) object = FindClusters(object = object, resolution = i)
  }
  
  if(is.null(object@commands$RunUMAP.RNA.pca) | reprocess) {
    object = RunUMAP(object = object, dims = 1:10, assay = assay)
  }
  
  if(is.null(object@commands$RunTSNE) | reprocess) {
    object = RunTSNE(object = object)
  }
  
  # Idents(object = object) = "orig.ident"
  
  return(object)
}

# Processes scRNA-seq data for spike-in gene-based batch correction
SpikeInProcess = function(counts, spikeinGenes,
                          uninterestingVars = NULL,
                          seed = 0) {
  # load("/root/ec2_volume/ernest_09102019/TE_vs_placenta_09102019/IS021.seurat.RData")
  # counts = as.matrix(IS021.seurat@assays$RNA@counts)
  # spikeinGenes = c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62',
  #                  'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX',
  #                  'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1',
  #                  'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14',
  #                  'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP',
  #                  'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A',
  #                  'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544',
  #                  'C1orf43', 'CHMP2A', 'EMC7', 'GPI', 'PSMB2', 'PSMB4', 'RAB7A', 'REEP5',
  #                  'SNRPD3', 'VCP', 'VPS29')
  # sampleNames = IS021.seurat@meta.data$orig.ident
  # sampleNames = gsub(pattern = "^.*D[03].*$", replacement = "Early", x = sampleNames)
  # sampleNames = gsub(pattern = "^.*D[69].*$", replacement = "Late", x = sampleNames)
  # 
  # uninterestingVars = data.frame(sample = sampleNames)
  # uninterestingVars = data.frame(sample = sampleNames,
  #                                batch = sample(x = c("batch1", "batch2", "batch3"),
  #                                               size = length(sampleNames),
  #                                               replace = T),
  #                                day = sample(x = c("day1", "day2", "day3"),
  #                                             size = length(sampleNames),
  #                                             replace = T))
  # seed = 0
  # counts = as.matrix(mergedSeurat@assays$RNA@counts)
  # spikeinGenes = scan(file = "TE_vs_placenta_09102019/spikeinGenes_09132019.txt",
  #                     what = "character")
  # spikeinGenes = spikeinGenes
  # batchNamesMerged = gsub(pattern = "^P.*", replacement = "placenta", 
  #                         x = mergedSeurat@meta.data$orig.ident)
  # batchNamesMerged = gsub(pattern = "^H.*", replacement = "TE", 
  #                         x = batchNamesMerged)
  # uninterestingVars = data.frame(batch = batchNamesMerged)
  # seed = 0
  
  spikeinGenes = unique(spikeinGenes[spikeinGenes %in% rownames(counts)])
  
  # batchNames = gsub(pattern = "^H.*", replacement = "SCTL", x = Idents(mergedSeurat))
  # batchNames = gsub(pattern = "^P.*", replacement = "Placenta", x = batchNames)
  
  # Identify cells whose read count for every spike-in gene is 0
  spikeinZero = colSums(x = counts[spikeinGenes, ]) == 0
  
  # Create SingleCellExperiment object
  counts = counts[, !spikeinZero]
  sceInputList = list(counts = counts)
  
  uninterestingVars = uninterestingVars[!spikeinZero, , drop = F]
  
  sceColData = as(uninterestingVars, "DataFrame")
  
  sceRowData = DataFrame(Feature = rownames(counts))
  
  sce = SingleCellExperiment(sceInputList,
                             colData = sceColData,
                             rowData = sceRowData)
  isSpike(sce, "spike") = rownames(counts) %in% spikeinGenes
  
  # Quality control and normalization
  sce = calculateQCMetrics(object = sce, compact = T)
  qc = sce$scater_qc
  lowLib = isOutlier(metric = qc$all$log10_total_counts, type = "lower", nmad = 3)
  lowGenes = isOutlier(metric = qc$all$log10_total_features_by_counts, type = "lower", nmad = 3)
  highSpike = isOutlier(metric = qc$feature_control_spike$pct_counts, type = "higher", nmad = 3)
  # data.frame(LowLib = sum(lowLib), LowNgenes = sum(lowGenes), 
  #            HighSpike = sum(highSpike, na.rm = TRUE))
  discard = lowLib | lowGenes | highSpike
  sce = sce[,!discard]
  # summary(discard)
  
  # Compute size factors
  set.seed(seed)
  clusters = quickCluster(x = sce, use.ranks = F, BSPARAM = IrlbaParam())
  #table(clusters)
  
  sce = computeSumFactors(x = sce, min.mean = 0.1, clusters = clusters)
  #summary(sizeFactors(sce))
  
  sce = computeSpikeFactors(x = sce, general.use = F)
  #summary(sizeFactors(sce, "spike"))
  
  # Normalize
  sceUnlogged = normalize(object = sce, return_log = F)
  
  sce = normalize(object = sce, return_log = T)
  
  # Identify highly variable genes, block to ensure that uninteresting
  # differences between batches and samples do not inflate the variance
  block = apply(X = as.data.frame(sce@colData[, colnames(uninterestingVars)]), 
                MARGIN = 1, FUN = paste, collapse = "_")
  #block = paste0(sce$Batch, "_", sce$Sample)
  
  fit = trendVar(x = sce, block = block, parametric = T) 
  dec = decomposeVar(x = sce, fit = fit)
  
  # plot(dec$mean, dec$total, xlab="Mean log-expression",
  #      ylab="Variance of log-expression", pch=16)
  # isSpikeIn = isSpike(sce)
  # points(dec$mean[isSpikeIn], dec$total[isSpikeIn], col="red", pch=16)
  # curve(fit$trend(x), col="dodgerblue", add=TRUE)
  
  # Order genes by biological component
  dec$Feature = rowData(sce)$Feature
  dec = dec[order(dec$bio, decreasing = T),]
  
  return(list(sce = sce, sceUnlogged = sceUnlogged, dec = dec, fit = fit))
  
}

# Returns vector of p-values based on 
require(limma)
ReductionGeneSetEnrichment = function(testGenes, reductionDF, alternative = "mixed") {
  # testGenes = as.character(ccGenes)
  # reductionDF = as.data.frame(mergedSeurat@reductions$pca@feature.loadings)
  # alternative = "mixed"
  
  testGenes = testGenes[testGenes %in% rownames(reductionDF)]
  returnVec = vector(mode = "numeric", length = ncol(reductionDF))
  names(returnVec) = colnames(reductionDF)
  testGenesIndices = which(rownames(reductionDF) %in% testGenes)
  for(i in 1:ncol(reductionDF)) {
    returnVec[i] = geneSetTest(index = testGenesIndices, statistics = reductionDF[, i],
                               alternative = alternative)
  }
  
  return(returnVec)
}

# Takes spikein genes and counts matrix and uses scran to return size factors for spikein genes
require(scran)
require(scater)
SpikeFactors = function(counts, spikeinGenes) {
  # counts = as.matrix(mergedSeuratOrig@assays$RNA@counts)
  # spikeinGenes = scan("TE_vs_placenta_09102019/spikeinGenes_09132019.txt", what = "character")
  
  spikeinGenes = unique(spikeinGenes[spikeinGenes %in% rownames(counts)])
  
  # Identify cells whose read count for every spike-in gene is 0
  spikeinZero = colSums(x = counts[spikeinGenes, ]) == 0
  
  # Create SingleCellExperiment object
  counts = counts[, !spikeinZero]
  sceInputList = list(counts = counts)
  
  sceRowData = DataFrame(Feature = rownames(counts))
  
  sce = SingleCellExperiment(sceInputList,
                             rowData = sceRowData)
  isSpike(sce, "spike") = rownames(counts) %in% spikeinGenes
  
  # Quality control and normalization
  sce = calculateQCMetrics(object = sce, compact = T)
  qc = sce$scater_qc
  lowLib = isOutlier(metric = qc$all$log10_total_counts, type = "lower", nmad = 3)
  lowGenes = isOutlier(metric = qc$all$log10_total_features_by_counts, type = "lower", nmad = 3)
  highSpike = isOutlier(metric = qc$feature_control_spike$pct_counts, type = "higher", nmad = 3)
  discard = lowLib | lowGenes | highSpike
  sce = sce[,!discard]
  
  spikeFactors = computeSpikeFactors(x = sce, general.use = F, sf.out = T)
  names(spikeFactors) = colnames(sce)
  return(spikeFactors)
}

# Runs UMAP on seurat object with given number of variable features
UMAP = function(object, hvgs = 2000) {
  object = RunPCA(object, features = VariableFeatures(object)[1:hvgs])
  object = FindNeighbors(object, features = VariableFeatures(object)[1:hvgs])
  object = FindClusters(object)
  object = RunUMAP(object, dims = 1:10)
  Idents(object) = "orig.ident"
  return(object)
}

# Sort output data frame from FindMarkers in a more biologically meaningful way
SortMarkers = function(x) {
  x[order(x$p_val, 1/abs(x$avg_logFC), 1/abs(x$pct.1 - x$pct.2), decreasing = F), ]
}

# Match "uninteresting" features, such as ribosomal protein genes, LINCs, etc
MatchUnknownFeats = function(x, value = F, returnBad = T) {
  ribosomal = grepl(pattern = "^RP[SL].*", x = x, ignore.case = F)
  linc = grepl(pattern = "^LINC.*", x = x, ignore.case = F)
  linc = linc | grepl(pattern = "^[A-Z][A-Z]+\\d+\\.\\d+", x = x, ignore.case = F)
  orf = grepl(pattern = "^.*orf.*", x = x, ignore.case = F)
  mitochondrial = grepl(pattern = "^MT-.*", x = x, ignore.case = F)
  if(returnBad) {
    if(value) return(x[ribosomal | linc | orf | mitochondrial])
    return(ribosomal | linc | orf | mitochondrial) 
  }
  if(value) return(x[!(ribosomal | linc | orf | mitochondrial)])
  return(!(ribosomal | linc | orf | mitochondrial))
}


# Creates sample table for singleCellNet from Seurat object
SeuratToSCNSampleTable = function(seuratObject, cell_type = Idents(seuratObject)) {
  # seuratObject = teSeurat
  # clusterIdent = "seurat_clusters"
  # description1 = Idents(seuratObject)
  # description2 = seuratObject@meta.data$orig.ident
  # df = as.data.frame(matrix(nrow = ncol(seuratObject), ncol = 4))
  df = as.data.frame(matrix(nrow = ncol(seuratObject), ncol = 2))
  rownames(df) = colnames(seuratObject)
  # colnames(df) = c("sample_name", "Cluster_Number", "description1", "description2")
  colnames(df) = c("sample_name", "cell_type")
  df[, "sample_name"] = rownames(df)
  # df[, "Cluster_Number"] = seuratObject[[clusterIdent]]
  df[, "cell_type"] = cell_type
  # df[, "description2"] = description2
  
  return(df)
}

# Edited version of plot_metrics from singleCellNet
require(gridExtra)
plot_metrics2 = function (assessed) {
  # assessed = tm_heldoutassessment
  metric <- matrix(0, ncol = 5, nrow = 1)
  colnames(metric) <-
    c("cohen's kappa",
      "accuracy",
      "multiLogLoss",
      "mean_AUPRC",
      "weighted-AUPRC")
  rownames(metric) <- "value"
  metric[, 1:5] <-
    c(
      assessed$kappa,
      assessed$accuracy,
      assessed$multiLogLoss,
      assessed$AUPRC_w,
      assessed$AUPRC_wc
    )
  metric <- as.data.frame(metric)
  p1 <- ggplot(metric, aes(x = "cohen's kappa", y = metric[1,
                                                           1])) + geom_bar(stat = "identity") + xlab("") + ylab("") +
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8)) +
    ylim(0, 1) + theme(legend.position = "none")
  p4 <- ggplot(metric, aes(x = "mean_AUPRC", y = metric[1,
                                                        4])) + geom_bar(stat = "identity") + xlab("") + ylab("") +
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8)) +
    ylim(0, 1) + theme(legend.position = "none")
  grid.arrange(p1, p4, ncol = 2)
  #(p1 | p4)
}


# Runs SingleCellNet training, validation, and prediction steps
RunSCN = function(trainExp, trainST, queryExp, queryST,
                  dLevelTrain = "cell_type", dLevelQuery = "cell_type",
                  sampleColTrain = "sample_name", sampleColQuery = "sample_name",
                  ncellsTrain = 100, nTopGenes = 10, nRandTrain = 70,
                  nTrees = 1000, nTopGenePairs = 25, nCellsVal = 100,
                  nRandPred = 50, nRandAssess = 50, numPairs = 20,
                  nRandQuery = 50, secondSplit = T) {
  # setwd("/root/ec2_volume/ernest_09102019/")
  # mergedData = readRDS("TE_vs_placenta_09102019/seurat_processed_cell-labels_10082019.rds")
  # trainExp = mergedData@assays$RNA@counts[, grep("^P", mergedData$orig.ident)]
  # trainST = read.csv("singleCellNet_testing_10072019/placentaTrainingData_10072019.csv",
  #                    row = 1, as.is = T)
  # queryExp = mergedData@assays$RNA@counts[, grep("^H", mergedData$orig.ident)]
  # queryST = SeuratToSCNSampleTable(mergedData[, grep("^H", mergedData$orig.ident)])
  
  # trainExp = allExp
  # trainST = allST
  # queryExp = teCounts
  # queryST = teST
  # trainExp = zhouIntCounts
  # trainST = zhouIntST
  # queryExp = teIntCounts
  # queryST = teIntST
  # dLevelTrain = "cell_type"; dLevelQuery = "cell_type";
  # sampleColTrain = "sample_name"; sampleColQuery = "sample_name";
  # ncellsTrain = 100; nTopGenes = 10; nRandTrain = 70;
  # nTrees = 1000; nTopGenePairs = 25; nCellsVal = 100;
  # nRandPred = 50; nRandAssess = 50; numPairs = 20;
  # nRandQuery = 50
  # secondSplit = F
  
  commonGenes = intersect(rownames(trainExp), rownames(queryExp))
  trainExp = trainExp[commonGenes, ]
  trainSTSplit = splitCommon(sampTab = trainST, ncells = ncellsTrain, dLevel = dLevelTrain)
  trainSTTrain = trainSTSplit[["train"]]
  trainExpTrain = trainExp[, rownames(trainSTTrain)]
  class_info = scn_train(stTrain = trainSTTrain, expTrain = trainExpTrain, 
                         nTopGenes = nTopGenes, 
                         nRand = nRandTrain, nTrees = nTrees, nTopGenePairs = nTopGenePairs, 
                         dLevel = dLevelTrain, colName_samp = sampleColTrain)
  
  #validate data
  if(secondSplit) {
    stTestSplit = splitCommon(sampTab = trainSTSplit[["val"]], ncells = nCellsVal, 
                              dLevel = dLevelTrain) 
    stTest = stTestSplit[["train"]]
    expTest = trainExp[commonGenes, rownames(stTest)]
  } else {
    stTest = trainSTSplit[["val"]]
    expTest = trainExp[commonGenes, rownames(stTest)]
  }
  
  #predict
  classRes_val_all = scn_predict(cnProc = class_info[['cnProc']], expDat = expTest, 
                                 nrand = nRandPred)
  
  # Assess classifier
  tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, 
                                     stTrain = trainSTTrain, 
                                     stQuery = stTest, dLevelSID = sampleColTrain, 
                                     classTrain = dLevelTrain, 
                                     classQuery = dLevelTrain, 
                                     nRand = nRandAssess)
  
  #create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
  valGroups = as.vector(stTest[, dLevelTrain])
  names(valGroups) = as.vector(stTest[, sampleColTrain])
  valGroupsRand = rep("rand", nRandPred) 
  names(valGroupsRand) = paste("rand_", 1:nRandPred, sep = "")
  valGroups = append(valGroups, valGroupsRand) #include in the random cells profile created
  
  gpTab = compareGenePairs(query_exp = expTest, training_exp = trainExpTrain, 
                           training_st = trainSTTrain, 
                           classCol = dLevelTrain, sampleCol = sampleColTrain, 
                           RF_classifier = class_info$cnProc$classifier, 
                           numPairs = numPairs, trainingOnly = T)
  
  train = findAvgLabel(gpTab = gpTab, stTrain = trainSTTrain, dLevel = dLevelTrain)
  
  # Apply to query data
  crQueryAll = scn_predict(class_info[['cnProc']], queryExp, nrand = nRandQuery)
  
  predGroups = as.vector(queryST[, dLevelQuery])
  names(predGroups) = as.vector(queryST[, sampleColQuery])
  predGroupsRand = rep("rand", nRandQuery)
  names(predGroupsRand) = paste("rand_", 1:nRandQuery, sep = "")
  predGroups = append(predGroups, predGroupsRand)
  
  return(list(tm_heldoutassessment = tm_heldoutassessment, 
              classRes_val_all = classRes_val_all,
              stTest = stTest,
              valGroups = valGroups,
              nRandQuery = nRandQuery, 
              gpTab = gpTab,
              train = train,
              class_info = class_info,
              predGroups = predGroups,
              crQueryAll = crQueryAll))
}

# Returns common values from any number of vectors
Common = function(...) {
  argList = list(...)
  commonValues = argList[[1]]
  for(i in 2:length(argList)) {
    commonValues = intersect(commonValues, argList[[i]])
  }
  return(commonValues)
}

# Returns values that are not unique to one vector
NonUnique = function(...) {
  argList = list(...)
  allVals = unlist(argList)
  dupIndices = duplicated(allVals)
  if(all(!dupIndices)) return(NULL)
  return(allVals[dupIndices])
}

# Takes character vector and returns integer vector with input values as names
Integers = function(x) {
  xUnique = unique(x)
  refVector = 1:length(xUnique)
  names(refVector) = xUnique
  return(refVector[x])
}

# Takes a list of Seurat objects and prepares data for BERMUDA
BermudaPrep = function(seuratList, outputCountFileNames, mnOutputFile) {
  # load("TE_vs_placenta_09102019/IS021.seurat.RData")
  # load("TE_vs_placenta_09102019/placenta.reorder.Rda")
  # teSeurat = IS021.seurat
  # teSeurat[["cellLabels"]] = teSeurat@meta.data$orig.ident
  # plSeurat = placenta.reorder
  # plSeurat[["cellLabels"]] = Idents(plSeurat)
  # plSeurat = subset(plSeurat, subset = orig.ident %in% c("P17","P23"))
  # seuratList = list(te = teSeurat, pl = plSeurat)
  # outputCountFileNames = c("TE_vs_placenta_09102019/bermuda_10152019/bermuda_data_PL_10162019.csv",
  #                          "TE_vs_placenta_09102019/bermuda_10152019/bermuda_data_TE_10162019.csv")
  # mnOutputFile = "TE_vs_placenta_09102019/bermuda_10152019/bermuda_metaneighbor_10162019.csv"
  
  for(i in 1:length(seuratList)) seuratList[[i]] = SeuratProcess(seuratList[[i]]) 
  
  varFeats = VariableFeatures(seuratList[[1]])
  for(i in 2:length(seuratList)) varFeats = intersect(varFeats, VariableFeatures(seuratList[[i]]))
  
  numCells = sum(sapply(seuratList, ncol))
  
  for(i in 1:length(seuratList)) seuratList[[i]] = seuratList[[i]][varFeats, ]
  
  logData = as.matrix(seuratList[[1]]@assays$RNA@data)
  dimnames(logData) = dimnames(seuratList[[1]])
  for(i in 2:length(seuratList)) {
    tmp = as.matrix(seuratList[[i]]@assays$RNA@data)
    colnames(tmp) = colnames(seuratList[[i]])
    logData = cbind(logData, tmp)
  }
  
  clusterLabelList = lapply(seuratList, function(i) 0)
  cellLabelList = lapply(seuratList, function(i) 0)
  datasetLabelList = lapply(seuratList, function(i) 0)
  for(i in 1:length(seuratList)) clusterLabelList[[i]] = as.integer(seuratList[[i]]$seurat_clusters) + max(unlist(clusterLabelList))
  for(i in 1:length(seuratList)) cellLabelList[[i]] = Integers(seuratList[[i]]$cellLabels) + max(unlist(cellLabelList))
  for(i in 1:length(seuratList)) datasetLabelList[[i]] = rep(i, ncol(seuratList[[i]]))
  
  for(i in 1:length(seuratList)) {
    write_dataset_cluster(filename = outputCountFileNames[i],
                          data = as.matrix(seuratList[[i]]@assays$RNA@counts),
                          sample_labels = cellLabelList[[i]],
                          cell_labels = cellLabelList[[i]],
                          cluster_labels = clusterLabelList[[i]])
  }
  
  pheno = data.frame(Celltype = as.character(unlist(clusterLabelList)),
                     Study_ID = as.character(unlist(datasetLabelList)),
                     stringsAsFactors = F)
  # run metaneighbor
  clusterLabels = unique(unlist(clusterLabelList))
  clusterSimilarity = run_MetaNeighbor_US(varFeats, logData, clusterLabels, pheno)
  
  # set cluster pairs from the same dataset to 0
  for (i in 1:length(seuratList)) { 
    clustIndices = unique(clusterLabelList[[i]])
    clusterSimilarity[clustIndices, clustIndices] = 0
  }
  
  # order rows and columns
  clusterSimilarity = clusterSimilarity[order(as.numeric(rownames(clusterSimilarity))),]
  clusterSimilarity = clusterSimilarity[,order(as.numeric(colnames(clusterSimilarity)))]
  
  # write out metaneighbor file
  write.table(clusterSimilarity, mnOutputFile, sep = ",", quote = F, col.names = T, row.names = F)
  
  return(list(seuratList = seuratList, logData = logData, clusterLabelList = clusterLabelList,
              cellLabelList = cellLabelList, datasetLabelList = datasetLabelList, clusterSimilarity = clusterSimilarity))
} 

# Run FindMarkers and SortMarkers on multiple comparisons
FindMarkersMulti = function(seuratObject, compList, minPct = -Inf, 
                            logFC = 0.25, minDiffPct = -Inf, verbose = T) {
  # seuratObject = mergedSeuratProcessed
  # compList = compList[1:2]
  
  origIdents = as.character(Idents(seuratObject))
  
  resultsList = vector("list", length = length(compList))
  names(resultsList) = names(compList)  

  for(group1 in names(compList)) {
    for(group2 in compList[[group1]]) {
      if(verbose) cat("Comparing ", group1, " vs. ", paste0(group2, collapse = ", "), "\n", sep = "")
      group2Label = paste0(group2, collapse = ".")
      tmpLabels = replace(origIdents, list = origIdents %in% group2, values = group2Label)
      Idents(seuratObject) = tmpLabels
      markers = FindMarkers(seuratObject, ident.1 = group1, ident.2 = group2Label,
                            logfc.threshold = logFC, min.pct = minPct, min.diff.pct = minDiffPct)
      markers = SortMarkers(markers)
      resultsList[[group1]][[group2Label]] = markers
    }
  }
  
  return(resultsList)    
}

# Return date in mmddyyyy format
Date = function() {
  format(Sys.Date(), "%m%d%Y")
}

# Output feature plots from FindMarkersMulti results
MarkerPlots = function(seuratObject, markerResults, outputBase, nTopFeats = 16) {
  # seuratObject = mergedSeuratProcessed
  # markerResults = multiCompMarkers
  # outputBase = "TE_vs_placenta_09102019/featplots_TE-PL"
  # nTopFeats = 16
  
  for(group1 in names(markerResults)) {
    for(group2 in names(markerResults[[group1]])) {
      markers = markerResults[[group1]][[group2]]
      fileName1 = paste0(outputBase, "_", group1, "-vs-", group2, "_", group1, "_", Date(), ".png")
      fileName2 = paste0(outputBase, "_", group1, "-vs-", group2, "_", group2, "_", Date(), ".png")
      group1Markers = rownames(markers)[markers$avg_logFC > 0]
      group2Markers = rownames(markers)[markers$avg_logFC < 0]
      png(fileName1, width = 800, height = 600)
      p = FeaturePlot(seuratObject, features = group1Markers[1:nTopFeats]) + 
        labs(title = paste0(group1, " vs. ", group2, ", ", group1, " markers")) + 
        theme(plot.title = element_text(hjust = 0.5, vjust = 1, face = "bold", size = 30),
              plot.margin = unit(c(1,1,1,1), "cm"))
      print(p)
      dev.off()
      png(fileName2, width = 800, height = 600)
      p = FeaturePlot(seuratObject, features = group2Markers[1:nTopFeats]) +
        labs(title = paste0(group1, " vs. ", group2, ", ", group2, " markers")) + 
        theme(plot.title = element_text(hjust = 0.5, vjust = 1, face = "bold", size = 30),
              plot.margin = unit(c(1,1,1,1), "cm"))
      print(p)
      dev.off()
    }
  }

}

# Write FindMarkers results from FindMarkersMulti to csv files
OutputMultiMarkers = function(markerResults, outputBase, sep = ",") {
  # markerResults = multiCompMarkers
  # outputBase = "TE_vs_placenta_09102019/markers_TE-PL"
  
  ext = ifelse(sep == ",", yes = ".csv", no = ".txt")
  
  for(group1 in names(markerResults)) {
    for(group2 in names(markerResults[[group1]])) {
      markers = markerResults[[group1]][[group2]]
      fileName = paste0(outputBase, "_", group1, "-vs-", group2, "_", Date(), ext)
      write.table(markers, file = fileName, quote = F, row.names = T, col.names = T, sep = sep)
      cat("Exported to ", fileName, "\n", sep = "")
    }
  }
}

# Run WGCNA workflow
RunWGCNA = function(expData, power, outputBase, blockSize = 5000,
                    dendPlotFileName = NULL,
                    dendTitle = NULL, 
                    minModuleSize = min(20, ncol(expData)/2 ),
                    verbosity = 3) {
  # expData = clean_te_pl_int$expData
  # power = 5
  # outputBase = "TE_vs_placenta_09102019/WGCNA_TE-PL_Int_AllVarFeats"
  # dendPlotFileName = NULL
  if(!is.null(blockSize)) {
    if(blockSize == "all") blockSize = ncol(expData)
  }
  
  net = blockwiseModules(expData, power = power, maxBlockSize = blockSize,
                         TOMType = "unsigned", minModuleSize = minModuleSize,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = outputBase, 
                         verbose = verbosity)
  
  
  mergedColors = labels2colors(net$colors)
  if(is.null(dendPlotFileName)) dendPlotFileName = paste0(outputBase, "_Dendrogram_", Date(), ".png")
  png(dendPlotFileName, width = 1200, height = 600)
  par(font = 2, cex.main = 2, cex.lab = 2, cex.axis = 2)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      cex.colorLabels = 1.5,
                      cex.dendroLabels = 0.5,
                      cex.rowText = 0.5,
                      marAll = c(1,8,3,1))
  if(!is.null(dendTitle)) {
    mtext(dendTitle, line = 10, cex = 2, font = 2)
  }
  dev.off()
  
  return(list(net = net, mergedColors = mergedColors))
  
}

# Clean data if necessary, run pickSoftThreshold, and plot results
CleanAndSFT = function(expData, outputFile, powers = c(1:10, seq(from = 12, to=20, by=2)), 
                       blockSize = NULL, verbose = 5) {
  # expData = mergedExpData
  
  gsg = goodSamplesGenes(expData, verbose = 3)
  if(!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(expData)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(expData)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    expData = expData[gsg$goodSamples, gsg$goodGenes]
  }

  if(!is.null(blockSize)) {
    if(blockSize == "all") blockSize = ncol(expData)
  }
  sft = pickSoftThreshold(expData, powerVector = powers, blockSize = blockSize, 
                          verbose = verbose)
  png(outputFile, width = 960, height = 600)
  par(mfrow = c(1,2))
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  return(list(expData = expData, sft = sft))
}

# Run Seurat integration approach
RunSeuratIntegration = function(seuratList, nAnchorFeatures = 2000, reductionMethod = "cca") {
  anchors = FindIntegrationAnchors(seuratList, dims = 1:30, anchor.features = nAnchorFeatures, reduction = reductionMethod)
  seuratIntegrated = IntegrateData(anchors, dims = 1:30)
  DefaultAssay(seuratIntegrated) = "integrated"
  
  seuratIntegrated = ScaleData(seuratIntegrated)
  seuratIntegrated = RunPCA(seuratIntegrated)
  seuratIntegrated = RunUMAP(seuratIntegrated, dims = 1:30)
  return(seuratIntegrated)
}

# Replaces values in vector matching original vector with those matching new vector
Replace = function(vec, origValues, newValues) {
  for(i in 1:length(origValues)) {
    vec[vec == origValues[i]] = newValues[i]
  }
  return(vec)
}

# Plot multiple dendrograms from blockwise module detection in one plot
PlotMultiDend = function(net) {
  nDends = length(net$net$dendrograms)
  nRows = ifelse(nDends %% 2 == 0, yes = nDends*2, no = (nDends + 1)*2)
  layout(matrix(1:nRows,ncol=2),  heights = rep.int(c(7,3), times = nDends))
  for(i in 1:nDends) {
    plotDendroAndColors(net$net$dendrograms[[i]], 
                        net$mergedColors[net$net$blockGenes[[i]]],
                        "Module colors",
                        setLayout = F,
                        main = paste0("Block ", i),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        cex.colorLabels = 1.5,
                        cex.dendroLabels = 0.5,
                        cex.rowText = 0.5,
                        marAll = c(1,10,3,1))
  }

}

combineResults <- function(enr) {
  col_names = c("libName","lib_rank", "gene_count", "term", "overlap","pval", "adjPval", "oldPval","oldAdjPval","Z_score","score","gene_list")
  fullRes <- data.frame(libName=character(), lib_rank=integer(), gene_count=integer(), term=character(), overlap=character(),pval=double(), adjPval=double(),
                        oldPval=double(),oldAdjPval=double(),Z_score=double(),score=double(),gene_list=character(), stringsAsFactors=FALSE)
  print(length(libs))
  for (i in 1:length(libs)) {
    res<-data.frame(enr[i])
    print(dim(res)[1])
    for(j in 1:(dim(res)[1])) {
      currdf <- data.frame(c(libs[i], j, as.integer(unlist(strsplit(res[j,2],"/"))[1]), res[j,]))
      names(currdf) <- col_names
      fullRes <- rbind(fullRes, currdf)
    }
  }
  return(fullRes)
}

getSigTerms <- function(enr) {
  col_names = c("libName","lib_rank", "gene_count", "term", "overlap","pval", "adjPval", "oldPval","oldAdjPval","Z_score","score","gene_list")
  fullSigRes <- data.frame(libName=character(), lib_rank=integer(), gene_count=integer(), term=character(), overlap=character(),pval=double(), adjPval=double(),
                           oldPval=double(),oldAdjPval=double(),Z_score=double(),score=double(),gene_list=character(), stringsAsFactors=FALSE)
  print(length(libs))
  for (i in 1:length(libs)) {
    res<-data.frame(enr[i])
    print(dim(res)[1])
    for(j in 1:(dim(res)[1])) {
      if(res[j,4]<0.1) {				
        currdf <- data.frame(c(libs[i], j, as.integer(unlist(strsplit(res[j,2],"/"))[1]), res[j,]))
        names(currdf) <- col_names
        fullSigRes <- rbind(fullSigRes, currdf)
      }
    }
  }
  return(fullSigRes)
}

getTop10SigTerms <- function(enr) {
  col_names = c("libName","lib_rank", "gene_count", "term", "overlap","pval", "adjPval", "oldPval","oldAdjPval","Z_score","score","gene_list")
  topSigRes <- data.frame(libName=character(), lib_rank=integer(), gene_count=integer(), term=character(), overlap=character(),pval=double(), adjPval=double(),
                          oldPval=double(),oldAdjPval=double(),Z_score=double(),score=double(),gene_list=character(), stringsAsFactors=FALSE)
  for (i in 1:length(libs)) {
    #loop over result, append lib name to front, save only adjP <0.1
    res<-data.frame(enr[i])
    rowCnt<-dim(res)[2]
    # only consider the top 10 scores, only take adjP < 0.1
    for(j in 1:min(10,rowCnt)) {
      print(paste("**",res[j,4],"**"))
      if((!is.na(res[j,4])) & (res[j,4]<0.1)) {
        currdf <- data.frame(c(libs[i], j, as.integer(unlist(strsplit(res[j,2],"/"))[1]), res[j,]))
        names(currdf) <- col_names
        topSigRes <- rbind(topSigRes, currdf)
      }
    }
  }
  return(topSigRes)
}

unpivotGeneLists <- function(fsr) {
  geneFrame<- data.frame(libName=character(), term=character(), adjP=double(), comb_score=double(), gene_symbol=character(), stringsAsFactors=FALSE)
  currRow <- 1
  mainRow <- 1
  for(r in 1:(dim(fsr)[1])){
    genes <- unlist(strsplit(as.character(fsr[r,"gene_list"]),';'))		
    for(gene in 1:length(genes)){
      geneFrame[mainRow,1]<-as.vector(fsr[currRow,"libName"])
      geneFrame[mainRow,2]<-as.vector(fsr[currRow,"term"])
      geneFrame[mainRow,3]<-fsr[currRow,"adjPval"]
      geneFrame[mainRow,4]<-fsr[currRow,"score"]
      geneFrame[mainRow,5]<-genes[gene]
      mainRow <- mainRow + 1
    }
    currRow <- currRow + 1
  }
  return(geneFrame)
}

EnrichR = function(genes, libFile, nGenes = NULL) {
  # genes = upMarkers
  # libFile = "rscripts_09102019/enrichr_libraries_KEGG2019_GObio2018_ARCHS4.txt"
  # nGenes = NULL
  
  if(is.null(nGenes)) {
    nGenes = length(genes)
  } else {
    nGenes = min(nGenes, length(genes))
  }
  
  libSet <- read.table(libFile, header=FALSE)
  libs <- as.vector(libSet[,1])
  
  enr<-enrichr(genes, libs)
  
  fsr <- getSigTerms(enr)
  
  return(fsr) 
}

# Return rows of edges between cell groups that are not similar as defined by similarList
KeepSimilarGroupEdges = function(seuratData, edges, similarList) {
  badRows = c()
  cellIndices = 1:ncol(seuratData)
  for(i in 1:length(similarList)) {
    cellGroup = names(similarList)[i]
    similarStages = c(cellGroup, similarList[[cellGroup]])
    indices = which(Idents(seuratData) == cellGroup)
    similarIndices = which(Idents(seuratData) %in% similarStages)
    badRowsTemp = which(edges$V1 %in% indices & !(edges$V2 %in% similarIndices) |
                          edges$V2 %in% indices & !(edges$V1 %in% similarIndices))
    badRows = c(badRows, badRowsTemp)
  }
  return(badRows)
}

# Split cells from Seurat object into two groups based on expression of two genes
# method = "quantiles" 
# method = "stdev"
# method = "default"
SplitCellsByMarkers = function(seuratObject, pos = NULL, neg = NULL, barcodes = colnames(seuratObject), 
                               method = "quantiles", quantileTHs = c(0.5, 0.75),
                               gt = 1, le = 1, fc = 2,
                               slot = "scale.data") {
  if(is.null(pos) & is.null(neg)) return(list())
  
  x1 = x2 = transition = NULL
  expData = GetAssayData(seuratObject, slot = slot)
  gt1 = gt2 = gt
  le1 = le2 = le
  
  if(!is.null(pos)) {
    exp1 = as.numeric(expData[pos, barcodes])
    names(exp1) = barcodes
    exp1All = as.numeric(expData[pos, ])
    names(exp1All) = colnames(seuratObject)
    if(method == "quantiles") {
      gt1 = quantile(exp1All, quantileTHs[2])
      le1 = quantile(exp1All, quantileTHs[1])
    }
  }
  
  if(!is.null(neg)) {
    exp2 = as.numeric(expData[neg, barcodes])
    names(exp2) = barcodes
    exp2All = as.numeric(expData[neg, ])
    names(exp2All) = colnames(seuratObject)
    if(method == "quantiles") {
      gt2 = quantile(exp2All, quantileTHs[2])
      le2 = quantile(exp2All, quantileTHs[1])
    }
  }
  
  if(!is.null(pos) & !is.null(neg)) {
    x1 = exp1[exp1 > gt1]
    x1 = x1[exp2[names(x1)] <= le2 | x1/exp2[names(x1)] >= fc]
    x2 = exp2[exp2 > gt2]
    x2 = x2[exp1[names(x2)] <= le1 | x2/exp1[names(x2)] >= fc]
  } else if(!is.null(pos)) {
    x1 = exp1[exp1 > gt1]
  } else if(!is.null(neg)) {
    x1 = exp2[exp2 <= le2]
  }
  x1 = names(x1)
  x2 = names(x2)
  
  if(!is.null(pos) & !is.null(neg)) {
    transition = barcodes[exp1 > gt1 & exp2 > gt2 & !(barcodes %in% c(x1,x2))]
    transition = transition[
      exp1[transition]/exp2[transition] < fc & 
        exp2[transition]/exp1[transition] < fc
      ]
  }
  
  other = barcodes[!(barcodes %in% c(x1,x2,transition))]
  if(!all(sort(c(x1,x2,transition,other)) == sort(barcodes))) {
    return(NULL)
  }
  return(list(x1 = x1, x2 = x2, transition = transition, other = other))    
}

# Takes two lists with identical structure and return a table containing 
# pairwise overlap
OverlapTable = function(list1, list2) {
  mat = matrix(nrow = (length(list1)), ncol = (length(list2)))
  rownames(mat) = names(list1)
  colnames(mat) = names(list2)
  for(i in 1:(length(list1))) {
    for(j in 1:(length(list2))) {
      mat[i,j] = length(intersect(list1[[i]], list2[[j]]))
    }
  }
  return(mat)
}

# Performs Fisher's exact test on matrix containing overlap numbers and returns
# a list with a matrix of odds ratios and a matrix of p-values
FishTable = function(mat) {
  # mat = overlapTable
  tot = sum(mat)
  row_sums = rowSums(mat)
  col_sums = colSums(mat)
  outList = list(or = matrix(nrow = nrow(mat), ncol = ncol(mat)), p = matrix(nrow = nrow(mat), ncol = ncol(mat)))
  dimnames(outList$or) = dimnames(mat)
  dimnames(outList$p) = dimnames(mat)
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      cTable = matrix(c(mat[i,j], row_sums[i]-mat[i,j], col_sums[j], tot-col_sums[j]), nrow = 2)
      fishTest = fisher.test(cTable, alternative = "greater")
      outList$or[i, j] = fishTest$estimate
      outList$p[i, j] = fishTest$p.value
    }
  }
  return(outList)
}

# Classifies cells from Seurat object into cell types in sigList based on combinations of
# positive and negative expression of genes. sigList is a named list of lists. The names
# are the cell identities the cells will be grouped into. Each cell identity should have a 
# corresponding list containing pos, neg, and cluster. pos is the gene(s) that a cell should
# express, neg is one gene that the cell should express at a relatively low level. cluster
# is the seurat_clusters cluster(s) that the function should look in. If "quantiles" is the 
# method, cells will be identified as a particular cell identity if any of the pos genes has
# expression higher than the second quantile in quantileTHs (0.75 by default) and if either 
# expression of the neg gene is at or below the first quantile or the ratio of expression of
# the two genes is at or above fc. If quantiles is not chosen, the function uses 
# the thresholds "gt" and "le" in place of quantile values.
IDCells = function(seuratObject, sigList, method = "quantiles", quantileTHs = c(0.5, 0.75),
                   gt = 1, le = 1, fc = 2, slot = "scale.data") {
  # seuratObject = markerSeurat
  newIDs = rep.int("Other", times = ncol(seuratObject))
  names(newIDs) = colnames(seuratObject)
  cellIDs = names(sigList)
  for(cellID in cellIDs) {
    x1 = NULL
    seuratTemp = seuratObject[, seuratObject$seurat_clusters %in% sigList[[cellID]]$cluster]
    if(length(sigList[[cellID]]$pos) >= 1 & length(sigList[[cellID]]$neg) <= 1) {
      markerCellsList = lapply(sigList[[cellID]]$pos, 
                               function(i) SplitCellsByMarkers(seuratObject, 
                                                               pos = i, 
                                                               neg = sigList[[cellID]]$neg,
                                                               barcodes = colnames(seuratTemp),
                                                               method = "quantiles",
                                                               quantileTHs = quantileTHs,
                                                               gt = gt, le = le, fc = fc))
      nTests = length(sigList[[cellID]]$pos)
    } else {
      markerCellsList = lapply(sigList[[cellID]]$neg, 
                               function(i) SplitCellsByMarkers(seuratObject, 
                                                               pos = sigList[[cellID]]$pos,
                                                               neg = i,
                                                               barcodes = colnames(seuratTemp),
                                                               method = "quantiles",
                                                               quantileTHs = quantileTHs,
                                                               gt = gt, le = le, fc = fc))
      nTests = length(sigList[[cellID]]$neg)
    } 
    x1NonUnique = unlist(lapply(markerCellsList, function(i) i[["x1"]]))
    x1 = unique(x1NonUnique)
    counts = sapply(x1, function(i) sum(x1NonUnique == i))
    x1 = x1[counts == nTests]
    
    x2 = unique(unlist(lapply(markerCellsList, function(i) i[["x2"]])))
    transition = unique(unlist(lapply(markerCellsList, function(i) i[["transition"]])))
    dups = NonUnique(x1, x2, transition)
    x1 = x1[!(x1 %in% dups)]
    
    newIDs[x1] = cellID
  }
  return(newIDs)
}

# Example from https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2018-05-21/labs/igraph.html
# for creating knn graph
make.knn.graph<-function(D,k){
  # calculate euclidean distances between cells
  dist<-as.matrix(dist(D))
  # make a list of edges to k nearest neighbors for each cell
  edges <- mat.or.vec(0,2)
  for (i in 1:nrow(dist)){
    # find closes neighbours
    matches <- setdiff(order(dist[i,],decreasing = F)[1:(k+1)],i)
    # add edges in both directions
    edges <- rbind(edges,cbind(rep(i,k),matches))  
    edges <- rbind(edges,cbind(matches,rep(i,k)))  
  }
  # create a graph from the edgelist
  graph <- graph_from_edgelist(edges,directed=F)
  V(graph)$frame.color <- NA
  # make a layout for visualizing in 2D
  set.seed(1)
  g.layout<-layout_with_fr(graph)
  return(list(graph=graph,layout=g.layout))        
}

# Create list of graph and layout for use in plot.igraph
# coordinates: matrix with 2 columns, with each row containing 
# x and y coordinates for each cell. 
# edges: matrix with 2 columns, with each row containing a pair
# of integers indicating edges. Integers should correspond to 
# rows in coordinates
CreateGraph = function(coordinates, edges, seed = 1) {
  graph = graph_from_edgelist(edges, directed = F)
  set.seed(seed)
  graphLayout = coordinates
  return(list(graph = graph, layout = graphLayout))
}

# Generate vector of colors based on input categories
AssignColors = function(x, name = "Dark2", cols = NULL) {
  xUnique = unique(x)
  if(!is.null(cols)) {
    colorsUnique = cols[names(cols) %in% x]
  } else {
    colorsUnique = brewer.pal(length(xUnique), name = name)
    names(colorsUnique) = xUnique
  }
  return(list(def = colorsUnique, colors = colorsUnique[x]))
}

# Just hclust with Pearson
HClust = function(data, method = "pearson") {
  
  if(method %in% c("pearson", "spearman")) {
    options(warn = -1)
    my_cor = stats::cor(as.matrix(t(data)), method = method, use = "pairwise.complete.obs")
    options(warn = 0)
    
    my_cor = replace(my_cor, abs(my_cor) == 1, -1)
    my_cor = replace(my_cor, is.na(my_cor), -1)
    
    my_dist = as.dist(1 - my_cor)
    
  } else {
    my_dist = dist(data, method = method)
    my_dist = replace(my_dist, is.na(my_dist), max(my_dist, na.rm = T))
  }
  
  my_clust = hclust(my_dist)
  return(my_clust)
}

# Returns a dendrogram based on clustering using Pearson and hierarchical clustering
Dendrogram = function(data, method = "pearson") {
  my_clust = HClust(data, method = method)
  my_dend = as.dendrogram(my_clust)
  my_dend = reorder(my_dend, rowMeans(data, na.rm = T))
  return(my_dend)
}

# Write R-squared on plot. If value supplied, write R-squared = "value".
WriteR2 = function(pre = "", r2 = NULL, post = "") {
  if(is.null(r2)) {
    bquote(.(pre) ~ italic(R)^2 ~ .(post))
  } else {
    bquote(.(pre) ~ italic(R)^2 == .(format(r2, digits = 3)) ~ .(post))
  }
}

# Return the row and column of a matrix or data frame containing the 
# values defined by func (default = max). If more than one row/column,
# returns data frame
RowCol = function(x, func = max, excludeDiag = T) {
  x = as.matrix(x)
  if(excludeDiag) {
    diag(x) = NA
    x[lower.tri(x)] = NA
  }
  val = func(x, na.rm = T)
  eqVal = which(x == val)
  rows = eqVal %% nrow(x)
  cols = ceiling(eqVal / nrow(x))
  if(length(eqVal) == 1) return(c(rows,cols))
  return(data.frame(x = rows, y = cols))
}

# Return mediod indices from distance matrix
Meds = function(distMat, groups) {
  # distMat[lower.tri(distMat, diag = T)] = NA 
  diag(distMat) = NA
  uniqueGroups = unique(groups)
  names(groups) = colnames(distMat)
  meds = vector("integer", length = length(uniqueGroups))
  sampleIndices = 1:ncol(distMat)
  names(sampleIndices) = colnames(distMat)
  for(i in 1:length(uniqueGroups)) {
    distTemp = distMat[names(groups)[groups == uniqueGroups[i]], 
                       names(groups)[groups == uniqueGroups[i]]]
    rowSumsTemp = rowSums(distTemp, na.rm = T)
    meds[i] = sampleIndices[rownames(distTemp)[rowSumsTemp == min(rowSumsTemp)][1]]
  }
  return(meds)
}
