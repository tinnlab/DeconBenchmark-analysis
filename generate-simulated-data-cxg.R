#################
### This code is for generating simulated data using singcle cell data from CELLxGENE

library(Seurat)
library(tidyverse)

source('./include/utils.R')

metaFiles <- readRDS("./include/metaFiles.rds")

file_ids <- metaFiles$dataset_id[metaFiles$dataset == "CELLxGENE"]
generateCXG <- function(nBulk, nGene, savePath, appendBulkCount = FALSE){
  dir.create(savePath, showWarnings = FALSE, recursive = TRUE)

  for (file_id in file_ids) {
    file_path <- paste0("./data/", file_id, ".rds")
    tissue <- metaFiles$tissue[metaFiles$dataset_id == file_id]
    tissue_in_publication <- metaFiles$tissue_in_publication[metaFiles$dataset_id == file_id]

    data <- readRDS(file_path)

    # Clean metadata
    metadata <- data@meta.data
    metadata$assay <- as.character(metadata$assay)
    metadata$organism <- as.character(metadata$organism)
    metadata$tissue <- as.character(metadata$tissue)
    metadata$disease <- as.character(metadata$disease) %>% tolower()
    metadata$cell_type <- as.character(metadata$cell_type) %>% tolower() %>%
              str_replace_all(pattern = " ", replacement = "_") %>%
              str_replace_all(pattern = "-", replacement = "_") %>%
              str_replace_all(pattern = ",", replacement = "")
    metadata$sampleID <- as.character(rownames(metadata))

    # Filter metadata
    # Use 10x only
    metadata <- metadata[metadata$assay == "10x 3' v3",]
    # Use human data
    metadata <- metadata[metadata$organism == "Homo sapiens",]
    # Use tissue only
    metadata <- metadata[metadata$tissue == as.character(tissue),]
    # Use normal donors
    metadata <- metadata[metadata$disease == "normal",]
    #Remove any cell types that has number of cells  < 10
    tmpLabel <- metadata$cell_type
    toUseCellTypes <- names(table(tmpLabel))[table(tmpLabel) >= 10]
    metadata <- metadata[metadata$cell_type %in% toUseCellTypes,]
    idx <- paste0(metadata$donor_id, "___", metadata$cell_type)
    tmp <- names(table(idx))[table(idx) >= 10]
    metadata <- metadata[idx %in% tmp,]

    countdata <- GetAssayData(object = data, layer = "counts")
    tmpData <- countdata[,metadata$sampleID]
    metaDat <- metadata
    rm(countdata, metadata)

    # Remove all zero genes and then select top nGene
    tmpData <- tmpData[rowSums(tmpData) > 0, ]
    tmpRowSum <- rowSums(tmpData)
    tmpData <- as.matrix(tmpData[order(tmpRowSum, decreasing = TRUE)[1:min(nGene, nrow(tmpData))], ])

    trainDonor <- levels(metaDat$donor_id)[1]

    # Some lables
    tmpLabel <- as.character(metaDat$cell_type)
    tmpDonor <- as.character(metaDat$donor_id)

    trainData <- tmpData[, tmpDonor == trainDonor]
    trainLabel <- tmpLabel[tmpDonor == trainDonor]

    testData <- tmpData[, tmpDonor != trainDonor]
    testLabel <- tmpLabel[tmpDonor != trainDonor]

    simData <- dataGenerate(trainData, trainLabel, testData, testLabel, n = nBulk, minCell = 10, minMarkers = 1)

    print(tissue_in_publication)
    print(dim(simData$bulk))
    print(dim(simData$singleCellExpr))
    print(dim(simData$cellTypeExpr))
    print(dim(simData$signature))

    saveRDS(simData, paste0(savePath, tissue_in_publication, ifelse(appendBulkCount, paste0("_", nBulk), ""), '.rds'))

  }
}

generateCXG(nBulk = 100, nGene = 1e4, savePath = './data/sim-cxg/')