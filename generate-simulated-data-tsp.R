#################
### This code is for generating simulated data using singcle cell data from Tabula Sapiens

library(Seurat)
library(tidyverse)

source('./include/utils.R')

metaFiles <- readRDS("./include/metaFiles.rds")

## Working with Tabula Sapiens
file_id <- metaFiles$dataset_id[metaFiles$dataset == "Tabula-Sapiens"]
file_path <- paste0("./data/", file_id, ".rds")
data <- readRDS(file_path)

metaData <- data@meta.data
metaData$sampleID <- rownames(metaData)
dataMatrix <- GetAssayData(object = data, layer = "counts")

# Use 10x only
idx <- which(metaData$assay == "10x 3' v3")
metaData <- metaData[idx, ]
dataMatrix <- dataMatrix[, idx]

metaData <- metaData %>% group_split(tissue_in_publication)

generateTSP <- function(nBulk, nGene, savePath, appendBulkCount = FALSE){
    dir.create(savePath, showWarnings = FALSE, recursive = TRUE)
    for (metaDat in metaData) {
        tmpData <- dataMatrix[, metaDat$sampleID]

        # Remove all zero genes and then select top nGene
        tmpData <- tmpData[rowSums(tmpData) > 0, ]
        tmpRowSum <- rowSums(tmpData)
        tmpData <- as.matrix(tmpData[order(tmpRowSum, decreasing = TRUE)[1:min(nGene, nrow(tmpData))], ])

        # Some lables
        tmpLabel <- as.character(metaDat$cell_type) %>% str_replace_all(pattern = " ", replacement = "_") %>%
            str_replace_all(pattern = "-", replacement = "_") %>%
            str_replace_all(pattern = ",", replacement = "") # to prevent error when using design.matrix
        tmpDonor <- as.character(metaDat$donor)

        if(length(unique(tmpDonor)) > 1) { #if there are multiple donors, choose one randomly for reference matrix
            trainDonor <- unique(tmpDonor)[1]
            trainData <- tmpData[, tmpDonor == trainDonor]
            trainLabel <- tmpLabel[tmpDonor == trainDonor]

            testData <- tmpData[, tmpDonor != trainDonor]
            testLabel <- tmpLabel[tmpDonor != trainDonor]

            simData <- dataGenerate(trainData, trainLabel, testData, testLabel, n = nBulk)
        } else { #if there is only one donor, use 50% for reference matrix, skip for now
            next
            idx <- sample(seq(length(tmpLabel)), round(length(tmpLabel)/2))
            trainData <- tmpData[, idx]
            trainLabel <- tmpLabel[idx]

            testData <- tmpData[, -idx]
            testLabel <- tmpLabel[-idx]

            simData <- dataGenerate(trainData, trainLabel, testData, testLabel, n = nBulk)
        }
        tissue <- as.character(metaDat$tissue_in_publication[1])

        print(tissue)
        print(dim(simData$bulk))
        print(dim(simData$singleCellExpr))
        print(dim(simData$cellTypeExpr))
        print(dim(simData$signature))

        saveRDS(simData, paste0(savePath, tissue, ifelse(appendBulkCount, paste0("_", nBulk), ""), '.rds'))
    }
}

generateTSP(nBulk = 100, nGene = 1e4, savePath = './data/sim-tsp/')
generateTSP(nBulk = 100, nGene = 1e4, savePath = './data/scalability/', appendBulkCount = TRUE)
generateTSP(nBulk = 250, nGene = 1e4, savePath = './data/scalability/', appendBulkCount = TRUE)
generateTSP(nBulk = 500, nGene = 1e4, savePath = './data/scalability/', appendBulkCount = TRUE)
generateTSP(nBulk = 1000, nGene = 1e4, savePath = './data/scalability/', appendBulkCount = TRUE)