#################
### This code is for generating simulated data using singcle cell data from CELLxGENE

library(Seurat)
library(tidyverse)

source('./include/utils.R')

metaFiles <- readRDS("./include/metaFiles.rds")

file_ids <- metaFiles$dataset_id[metaFiles$dataset == "CELLxGENE"]

generateCXGMissing <- function(nBulk, nGene, savePath, appendBulkCount = FALSE, percentage_missing = 0.75){
  dir.create(savePath, showWarnings = FALSE, recursive = TRUE)

  for (file_id in file_ids) {
    file_path <- paste0("./data/Upload/", file_id, ".rds")
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

    # Remove Cell types:
    df <- lapply(unique(as.character(metadata$donor_id)), function(did) {
      metaTMP <- metadata[metadata$donor_id == did,]
      ctCount <- table(metaTMP$cell_type)
      data.frame(
        donor_id = rep(did, length(unique(metaTMP$cell_type))),
        cell_type = names(ctCount),
        count = as.numeric(ctCount),
        stringsAsFactors = FALSE
      )
    }) %>% do.call(what = bind_rows)

    # 1. Select donor that have a least of `percentage` of cell types
    allCellTypes <- unique(df$cell_type)
    donorCellTypes <- df %>% group_by(donor_id) %>%
      group_split() %>% lapply(., function(xdf) {
      list(
        donor_id = as.character(xdf$donor_id[1]),
        cell_type = as.character(unique(xdf$cell_type))
      )
    })

    toUseDonors <- lapply(donorCellTypes, function(xls) {
      commonCellTypes <- intersect(xls$cell_type, allCellTypes)
      if (length(commonCellTypes) == length(allCellTypes)) {
        xls$donor_id
      } else {
        NULL
      }
    }) %>% unlist()

    # 2. Add more donor until we have at least 75% of cell types all donors
    calCellTypes <-  lapply(donorCellTypes, function(xls) {
      if (xls$donor_id %in% toUseDonors) {
        NULL
      } else {
        retx <- c(length(xls$cell_type))
        names(retx) <- xls$donor_id
        retx
      }
    }) %>% unlist() %>% sort(decreasing=TRUE)

    toCheckDonorCTls <- allCellTypes
    for (nameDN in names(calCellTypes)) {
      toCheckDonorCT <- donorCellTypes[lapply(donorCellTypes, function(xls) {xls$donor_id == nameDN}) %>% unlist()]
      toCheckDonorCTls <- intersect(toCheckDonorCTls, toCheckDonorCT[[1]]$cell_type)
      if ((length(toCheckDonorCTls) / length(allCellTypes)) >= percentage_missing) {
        toUseDonors <- c(toUseDonors, nameDN)
      } else {
        break
      }
    }

    # 3. Intersect cell type
    toIntersectCT <- lapply(donorCellTypes, function(xls) {
      if (xls$donor_id %in% toUseDonors) {
        xls$cell_type
      } else {
        NULL
      }
    }) %>% .[lengths(.) > 0] %>% Reduce(f = intersect)

    # 4. Calculate number of cells in each cell type
    callCellsCT <- df[df$donor_id %in% toUseDonors,] %>%
      dplyr::filter(cell_type %in% toIntersectCT) %>%
      group_by(cell_type) %>%
      summarise(count = sum(count)) %>%
      arrange(desc(count))
    toUseCellTypes <- as.character(callCellsCT$cell_type[1:round(18*percentage_missing)])


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

      # Use only toUseCellTypes
      toUseidx <- which(trainLabel %in% toUseCellTypes)
      trainData <- trainData[,toUseidx]
      trainLabel <- trainLabel[toUseidx]


      testData <- tmpData[, tmpDonor != trainDonor]
      testLabel <- tmpLabel[tmpDonor != trainDonor]

    simData <- dataGenerateMissing(trainData, trainLabel, testData, testLabel, n = nBulk, minCell = 10)

    print(tissue_in_publication)
    print(dim(simData$bulk))
    print(dim(simData$singleCellExpr))
    print(dim(simData$cellTypeExpr))
    print(dim(simData$signature))

    saveRDS(simData, paste0(savePath, tissue_in_publication, ifelse(appendBulkCount, paste0("_", nBulk), ""), '.rds'))

  }
}

generateCXGMissing(nBulk = 100, nGene = 1e4, percentage_missing = 0.75, savePath = './data/sim-cxg-missing25/')
generateCXGMissing(nBulk = 100, nGene = 1e4, percentage_missing = 0.50, savePath = './data/sim-cxg-missing50/')