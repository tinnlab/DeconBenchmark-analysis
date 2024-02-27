dataGenerate <- function(trainData, trainLabel, testData, testLabel, n = 10, minCell = 20, minMarkers = 10) {
  set.seed(1)
  trainLabel <- as.character(trainLabel)
  testLabel <- as.character(testLabel)

  # Choose 1000 random cells for each cell type
  idx <- NULL
  for (cellType in unique(trainLabel)) {
    if (length(which(trainLabel == cellType)) > 1000) {
      tmpIdx <- which(trainLabel == cellType)
      idx <- c(idx, sample(tmpIdx, 1000))
    } else {
      idx <- c(idx, which(trainLabel == cellType))
    }
  }
  trainData <- trainData[, idx]
  trainLabel <- trainLabel[idx]

  # Remove genes with only zeros in both train and test data
  idx <- which(rowSums(trainData) > 0 & rowSums(testData) > 0)
  trainData <- trainData[idx, ]
  testData <- testData[idx, ]

  set.seed(1)
  markers <- generateSignLimma(trainData, trainLabel, log2.threshold = 0.5)

  # Keep genes in marker list, still thousands
  trainData <- trainData[unique(markers$gene), ]
  testData <- testData[unique(markers$gene), ]

  # Keep top 200 markers for each cell type
  markers <- markers %>% dplyr::filter(abs(log2FC) >= 1)
  markers <- markers %>% dplyr::arrange(CT,desc(abs(log2FC)))
  markers <- markers %>% group_by(CT) %>% slice_head(n=200) %>% data.frame()

  # Keep cell types with more than minMarkers markers
  keepCellTypes <- table(markers$CT) %>% '>'(minMarkers) %>% names(.)[.]
  markers <- markers[markers$CT %in% keepCellTypes,]
  idx <- trainLabel %in% keepCellTypes
  trainData <- trainData[, idx]
  trainLabel <- trainLabel[idx]

  # Keep only valid cell types
  validCellType <- intersect(names(table(testLabel)[table(testLabel) > minCell]),
                             names(table(trainLabel)[table(trainLabel) > minCell]))
  idx <- which(trainLabel %in% validCellType)
  trainData <- trainData[,idx]
  trainLabel <- trainLabel[idx]

  idx <- which(testLabel %in% validCellType)
  testData <- testData[,idx]
  testLabel <- testLabel[idx]

  markers <- markers[markers$CT %in% validCellType, ]

  # Generate reference matrix
  refMatrix <- matrix(ncol = length(unique(trainLabel)), nrow = nrow(trainData))
  colnames(refMatrix) <- unique(trainLabel)
  rownames(refMatrix) <- rownames(trainData)
  for (cellType in unique(trainLabel)) {
    tmp <- rowMeans(trainData[,trainLabel==cellType])
    # tmp <- tmp/sum(tmp)*1e6
    refMatrix[,cellType] <- tmp
  }

  # Generate simulated data
  testSet <- list()
  targetN <- min(5000, length(testLabel)*5)
  cellTypes <- unique(testLabel)
  cellTypesRatios <- table(testLabel)/length(testLabel)
  for (i in seq(n)) {
    set.seed(i)
    nCellTypesTest <- seq(from=length(cellTypes)-(length(cellTypes)%/%5), to=length(cellTypes))
    if(length(nCellTypesTest) > 1) nCellTypesTest <- sample(nCellTypesTest, 1)
    cellTypesTest <- sample(cellTypes, nCellTypesTest)
    cellTypesTestRatios <- cellTypesRatios[cellTypesTest]/sum(cellTypesRatios[cellTypesTest]) # Normalize to 1
    cellTypesTestRatios <- cellTypesTestRatios^(1/2) # Reduce the extreme ratios
    cellTypesTestRatios <- cellTypesTestRatios*runif(nCellTypesTest, min = 0.5, max = 1.5) # Randomize the ratios
    cellTypesTestRatios <- cellTypesTestRatios/sum(cellTypesTestRatios) # Normalize to 1
    names(cellTypesTestRatios) <- cellTypesTest

    count <- 0
    cellTypeIdx <- NULL
    for (cellType in cellTypesTest) {
      cellTypeN <- round(targetN*cellTypesTestRatios[cellType])
      count <- count + cellTypeN
      tmp <- which(testLabel == cellType)
      tmp <- sample(tmp, round(length(tmp)/2))
      cellTypeIdx <- c(cellTypeIdx, sample(tmp, cellTypeN, replace = TRUE))
    }

    mixture <- rowMeans(testData[, cellTypeIdx])
    # mixture <- mixture/sum(mixture)*1e6
    testSet[[i]] <- list()
    ratio <- rep(0, length(cellTypes))
    names(ratio) <- cellTypes
    ratio[cellTypesTest] <- cellTypesTestRatios
    testSet[[i]]$ratio <- ratio
    testSet[[i]]$mixture <- mixture
  }

  testSetRatio <- do.call(cbind, lapply(testSet, function(x) x$ratio))
  testSet <- do.call(cbind, lapply(testSet, function(x) x$mixture))
  colnames(testSet) <- paste0("Mixture", seq_len(n))
  colnames(testSetRatio) <- paste0("Mixture", seq_len(n))

  # Keep valid genes
  validGenes <- Reduce(intersect, list(rownames(testSet[rowSums(testSet) > 0,]),
                                       rownames(trainData[rowSums(trainData) > 0,]),
                                       rownames(refMatrix[rowSums(refMatrix) > 0,])))
  testSet <- testSet[validGenes,]
  refMatrix <- refMatrix[validGenes,]
  trainData <- trainData[validGenes,]
  markers <- markers[markers$gene %in% validGenes,]

  tmp  <- lapply(unique(markers$CT), function (CT) markers$gene[markers$CT == CT])
  names(tmp) <- unique(markers$CT)
  print(sapply(tmp, function(x) length(x)))
  list(
    bulk = testSet,
    bulkRatio = testSetRatio,
    nCellTypes = length(unique(trainLabel)),
    markers = tmp,
    cellTypeExpr = refMatrix,
    sigGenes = unique(markers$gene),
    signature = refMatrix[unique(markers$gene),],
    singleCellExpr = trainData,
    singleCellLabels = trainLabel,
    singleCellSubjects = rep("Subject1", length(trainLabel))
  )
}

dataGenerateMissing <- function(trainData, trainLabel, testData, testLabel, n = 10, minCell = 20) {
  set.seed(1)
  trainLabel <- as.character(trainLabel)
  testLabel <- as.character(testLabel)

  # Choose 1000 random cells for each cell type
  idx <- NULL
  for (cellType in unique(trainLabel)) {
    if (length(which(trainLabel == cellType)) > 1000) {
      tmpIdx <- which(trainLabel == cellType)
      idx <- c(idx, sample(tmpIdx, 1000))
    } else {
      idx <- c(idx, which(trainLabel == cellType))
    }
  }
  trainData <- trainData[, idx]
  trainLabel <- trainLabel[idx]

  # Remove genes with only zeros in both train and test data
  idx <- which(rowSums(trainData) > 0 & rowSums(testData) > 0)
  trainData <- trainData[idx, ]
  testData <- testData[idx, ]

  set.seed(1)
  markers <- generateSignLimma(trainData, trainLabel, log2.threshold = 0.5)

  # Keep genes in marker list, still thousands
  trainData <- trainData[unique(markers$gene), ]
  testData <- testData[unique(markers$gene), ]

  # Keep top 200 markers for each cell type
  markers <- markers %>% dplyr::filter(abs(log2FC) >= 1)
  markers <- markers %>% dplyr::arrange(CT,desc(abs(log2FC)))
  markers <- markers %>% group_by(CT) %>% slice_head(n=200) %>% data.frame()
  # markers <- markers %>% group_by(CT) %>% data.frame()

  # Keep cell types with more than 10 markers
  # keepCellTypes <- table(markers$CT) %>% '>'(10) %>% names(.)[.]
  keepCellTypes <- table(markers$CT) %>% '>'(1) %>% names(.)[.]
  markers <- markers[markers$CT %in% keepCellTypes,]
  idx <- trainLabel %in% keepCellTypes
  trainData <- trainData[, idx]
  trainLabel <- trainLabel[idx]

  # Keep only valid cell types
  validCellType <- intersect(names(table(testLabel)[table(testLabel) > minCell]),
                             names(table(trainLabel)[table(trainLabel) > minCell]))
  idx <- which(trainLabel %in% validCellType)
  trainData <- trainData[,idx]
  trainLabel <- trainLabel[idx]

  # idx <- which(testLabel %in% validCellType)
  # testData <- testData[,idx]
  # testLabel <- testLabel[idx]

  markers <- markers[markers$CT %in% validCellType, ]

  # Generate reference matrix
  refMatrix <- matrix(ncol = length(unique(trainLabel)), nrow = nrow(trainData))
  colnames(refMatrix) <- unique(trainLabel)
  rownames(refMatrix) <- rownames(trainData)
  for (cellType in unique(trainLabel)) {
    tmp <- rowMeans(trainData[,trainLabel==cellType])
    # tmp <- tmp/sum(tmp)*1e6
    refMatrix[,cellType] <- tmp
  }

  # Generate simulated data
  testSet <- list()
  targetN <- min(5000, length(testLabel)*5)
  cellTypes <- unique(testLabel)
  cellTypesRatios <- table(testLabel)/length(testLabel)
  for (i in seq(n)) {
    set.seed(i)
    nCellTypesTest <- seq(from=length(cellTypes)-(length(cellTypes)%/%5), to=length(cellTypes))
    if(length(nCellTypesTest) > 1) nCellTypesTest <- sample(nCellTypesTest, 1)
    cellTypesTest <- sample(cellTypes, nCellTypesTest)
    cellTypesTestRatios <- cellTypesRatios[cellTypesTest]/sum(cellTypesRatios[cellTypesTest]) # Normalize to 1
    cellTypesTestRatios <- cellTypesTestRatios^(1/2) # Reduce the extreme ratios
    cellTypesTestRatios <- cellTypesTestRatios*runif(nCellTypesTest, min = 0.5, max = 1.5) # Randomize the ratios
    cellTypesTestRatios <- cellTypesTestRatios/sum(cellTypesTestRatios) # Normalize to 1
    names(cellTypesTestRatios) <- cellTypesTest

    count <- 0
    cellTypeIdx <- NULL
    for (cellType in cellTypesTest) {
      cellTypeN <- round(targetN*cellTypesTestRatios[cellType])
      count <- count + cellTypeN
      tmp <- which(testLabel == cellType)
      tmp <- sample(tmp, round(length(tmp)/2))
      cellTypeIdx <- c(cellTypeIdx, sample(tmp, cellTypeN, replace = TRUE))
    }

    mixture <- rowMeans(testData[, cellTypeIdx])
    # mixture <- mixture/sum(mixture)*1e6
    testSet[[i]] <- list()
    ratio <- rep(0, length(cellTypes))
    names(ratio) <- cellTypes
    ratio[cellTypesTest] <- cellTypesTestRatios
    testSet[[i]]$ratio <- ratio
    testSet[[i]]$mixture <- mixture
  }

  testSetRatio <- do.call(cbind, lapply(testSet, function(x) x$ratio))
  testSet <- do.call(cbind, lapply(testSet, function(x) x$mixture))
  colnames(testSet) <- paste0("Mixture", seq_len(n))
  colnames(testSetRatio) <- paste0("Mixture", seq_len(n))

  # Keep valid genes
  validGenes <- Reduce(intersect, list(rownames(testSet[rowSums(testSet) > 0,]),
                                       rownames(trainData[rowSums(trainData) > 0,]),
                                       rownames(refMatrix[rowSums(refMatrix) > 0,])))
  testSet <- testSet[validGenes,]
  refMatrix <- refMatrix[validGenes,]
  trainData <- trainData[validGenes,]
  markers <- markers[markers$gene %in% validGenes,]

  tmp  <- lapply(unique(markers$CT), function (CT) markers$gene[markers$CT == CT])
  names(tmp) <- unique(markers$CT)
  print(sapply(tmp, function(x) length(x)))
  list(
    bulk = testSet,
    bulkRatio = testSetRatio,
    nCellTypes = length(unique(trainLabel)),
    markers = tmp,
    cellTypeExpr = refMatrix,
    sigGenes = unique(markers$gene),
    signature = refMatrix[unique(markers$gene),],
    singleCellExpr = trainData,
    singleCellLabels = trainLabel,
    singleCellSubjects = rep("Subject1", length(trainLabel))
  )
}



library(tidyverse)
generateSignLimma <- function(data, label, log2.threshold = 1) {
  colnames(data) <- as.character(label)
  train <- data

  #for marker selection, keep genes where at least 30% of cells within a cell type have a read/UMI count different from 0
  cellType <- colnames(train)
  keep <- sapply(unique(cellType), function(x) {
    CT_hits <- which(cellType %in% x)
    size <- ceiling(0.3*length(CT_hits))
    Matrix::rowSums(train[,CT_hits,drop=FALSE] != 0) >= size
  })
  train <- train[Matrix::rowSums(keep) > 0,]
  train2 <- Normalization(train)

  # INITIAL CONTRASTS for marker selection WITHOUT taking correlated CT into account
  #[compare one group with average expression of all other groups]
  annotation <- factor(colnames(train2))
  design <- model.matrix(~0+annotation)
  colnames(design) <- unlist(lapply(strsplit(colnames(design),"annotation"), function(x) x[2]))
  cont.matrix <- matrix((-1/ncol(design)),nrow=ncol(design),ncol=ncol(design))
  colnames(cont.matrix) <- colnames(design)
  diag(cont.matrix) <- (ncol(design)-1)/ncol(design)

  v <- limma::voom(train2, design=design, plot=FALSE)
  fit <- limma::lmFit(v, design)
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2, trend=TRUE)

  markers <- marker.fc(fit2, cont.matrix, log2.threshold = log2.threshold)
  markers
}

Normalization <- function(data){

  data <- edgeR::DGEList(data)

  data <- edgeR::calcNormFactors(data, method = "TMM")

  return(data)

}

marker.fc <- function(fit2, cont.matrix, log2.threshold = 1, output_name = "markers"){

  topTable_RESULTS <- limma::topTable(fit2, coef = 1:ncol(cont.matrix), number = Inf, adjust.method = "BH", p.value = 0.05, lfc = log2.threshold)
  AveExpr_pval <- topTable_RESULTS[,(ncol(topTable_RESULTS)-3):ncol(topTable_RESULTS)]
  topTable_RESULTS <- topTable_RESULTS[,1:(ncol(topTable_RESULTS)-4)]

  if(length(grep("ERCC-",topTable_RESULTS$gene)) > 0){ topTable_RESULTS <- topTable_RESULTS[-grep("ERCC-",topTable_RESULTS$gene),] }

  markers <- apply(topTable_RESULTS,1,function(x){
    temp <- sort(x)
    ((temp[ncol(topTable_RESULTS)] - temp[ncol(topTable_RESULTS)-1]) >= log2.threshold) | (abs(temp[1] - temp[2]) >= log2.threshold)

  })

  topTable_RESULTS <- topTable_RESULTS[markers,]

  markers <- cbind.data.frame(rownames(topTable_RESULTS),
                              t(apply(topTable_RESULTS, 1, function(x){
                                temp <- max(x)
                                if(temp < log2.threshold){
                                  temp <- c(min(x), colnames(topTable_RESULTS)[which.min(x)])
                                } else {
                                  temp <- c(max(x), colnames(topTable_RESULTS)[which.max(x)])
                                }
                                temp
                              })))

  colnames(markers) <- c("gene","log2FC","CT")
  markers$log2FC <- as.numeric(as.character(markers$log2FC))
  markers <- markers %>% dplyr::arrange(CT,desc(log2FC))

  markers$AveExpr <- AveExpr_pval$AveExpr[match(markers$gene,rownames(AveExpr_pval))]
  markers$gene <- as.character(markers$gene)
  markers$CT <- as.character(markers$CT)

  #write.table(markers, file = output_name, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

  return(markers)

}

writeArgs <- function(h5file,
                      bulk,
                      nCellTypes = NULL,
                      markers = NULL,
                      isMethylation = F,
                      seed = 1,
                      singleCellExpr = NULL,
                      singleCellLabels = NULL,
                      singleCellSubjects = NULL,
                      cellTypeExpr = NULL,
                      sigGenes = NULL,
                      signature = NULL
) {

  writeListOrVector <- function(v, name) {
    h5createGroup(h5file, name)
    h5write(v, h5file, paste0(name, "/values"))
    if (!is.null(names(v))) {
      h5write(names(v), h5file, paste0(name, "/names"))
    }
  }

  writeMatrix <- function(m, name) {
    h5createGroup(h5file, name)
    h5write(m, h5file, paste0(name, "/values"))
    if (!is.null(rownames(m))) {
      h5write(rownames(m), h5file, paste0(name, "/rownames"))
    }
    if (!is.null(colnames(m))) {
      h5write(colnames(m), h5file, paste0(name, "/colnames"))
    }
  }

  unlink(h5file)
  h5createFile(h5file)

  if (!is.null(seed))               writeListOrVector(seed, "seed")
  if (!is.null(bulk))               writeMatrix(bulk, "bulk")
  if (!is.null(nCellTypes))         writeListOrVector(nCellTypes, "nCellTypes")
  if (!is.null(markers))            writeListOrVector(markers, "markers")
  if (!is.null(cellTypeExpr))      writeMatrix(cellTypeExpr, "cellTypeExpr")
  if (!is.null(sigGenes))           writeListOrVector(sigGenes, "sigGenes")
  if (!is.null(signature))          writeMatrix(signature, "signature")
  if (!is.null(singleCellExpr))     writeMatrix(singleCellExpr, "singleCellExpr")
  if (!is.null(singleCellLabels))   writeListOrVector(singleCellLabels, "singleCellLabels")
  if (!is.null(singleCellSubjects)) writeListOrVector(singleCellSubjects, "singleCellSubjects")
  if (!is.null(isMethylation))      writeListOrVector(isMethylation, "isMethylation")
}

# runDeconvolution <- function(methods,
#                              bulk,
#                              nCellTypes = NULL,
#                              markers = NULL,
#                              isMethylation = F,
#                              singleCellExpr = NULL,
#                              singleCellLabels = NULL,
#                              singleCellSubjects = NULL,
#                              cellTypeExpr = NULL,
#                              sigGenes = NULL,
#                              signature = NULL,
#                              seed = 1,
#                              matlabLicenseFile = NULL,
#                              timeout = 60*60*12,
#                              containerEngine = "docker",
#                              dockerArgs = c(
#                                  '--cpus=8.0',
#                                  '-m=32G',
#                                  '--memory-reservation=4G'
#                              ),
#                              verbose = T) {
#
#   tmpDir <- file.path(tempdir(), paste0(sample(c(LETTERS, letters), 10, TRUE), collapse = ""))
#   dir.create(tmpDir, recursive = TRUE, showWarnings = F)
#   tmpH5File <- file.path(tmpDir, "args.h5")
#
#   writeArgs(h5file = tmpH5File,
#             bulk = bulk,
#             nCellTypes = nCellTypes,
#             markers = markers,
#             isMethylation = isMethylation,
#             singleCellExpr = singleCellExpr,
#             singleCellLabels = singleCellLabels,
#             singleCellSubjects = singleCellSubjects,
#             cellTypeExpr = cellTypeExpr,
#             sigGenes = sigGenes,
#             signature = signature,
#             seed = seed)
#
#   .message <- function(...) {
#     if (verbose) message(paste(...))
#   }
#
#   allResults <- list()
#
#   for (method in methods) {
#     if (containerEngine == "docker"){
#       dockerName <- paste0(sample(c(LETTERS, letters), 20, TRUE), collapse = "")
#       params <- c("run", "--rm", dockerArgs,  paste0('--name=', dockerName, '_', method),
#                   "-v", paste0(tmpH5File, ":", "/input/args.h5"),
#                   "-v", paste0(tmpDir, ":", "/output"),
#                   '-e', 'INPUT_PATH=/input/args.h5',
#                   '-e', paste0('OUTPUT_PATH=/output/', method,'-results.h5')
#       )
#
#       if (!is.null(matlabLicenseFile)) {
#         params <- c(params, "-v", paste0(matlabLicenseFile, ':', '/licenses/license.lic'), '--network=host')
#       }
#       params <- c(params, paste0(DOCKER_TAG_PREFIX, tolower(method), DOCKER_TAG_SUFFIX))
#
#       .message("Running docker ", paste(params, collapse = " "))
#
#       res <- tryCatch(
#           processx::run("docker",
#                         params,
#                         timeout = timeout,
#                         error_on_status = FALSE,
#                         stdout_callback = function(newout, proc) {
#                           .message(newout)
#                         },
#                         stderr_callback = function(newerr, proc) {
#                           message(newerr)
#                         }),
#           error = function(e) {
#             list(status = T)
#           }
#       )
#     } else if (containerEngine == "singularity"){
#       params <- c("run",
#                   '--env', paste0('INPUT_PATH=', tmpH5File),
#                   '--env', paste0('OUTPUT_PATH=', file.path(tmpDir, paste0(method, "-results.h5")))
#       )
#
#       if (!is.null(matlabLicenseFile)) {
#         params <- c(params, "--env", paste0("MLM_LICENSE_FILE=", matlabLicenseFile))
#       }
#       params <- c(params, paste0("docker://", .getMethodDockerRepos(method)))
#
#       .message("Running singularity ", paste(params, collapse = " "))
#
#       res <- tryCatch(
#           processx::run("singularity",
#                         params,
#                         timeout = timeout,
#                         error_on_status = FALSE,
#                         stdout_callback = function(newout, proc) {
#                           .message(newout)
#                         },
#                         stderr_callback = function(newerr, proc) {
#                           message(newerr)
#                         }),
#           error = function(e) {
#             message(e)
#             list(status = T)
#           }
#       )
#     } else {
#       stop("Unsupported container engine", containerEngine)
#     }
#
#     if (res$status) {
#       if (res$timeout) {
#         params <- c('rm', '-f', paste0(dockerName, '_', method))
#         processx::run("docker",
#                       params)
#       }
#
#       allResults[[method]] <- res
#       next()
#     }
#
#     resFile <- file.path(tmpDir, paste0(method, "-results.h5"))
#     if (length(resFile) == 0) {
#       allResults[[method]] <- NULL
#       next()
#     }
#
#     P <- tryCatch(
#         rhdf5::h5read(resFile, "P"),
#         error = function(e) {
#           # message(e)
#           NULL
#         }
#     )
#
#     S <- tryCatch(
#         rhdf5::h5read(resFile, "S"),
#         error = function(e) {
#           # message(e)
#           NULL
#         }
#     )
#
#     if (!is.null(P)) {
#       .P <- P$values
#       if (!is.null(P$rownames)) {
#         rownames(.P) <- P$rownames
#       }
#       if (!is.null(P$colnames)) {
#         colnames(.P) <- P$colnames
#       }
#       P <- .P
#     }
#
#     if (!is.null(S)) {
#       .S <- S$values
#       if (!is.null(S$rownames)) {
#         rownames(.S) <- S$rownames
#       }
#       if (!is.null(S$colnames)) {
#         colnames(.S) <- S$colnames
#       }
#       S <- .S
#     }
#
#     allResults[[method]] <- list(P = P, S = S, dockerName = paste0(dockerName, '_', method))
#   }
#
#   unlink(tmpDir, recursive = TRUE)
#
#   allResults
# }

detectCellType <- function(P, groundTruth) {
  tmp              <- suppressWarnings(cor(P, t(groundTruth)))
  tmp[is.na(tmp)]  <- -1
  tmp              <- 1 - tmp
  mapRes           <- data.frame(HungarianSolver(tmp)$pairs)
  mapRes[, 1]      <- rownames(tmp)[mapRes[, 1]]
  mapRes[, 2]      <- colnames(tmp)[mapRes[, 2]]
  colnames(mapRes) <- c("from", "to")
  rownames(mapRes) <- mapRes[, 1]
  return(mapRes)
}

evaluateResult <- function(P, groundTruth, method = "sample", metric = "Corr") {
  P <- apply(P, 1, function(x) {
    x <- x - min(x)
    x <- x / sum(x)
    x
  }) # Normalize to sum to 1
  P <- t(P)

  if (!rownames(groundTruth)[1] %in% colnames(P)) {
    mapRes      <- detectCellType(P, groundTruth)
    colnames(P) <- mapRes[colnames(P),]$to
  }

  P <- P[, rownames(groundTruth)]
  if (metric == "Corr") {
    if (method == "sample") {
      tmp             <- sapply(seq_len(ncol(groundTruth)), function(i) {
        suppressWarnings(cor(unlist(P[i,]), unlist(groundTruth[, i]), method = "spearman"))
      })
      tmp[is.na(tmp)] <- 0
      return(tmp)
    }
    if (method == "CT") {
      tmp             <- sapply(rownames(groundTruth), function(i) {
        suppressWarnings(cor(unlist(P[, i]), unlist(groundTruth[i,]), method = "spearman"))
      })
      tmp[is.na(tmp)] <- 0
      return(tmp)
    }
  }
  if (metric == "MAE") {
    tmp             <- sapply(seq_len(ncol(groundTruth)), function(i) {
      suppressWarnings(mean(abs(unlist(P[i,]) - unlist(groundTruth[, i]))))
    })
    tmp[is.na(tmp)] <- 1 / nrow(groundTruth)
    return(tmp)
  }

  if (metric == "samplePairwise") {
    corrP  <- cor(t(P), method = "pearson")
    corrGT <- cor(groundTruth, method = "pearson")
    return(mean(abs(corrGT - corrP), na.rm = TRUE))
  }

  return(NA)  # Unknown method
}
