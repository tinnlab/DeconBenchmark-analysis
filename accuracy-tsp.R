library(babelwhale)
library(rhdf5)
library(parallel)
library(RcppHungarian)
library(stringr)
library(tidyverse)

# Load DeconBenchmark package
require(DeconBenchmark)
# If not install, run the following lines
# if (!requireNamespace("devtools", quietly = TRUE)) {
#    install.packages("devtools")
# }
# devtools::install_github("tinnlab/DeconBenchmark")

## Utilities functions
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


EX_SLOW_METHODS <- c("BayesCCE", "DAISM", "BayCount", "debCAM", "DeCompress", "BayesPrism")
SLOW_METHODS <- c("EMeth", "BayICE", "deconvSeq", "CPM", "DWLS", "DESeq2", "spatialDWLS", "deconf", "DecOT", "DeMixT", "digitalDLSorter", "ImmuCellAI", "Linseed", "MOMF", "AdRoit", "NITUMID", "quanTIseq", "scaden", "ARIC")
FAST_METHODS <- c("AutoGeneS", "BisqueMarker", "BisqueRef", "CellDistinguisher", "CIBERSORT", "Deblender", "DeconICA", "DeconPeaker", "DeconRNASeq", "DSA", "dtangle", "EPIC", "FARDEEP", "LinDeconSeq", "MCPcounter", "MethylResolver", "MIXTURE", "MuSic", "MySort", "PREDE", "ReFACTor", "RNA-Sieve", "SCDC", "TOAST", "BseqSC")
ALL_METHODS <- c(FAST_METHODS, SLOW_METHODS, EX_SLOW_METHODS)

PATH_TO_MATLAB_LICENSE <- "/path/to/matlab/license.lic"

DOCKER_TAG_PREFIX <- "deconvolution/"
DOCKER_TAG_SUFFIX <- ":latest"

methods <- ALL_METHODS

dataFiles <- list.files('./data/sim-tsp')
allConfig <- expand.grid(dataFiles, methods[2])

# 8 cores, 16GB per job, 12GB softlimit, 12 concurrent jobs, max runtime of 12 hours

dir.create('./results/accuracy-tsp', showWarnings = FALSE, recursive = TRUE)
lapply(seq_len(nrow(allConfig)), function(i) {
    dataFile <- paste0('./data/sim-tsp/', as.character(allConfig[i,]$Var1))
    method <- as.character(allConfig[i,]$Var2)

    resFile <- paste0('./results/accuracy-tsp/', method, '_', as.character(allConfig[i,]$Var1))
    if (file.exists(resFile)) {
        return(NULL)
    }

    data <- readRDS(dataFile)
    groundTruth <- data$bulkRatio
    data <- data[names(data) != "bulkRatio"]

    start <- Sys.time()
    res <- do.call(
        runDeconvolution,
        c(
            list(methods = method, verbose = TRUE,
                 dockerArgs = c(
                     '--cpus=8.0',
                     '-m=16G',
                     '--memory-reservation=12G'
                 ),
                 timeout = 12*3600,
                 matlabLicenseFile=PATH_TO_MATLAB_LICENSE),
            data
        )
    )
    runningTime <- Sys.time() - start

    res <- res[[method]]
    res$groundTruth <- groundTruth
    res$runningTime <- runningTime

    saveRDS(res, file = resFile)

    NULL
})

# gather results

dataPath <- "./data/sim-tsp"
datasets <- list.files(dataPath) %>%
    str_match('([^/]+).rds$') %>%
    { function(x) x[, 2] }()

methods <- ALL_METHODS
rawRes <- lapply(datasets, function(dataset) {
    groundTruthFile <- paste0('./data/sim-tsp/', dataset, '.rds')
    groundTruth     <- readRDS(groundTruthFile)$bulkRatio %>% t()

    res <- lapply(methods, function(method) {
        # message(dataset, method)
        resFile <- paste0('./results/accuracy-tsp/', method, '_', dataset, '.rds')
        if (!file.exists(resFile)) return(NULL)
        print(resFile)

        res <- readRDS(resFile)

        if (is.null(res$P)) return(NULL)

        if (!is.null(res$stderr)) return(NULL)

        if (!rownames(res$groundTruth)[1] %in% colnames(res$P)) {
            if (ncol(res$P) < nrow(res$groundTruth))
            {
                return(NULL)
            }
            if (is.null(colnames(res$P))) {
                colnames(res$P) <- paste0('CT', seq_len(ncol(res$P)))
            }
        }

        P <- res$P
        P <- apply(P, 1, function(x) {
            x <- x - min(x)
            x <- x / sum(x)
            x
        }) %>% t()

        if (!colnames(groundTruth)[1] %in% colnames(P)) {
            mapRes      <- detectCellType(P, t(groundTruth))
            colnames(P) <- mapRes[colnames(P),]$to
        }

        list(
            P           = P[, colnames(groundTruth)] %>%
                data.frame(
                    method  = method,
                    dataset = dataset,
                    sample  = rownames(P)
                ) %>%
                tidyr::gather("cellType", "value", -method, -dataset, -sample),
            groundTruth = groundTruth %>%
                data.frame(
                    method  = method,
                    dataset = dataset,
                    sample  = rownames(groundTruth)
                ) %>%
                tidyr::gather("cellType", "value", -method, -dataset, -sample)
        )
    })

    list(
        P           = lapply(res, function(x) x$P) %>% do.call(what = rbind),
        groundTruth = lapply(res, function(x) x$groundTruth) %>% do.call(what = rbind)
    )
}) %>% { function(res) {
    list(
        P           = lapply(res, function(x) x$P) %>% do.call(what = rbind),
        groundTruth = lapply(res, function(x) x$groundTruth) %>% do.call(what = rbind)
    )
} }()

allRes <- lapply(datasets, function(dataset) {
    groundTruthFile <- paste0('./data/sim-tsp/', dataset, '.rds')
    groundTruth     <- readRDS(groundTruthFile)$bulkRatio %>% t()

    # method <- methods[1]
    lapply(methods, function(method) {
        resFile <- paste0('./results/accuracy-tsp/', method, '_', dataset, '.rds')
        if (!file.exists(resFile)) return(NULL)

        res <- readRDS(resFile)
      
        if (is.null(res$P)) return(NULL)

        if (!is.null(res$stderr)) return(NULL)

        if (!rownames(res$groundTruth)[1] %in% colnames(res$P)) {
            if (ncol(res$P) < nrow(res$groundTruth))
            {
                return(NULL)
            }
            if (is.null(colnames(res$P))) {
                colnames(res$P) <- paste0('CT', seq_len(ncol(res$P)))
            }
        }

        SCorr <- mean(evaluateResult(res$P, res$groundTruth, method = "sample"), na.rm = TRUE)
        CCorr <- mean(evaluateResult(res$P, res$groundTruth, method = "CT"), na.rm = TRUE)
        MAE   <- mean(evaluateResult(res$P, res$groundTruth, metric = "MAE"), na.rm = TRUE)
        MAECorr   <- evaluateResult(res$P, res$groundTruth, metric = "samplePairwise")

        data.frame(
            dataset = dataset,
            method  = method,
            SCorr   = SCorr,
            CCorr   = CCorr,
            MAE     = MAE,
            MAECorr     = MAECorr
        )
    }) %>% do.call(what = rbind)

}) %>% do.call(what = rbind)

saveRDS(allRes, file = './results/accuracy-tsp/allRes.rds')