library(babelwhale)
library(rhdf5)
library(parallel)

source('./include/utils.R')


EX_SLOW_METHODS <- c("BayesCCE", "DAISM", "BayCount", "debCAM", "DeCompress", "BayesPrism")
SLOW_METHODS <- c("EMeth", "BayICE", "deconvSeq", "CPM", "DWLS", "DESeq2", "spatialDWLS", "deconf", "DecOT", "DeMixT", "digitalDLSorter", "ImmuCellAI", "Linseed", "MOMF", "AdRoit", "NITUMID", "quanTIseq", "scaden", "ARIC")
FAST_METHODS <- c("AutoGeneS", "BisqueMarker", "BisqueRef", "CellDistinguisher", "CIBERSORT", "Deblender", "DeconICA", "DeconPeaker", "DeconRNASeq", "DSA", "dtangle", "EPIC", "FARDEEP", "LinDeconSeq", "MCPcounter", "MethylResolver", "MIXTURE", "MuSic", "MySort", "PREDE", "ReFACTor", "RNA-Sieve", "SCDC", "TOAST", "BseqSC")
ALL_METHODS <- c(FAST_METHODS, SLOW_METHODS, EX_SLOW_METHODS)

DOCKER_TAG_PREFIX <- "deconvolution/"
DOCKER_TAG_SUFFIX <- ":latest"

PATH_TO_MATLAB_LICENSE <- "/path/to/matlab/license.lic"

nRep <- 10

methods <- ALL_METHODS

dataFiles <- list.files('./data/sim-tsp') %>%
    str_match('([^/]+).rds$') %>%
    { function(x) x[, 2] }()

allConfig <- expand.grid(dataFiles, methods, c(100, 250, 500, 1000), stringsAsFactors = FALSE)
allConfig <- allConfig[order(factor(allConfig$Var1, levels= dataFiles)),]
rownames(allConfig) <- NULL

dir.create('./results/scalability', showWarnings = FALSE, recursive = TRUE)
# 8 cores, 16GB per job, 12GB softlimit, 10 concurrent jobs, max runtime of 12 hours
res <- lapply(seq_len(nrow(allConfig)), function(i) {
  method <- as.character(allConfig[i,]$Var2)
  nBulk <- as.integer(allConfig[i,]$Var3)
  dataFile <- paste0('./data/scalability/', as.character(allConfig[i,]$Var1), '_', nBulk, '.rds')

  resFile <- paste0('./results/scalability/', method, '_', as.character(allConfig[i,]$Var1), '_', nBulk, '.rds')
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

  res
})

# gather results

methods        <- ALL_METHODS
tissues <- list.files("./data/sim-tsp") |>
    str_match('([^/]+).rds$') |>
    { \(x) x[, 2] }()

allConfig      <- expand.grid(tissues, methods, c(100, 250, 500, 1000))
f              <- "./results/scalability/mem-monitor.rds"
allRecord      <- readRDS(f)
methodValid    <- NULL
methodError    <- NULL
methodNotFound <- NULL
allRes         <- NULL
for (i in seq_len(nrow(allConfig))) {
  method   <- as.character(allConfig[i,]$Var2)
  nBulk    <- as.integer(allConfig[i,]$Var3)
  tissue   <- as.character(allConfig[i,]$Var1)
  dataFile <- paste0('./results/scalability/', method, '_', as.character(allConfig[i,]$Var1), '_', nBulk, '.rds')
  if (!file.exists(dataFile)) {
    methodNotFound <- c(methodNotFound, method)
    next
  }
  tryCatch(
  {
    res <- readRDS(dataFile)
  },
      error = function(e) {
        res <- NULL
      }
  )
  if (is.null(res)) {
    methodError <- c(methodError, method)
    next
  }

  if (!is.null(res$stderr))
  {
    allRes      <- rbind(allRes, data.frame(method = method, tissue = tissue, nBulk = nBulk, rt = NA, mem = NA))
    methodError <- c(methodError, method)
  } else {
    if (!rownames(res$groundTruth)[1] %in% colnames(res$P)) {
      if (ncol(res$P) < nrow(res$groundTruth))
      {
        allRes      <- rbind(allRes, data.frame(method = method, tissue = tissue, nBulk = nBulk, rt = NA, mem = NA))
        methodError <- c(methodError, method)
        next
      }
      if (is.null(colnames(res$P))) {
        colnames(res$P) <- paste0('CT', seq_len(ncol(res$P)))
      }
    }
    if (res$dockerName %in% names(allRecord)) {
      allRes <- rbind(allRes,
                      data.frame(method = method, tissue = tissue, nBulk = nBulk, rt = as.numeric(res$runningTime, units = 'secs'), mem = as.numeric(unlist(allRecord[res$dockerName]))))
    } else {
      allRes <- rbind(allRes,
                      data.frame(method = method, tissue = tissue, nBulk = nBulk, rt = as.numeric(res$runningTime, units = 'secs'), mem = 0.1))
    }

    methodValid <- c(methodValid, method)
  }
}

save(allRes, file = "./results/scalability/allRes.rds")


