library(babelwhale)
library(rhdf5)
library(parallel)
library(RcppHungarian)

EX_SLOW_METHODS <- c("BayesCCE", "DAISM", "BayCount", "debCAM", "DeCompress", "BayesPrism")
SLOW_METHODS <- c("EMeth", "BayICE", "deconvSeq", "CPM", "DWLS", "DESeq2", "spatialDWLS", "deconf", "DecOT", "DeMixT", "digitalDLSorter", "ImmuCellAI", "Linseed", "MOMF", "AdRoit", "NITUMID", "quanTIseq", "scaden", "ARIC")
FAST_METHODS <- c("AutoGeneS", "BisqueMarker", "BisqueRef", "CellDistinguisher", "CIBERSORT", "Deblender", "DeconICA", "DeconPeaker", "DeconRNASeq", "DSA", "dtangle", "EPIC", "FARDEEP", "LinDeconSeq", "MCPcounter", "MethylResolver", "MIXTURE", "MuSic", "MySort", "PREDE", "ReFACTor", "RNA-Sieve", "SCDC", "TOAST", "BseqSC")
ALL_METHODS <- c(FAST_METHODS, SLOW_METHODS, EX_SLOW_METHODS)

PATH_TO_MATLAB_LICENSE <- "/path/to/matlab/license.lic"

DOCKER_TAG_PREFIX <- "deconvolution/"
DOCKER_TAG_SUFFIX <- ":latest"
nRep <- 10

dataFiles <- list.files('./data/sim-tsp') %>%
    str_match('([^/]+).rds$') %>%
    { function(x) x[, 2] }()
allConfig <- expand.grid(dataset = dataFiles, method = methods, seed = seq(nRep), stringsAsFactors = F)
# 8 cores, 16GB per job, 12GB softlimit, 12 concurrent jobs, max runtime of 12 hours

dir.create('./results/consistency', showWarnings = F, recursive = T)
res <- lapply(seq_len(nrow(allConfig)), function(i) {
    config <- allConfig[i,]
    dataFile <- paste0('./data/sim-tsp/', as.character(config$dataset), '.rds')
    method <- config$method
    seed <- config$seed

    resFile <- paste0('./results/consistency/', method, '_', config$dataset, '_', seed, '.rds')
    if (file.exists(resFile)) {
        res <- readRDS(resFile)
        if (!is.null(res$P)) {
            return(NULL)
        }
    }

    data <- readRDS(dataFile)
    groundTruth <- data$bulkRatio
    data <- data[names(data) != "bulkRatio"]

    #randomly remove markers and add noise
    set.seed(seed)
    data$markers <- lapply(data$markers, function(x) sample(x, size = round(length(x) * 0.95), replace = FALSE))
    data$sigGenes <- unique(unlist(data$markers))
    data$signature <- data$signature[data$sigGenes,]

    data$cellTypeExpr <- apply(data$cellTypeExpr, 1, function(x) {
        x <- x + rnorm(length(x), mean = 0, sd = sd(x) * 0.01)
        x[x < 0] <- 0
        x
    }) %>% t()

    data$signature <- apply(data$signature, 1, function(x) {
        x <- x + rnorm(length(x), mean = 0, sd = sd(x) * 0.01)
        x[x < 0] <- 0
        x
    }) %>% t()

    data$singleCellExpr <- apply(data$singleCellExpr, 1, function(x) {
        x <- x + rnorm(length(x), mean = 0, sd = sd(x) * 0.01)
        x[x < 0] <- 0
        round(x)
    }) %>% t()

    data$bulk <- apply(data$bulk, 1, function(x) {
        x <- x + rnorm(length(x), mean = 0, sd = sd(x) * 0.01)
        x[x < 0] <- 0
        x
    }) %>% t()

    set.seed(i)
    start <- Sys.time()
    res <- do.call(
        runDeconvolution,
        c(
            list(methods = method, verbose = T, seed = seed,
                 dockerArgs = c(
                     '--cpus=8.0',
                     '-m=16G',
                     '--memory-reservation=12G'
                 ),
                 timeout = 12 * 3600,
                 matlabLicenseFile = PATH_TO_MATLAB_LICENSE),
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

allRes <- lapply(seq_len(nrow(allConfig)), function(i) {
    config <- allConfig[i,]
    method <- config$method
    seed <- config$seed

    resFile <- paste0('./results/consistency/', method, '_', config$dataset, '_', seed, '.rds')
    dataFile <- paste0('./data/sim-tsp/', as.character(config$dataset), '.rds')

    if (!file.exists(resFile)) {
        return(NULL)
    }

    res <- readRDS(resFile)

    if (is.null(res$P)) return(NULL)

    data <- readRDS(dataFile)
    groundTruth <- data$bulkRatio

    res$MAE <- evaluateResult(res$P, groundTruth, metric="MAE")
    res$Corr <- evaluateResult(res$P, groundTruth, metric="Corr")
    data.frame(
        method = method,
        dataset = config$dataset,
        seed = seed,
        MAE = res$MAE,
        Corr = res$Corr,
        sample = 1:nrow(res$P)
    )
}) %>% do.call(what = rbind)

saveRDS(allRes, file = './results/consistency/allRes.rds')
