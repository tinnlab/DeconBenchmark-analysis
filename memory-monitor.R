library(processx)
library(stringr)

getRecord <- function(data, allRecord) {
    data <- strsplit(data, "\n")[[1]]
    data <- lapply(data, function(x) {
        tmp <- strsplit(x, "\t")[[1]]
        tmp[2] <- as.numeric(str_replace_all(tmp[2], "%", ""))
        tmp
    } )
    tmp <- lapply(data, function(x) as.numeric(x[2]))
    names(tmp) <- lapply(data, function(x) x[1])
    data <- tmp
    for (id in names(data)) {
        if(id %in% allRecord) {
            if(data[id] > allRecord[id]) {
                allRecord[id] <- data[id]
            }
        } else {
            allRecord[id] <- data[id]
        }
    }
    allRecord
}

f <- "./results/scalability/mem-monitor.rds"

if(file.exists(f)) {
    allRecord <- readRDS(f)
} else {
    allRecord <- list()
}
lastTime <- Sys.time()
while(TRUE)
{
    data <- tryCatch(processx::run(command = "docker",
                                   args = c('stats', '--no-stream', '--format', '{{.Name}}\t{{.MemPerc}}')),
                     error = function(e) {
                         list(status = T)
                     })
    if (data$status) {
        next()
    }
    data <- data$stdout
    tmpRecord <- try(getRecord(data, allRecord))
    if(is(tmpRecord, "try-error")) {
        next()
    }
    allRecord <- tmpRecord
    if(Sys.time() - lastTime > 30) {
        lastTime <- Sys.time()
        saveRDS(allRecord, f)
    }
}

saveRDS(allRecord, f)



