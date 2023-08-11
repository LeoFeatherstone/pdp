###########################################################
####### Gets posteriors from logs files on cluster ########
###########################################################

log <- dir(path = ".", pattern = ".+[.]log")
log <- log[which(file.size(log) > 0)]

data <- lapply(
    log,
    function(x) {
        print(paste0("Read ", x))
        return(read.table(paste0(x), header = TRUE))
    }
)
names(data) <- gsub(
    log,
    pattern = "[.]log",
    replacement = ""
)

save(data, file = "posteriors.RData")
