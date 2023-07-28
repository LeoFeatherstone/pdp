###########################################################
####### Gets posteriors from logs files on cluster ########
###########################################################

log <- dir(path = ".", pattern = ".+[.]log")

data <- lapply(
    sim_log,
    function(x) {
        print(paste0("Read ", x))
        return(read.table(paste0(x), header = TRUE))
    }
)
names(data) <- gsub(
    sim_log,
    pattern = "[.]log",
    replacement = ""
)

save(data, file = "posteriors.RData")