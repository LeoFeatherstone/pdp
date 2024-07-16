#########################################
####### Find Sampling Spans #############
#########################################

library(ape)
library(tidyverse)

files <- list.files(
    path = "analysis/empirical_data",
    pattern = "(h1n1|sars-cov-2)[.]fasta",
    full.names = TRUE
)

aln <- lapply(
    files,
    function(x) read.dna(x, format = "fasta")
)

names(aln) <- files

date_range <- lapply(
    aln,
    function(x) {
        range(as.Date(gsub(pattern = ".+_", replacement = "", rownames(x))))
    }
)

print(date_range)
