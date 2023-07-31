###########################################################
###### Gets alignment site patterns, sampling times #######
###########################################################

library(ape)
library(phangorn)

aln_files <- dir(path = ".", pattern = ".+[.]fasta")
aln_files <- log[which(file.size(aln_files) > 0)]

aln <- lapply(
    aln_files,
    function(x) {
        read.dna(x, format = "fasta")
    }
)
print("Loaded alignments")
data <- lapply(
    aln,
    function(x) {
        c(
            attr(phangorn::phyDat(x, return.index = TRUE, type = "DNA"), "nr"),
            length(rownames(x))
        )
    }
)
print("Calculated Site Patterns")
dates <- lapply(
    aln,
    function(x) {
        return(
            as.Date(gsub(patterh = ".+_", replacement = "", rownames(x)))
        )
    }
)
print("Got Sampling Times")

names(data) <- aln_files
names(dates) <- aln_files

save(data, dates, file = "alignment_data.RData")