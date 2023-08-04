###########################################################
###### Gets alignment site patterns, sampling times #######
###########################################################

library(ape)
library(phangorn)

aln_files <- dir(path = ".", pattern = ".+[.]fasta")

data <- data.frame()
dates <- list()
for (file in aln_files) {

    aln <- read.dna(file, format = "fasta")

    data <- rbind(
        data,
        c(
            attr(phyDat(aln, return.index = TRUE, type = "DNA"), "nr"),
            length(rownames(aln))
        )
    )

    dates <- c(
        dates,
        as.Date(gsub(patterh = ".+_", replacement = "", rownames(x)))
    )

    print(paste(file, " Done"))

}

names(data) <- aln_files
names(dates) <- aln_files

save(data, dates, file = "alignment_data.RData")