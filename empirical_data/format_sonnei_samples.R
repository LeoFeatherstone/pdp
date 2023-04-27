## One-off formatting sonnei samples from source
library(ape)
library(tidyverse)

aln <- read.dna(
    "empirical_data/Sson_NC_007384_BAPS3_CIDall.fasta",
    format = "fasta"
)

metadata <- read.csv(
    "empirical_data/ShigSon_BAPS3_08042023_dates_genotypes.csv",
    header = TRUE)

# subset 3.7.29.1.2.1 , format tips
metadata <- metadata %>%
    filter(final.genotype == "3.7.29.1.2.1") %>%
    mutate(new_name = paste0(
        Isolate,
        "_",
        year, "-01-01_",
        year, "-", sprintf("%02d", month), "-01", "_",
        year, "-", sprintf("%02d", month), "-", sprintf("%02d", day)
    ))

aln <- aln[which(rownames(aln) %in% metadata$Isolate), ]

ordered_new_tips <- unlist(
    sapply(
        rownames(aln),
        function(x) {
            indx <- grep(pattern = x, metadata$Isolate)
            return(metadata$new_name[indx])
        }
    )
)

rownames(aln) <- ordered_new_tips

# write new aln
# Consult Danielle on whether this needs aligning since its a SNP alignment
write.dna(
    aln,
    format = "fasta",
    file = "empirical_data/shigella.fasta"
)