## One-off formatting TB samples from source
library(ape)
library(tidyverse)
library(lubridate)

aln <- read.dna(
    "empirical_data/st_usa300_dates_no_outgroups.fasta",
    format = "fasta"
)

dates <- data.frame(
    whole_date = as.Date(
        date_decimal(
            as.numeric(
                gsub(pattern = ".+_", replacement = "", rownames(aln))
            ))),
    tip_name = gsub(pattern = "_([^_]*$)", replacement = "", rownames(aln))
)


dates <- dates %>%
    mutate(year = format(whole_date, "%Y-06-15")) %>%
    mutate(month = format(whole_date, "%Y-%m-15")) %>%
    mutate(day = format(whole_date, "%Y-%m-%d")) %>%
    mutate(new_name = paste0(
       tip_name, "_", year, "_", month, "_", day
    ))

rownames(aln) <- dates$new_name

# write new aln
write.dna(
    aln,
    format = "fasta",
    file = "empirical_data/saureus.fasta"
)
