## One-off formatting H3N2 samples from source
library(ape)
library(tidyverse)
library(lubridate)

aln <- read.dna(
    "analysis/empirical_data/h3n2_2deme.fa",
    format = "fasta"
)

dates <- data.frame(
    dec_date = as.numeric(
    gsub(
        pattern = "^.*\\_", replacement = "", rownames(aln)
    )
    ),
    whole_date = as.Date(
    date_decimal(
    as.numeric(
    gsub(
        pattern = "^.*\\_", replacement = "", rownames(aln)
    )
    ))),
    tip_name = gsub(
        pattern = "\\|([^_]*)$", replacement = "", rownames(aln)
    )
)


dates <- dates %>%
    mutate(year = format(whole_date, "%Y-06-15")) %>%
    mutate(month = format(whole_date, "%Y-%m-15")) %>%
    mutate(day = format(whole_date, "%Y-%m-%d")) %>%
    mutate(new_name = paste0(
       tip_name, "_", year, "_", month, "_", day
    ))

# min sample date is 2009-04-01
# Include up to end of june: decimal_date(as.Date("2009-06-30"))
rownames(aln) <- dates$new_name


# write new aln
write.dna(
    aln,
    format = "fasta",
    file = "analysis/empirical_data/h3n2.fasta"
)
