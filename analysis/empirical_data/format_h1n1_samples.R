## One-off formatting TB samples from source
library(ape)
library(tidyverse)
library(lubridate)

aln <- read.dna(
    "empirical_data/NorthAm.Nov.fasta",
    format = "fasta"
)

dates <- data.frame(
    dec_date = as.numeric(
    gsub(
        pattern = "^.*\\|", replacement = "", rownames(aln)
    )
    ),
    whole_date = as.Date(
    date_decimal(
    as.numeric(
    gsub(
        pattern = "^.*\\|", replacement = "", rownames(aln)
    )
    ))),
    tip_name = gsub(
        pattern = "\\|([^|]*)$", replacement = "", rownames(aln)
    )
)


dates <- dates %>%
    mutate(year = format(whole_date, "%Y-06-15")) %>%
    mutate(month = format(whole_date, "%Y-%m-15")) %>%
    mutate(day = format(whole_date, "%Y-%m-%d")) %>%
    mutate(new_name = paste0(
       row_number(), "_", year, "_", month, "_", day
    ))

# min sample date is 2009-04-01
# Include up to end of june: decimal_date(as.Date("2009-06-30"))
rownames(aln) <- dates$new_name

exp_tips <- dates %>%
    filter(dec_date <= decimal_date(as.Date("2009-06-30"))) %>%
    select(new_name)


# write new aln
write.dna(
    aln[which(rownames(aln) %in% exp_tips$new_name),],
    format = "fasta",
    file = "empirical_data/h1n1.fasta"
)
