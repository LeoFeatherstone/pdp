#########################################
####### Find Sampling Spans #############
#########################################

library(ape)
library(tidyverse)

files <- list.files(
    path = "analysis/empirical_data",
    pattern = "(h1n1|sars-cov-2|tb|saureus)[.]fasta",
    full.names = TRUE
)

aln <- lapply(
    files,
    function(x) read.dna(x, format = "fasta")
)

names(aln) <- files

dates <- lapply(
    aln,
    function(x) {
        as.Date(gsub(pattern = ".+_", replacement = "", rownames(x)))
    }
)

date_range <- lapply(
    dates,
    function(x) {
        range(dates)
    }
)

## Plot sampling span

sampling_hist <- lapply(dates, function(x) as_tibble(x)) %>%
    bind_rows(.id = "file") %>%
    mutate(
        dataset = gsub(
            file,
            pattern = "analysis/empirical_data/|[.]fasta",
            replacement = ""
        )
    ) %>%
    mutate(dataset = case_when(
        dataset == "tb" ~ "italic(M.~tuberculosis)",
        dataset == "saureus" ~ "italic(S.~aureus)",
        dataset == "h1n1" ~ "H1N1",
        dataset == "sars-cov-2" ~ "SARS-CoV-2"
    )) %>%
    ggplot(aes(x = value)) +
    geom_histogram(alpha = 0.5, col = "black", fill = "dodgerblue") +
    xlab("Date") + ylab("Number of Samples") +
    facet_wrap(~dataset, scales = "free", labeller = label_parsed) +
    theme(
        legend.position = "none",
        text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(
    plot = sampling_hist,
    "figures/empirical_sampling_times.pdf",
    dpi = 300
)
