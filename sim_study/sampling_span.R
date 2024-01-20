####################################################
######## Simulation sampling span confound #########
####################################################

library(ape)
library(tidyverse)

t_files <- dir("./sim_study", pattern = ".tree", full.names = TRUE)

trees <- lapply(
    t_files,
    function(x) read.tree(x)
)

names(trees) <- gsub(
    gsub(t_files, pattern = "[.]newick[.]tree", replacement = ""),
    pattern = "./sim_study/", replacement = ""
)

get_sampling_span <- function(tree) {
    height_above_root <- diag(ape::vcv.phylo(tree))
    return(max(height_above_root) - min(height_above_root))
}

span <- unlist(lapply(trees, function(x) get_sampling_span(x)))
df <- data.frame(span = span, treatment = names(span))

df %>%
    separate_wider_delim(treatment, delim = "_", names = c("microbe", "id")) %>%
    ggplot() +
    geom_histogram(aes(x = span * 365.25, fill = microbe)) +
    facet_wrap(~microbe, scales = "free", ncol = 1)
