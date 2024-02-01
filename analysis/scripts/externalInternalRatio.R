########################################################
######## External:Internal Branch Length ###############
########################################################

library(ape)
library(tidyverse)

ei_ratio <- function(tree) {
    index_external <- which(tree$edge[, 2] %in% seq_along(tree$tip.label))

    sum_external <- sum(tree$edge.length[index_external])
    sum_internal <- sum(sum(tree$edge.length[-index_external]))

    external_internal_ratio <- sum_external / sum_internal

    return(external_internal_ratio)
}


# saureus case
day_trees <- read.nexus("analysis/empirical_data/saureus_BD_SC_Day-Fixed.trees")
month_trees <- read.nexus("analysis/empirical_data/saureus_BD_SC_Month-Fixed.trees")
year_trees <- read.nexus("analysis/empirical_data/saureus_BD_SC_Year-Fixed.trees")

# burnin
day_trees <- day_trees[-c(1:100)]
month_trees <- month_trees[-c(1:100)]
year_trees <- year_trees[-c(1:100)]

# ei ratio
ei_day <- unlist(lapply(day_trees, function(x) ei_ratio(x)))
ei_month <- unlist(lapply(month_trees, function(x) ei_ratio(x)))
ei_year <- unlist(lapply(year_trees, function(x) ei_ratio(x)))

df <- data.frame("month" = ei_month, "day" = ei_day, "year" = ei_year)

df %>%
    pivot_longer(everything(), names_to = "resolution", values_to = "eiRatio") %>%
    ggplot() +
    geom_histogram(aes(x = eiRatio, fill = resolution), position = "identity") +
    scale_fill_manual(values = alpha(c("red", "blue", "black"), 0.3))
