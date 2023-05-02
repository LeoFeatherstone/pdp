# Adding year, month and date resolution to simulated tree tips
# Also add relaxed clock to TB trees
library(ape)
library(tidyverse)

tree_files <- dir(path = "./sim_study", pattern = ".+.newick.tree")
trees <- lapply(
    tree_files,
    function(x) {
        read.tree(paste0("./sim_study/", x))
    }
)

rename_tips <- function(t) {
    heights <- diag(vcv.phylo(t))

    date_data <- data.frame(
        day = as.Date(date_decimal(heights + 2000))
    )

    date_data <- date_data %>%
        mutate(month = format(day, "%Y-%m-01")) %>%
        mutate(year = format(day, "%Y-01-01")) %>%
        mutate(
            new_name = paste0(
                "t", row_number(), "_", year, "_", month, "_", day
            )
        )

    t$tip.label <- date_data$new_name

    return(t)
}

trees <- lapply(trees, function(x) rename_tips(x))
names(trees) <- tree_files

for (i in seq_along(trees)) {
    write.tree(
        trees[[i]],
        file = paste0("./sim_study/", tree_files[i])
    )
}
