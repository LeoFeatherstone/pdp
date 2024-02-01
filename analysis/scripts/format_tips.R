# Adding year, month and date resolution to simulated tree tips
require(ape)
require(lubridate)

round_tip_dates <- function(tree) {
    heights <- diag(vcv.phylo(tree))

    day <- as.Date(date_decimal(heights + 2000))
    month <- format(day, "%Y-%m-15")
    year <- format(day, "%Y-06-15")

    new_labels <- paste0("t", seq_along(day), "_", day, "_", month, "_", year)
    return(new_labels)
}

round_all_tree_dates <- function(data_path, out_path, threads, myparam) {
    tree <- ape::read.nexus(data_path)
    tree$tip.label <- round_tip_dates(tree)

    write.tree(tree, file = out_path)
}

round_all_tree_dates(
    snakemake@input[[1]], snakemake@output[[1]],
    snakemake@threads, snakemake@config[["myparam"]]
)
