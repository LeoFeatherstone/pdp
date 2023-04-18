#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(TreeSim))
suppressPackageStartupMessages(require(lubridate))

args <- commandArgs(trailingOnly = TRUE)
set.seed(1234)

main <- function() {

    r_naught <- as.numeric(args[1])
    delta <- as.numeric(args[2])
    p_samp <- as.numeric(args[3])
    time <- as.numeric(args[4])
    num_trees <- as.numeric(args[5])
    prefix <- args[6]

    print(time)

    trees <- list()
    class(trees) <- "multiPhylo"
    counter <- 1
    num_passed_sims <- 0
    min_tips <- 10

    while (
        (counter < 10 * num_trees)
        &&
        (length(trees) < num_trees)
        ) {
            sim_tree <- sim.bdsky.stt(
                n = 0,
                timesky = 0,
                lambdasky = r_naught * delta,
                deathsky = delta,
                sampprobsky = p_samp,
                timestop = time
            )[[1]]

        if ((class(sim_tree) == "phylo")) {
            if (length(sim_tree$tip.label) >= min_tips) {
                num_passed_sims <- num_passed_sims + 1
                trees[[num_passed_sims]] <- sim_tree
            }
        }
        counter <- counter + 1
    }

    for (i in seq_along(trees)) {

        tip_dates <- as.Date(date_decimal(2000 + diag(vcv.phylo(trees[[i]]))))
        trees[[i]]$tip.label <- paste0(
            trees[[i]]$tip.label,
            "_",
            tip_dates,
            "_",
            gsub(".{3}$", "", tip_dates)
        )
    }

    write.tree(trees, file = paste0(prefix, ".trees"))
}
main()