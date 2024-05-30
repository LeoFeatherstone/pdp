########################################################
############## Calculate tree imbalance ################
########################################################
# Check if treebalance & parallel loaded and/or load
if (!requireNamespace("treebalance", quietly = TRUE)) {
  install.packages("treebalance", repos = "https://mirror.aarnet.edu.au/pub/CRAN/")
}
if (!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel", repos = "https://mirror.aarnet.edu.au/pub/CRAN/")
}

library(treebalance)
library(parallel)
library(ape)
library(tidyverse)

# Function to calculate normalized coless imbalance over tree posts 
get_tree_stats <- function(out_path) {

  # Step 0: Get base Tree data
  base_tree_files <- list.files(path = "posterior_simulations/", pattern = "\\.tree$", full.names = TRUE)
  base_trees <- lapply(
    base_tree_files,
    function(x) read.tree(x)
  )
  
  base_tree_stats <- lapply(
    base_trees,
    function(x) {
        data.frame(
            n_tips = length(x$tip.label),
            sampling_span = max(diag(vcv.phylo(x))) - min(diag(vcv.phylo(x)))
        )
    }
  )

  names(base_tree_stats) <- gsub(
    pattern = ".+\\/", 
    replacement = "", 
    x = gsub(base_tree_files, pattern = "[.]tree", replacement = "")
  )

  base_tree_stats_df <- bind_rows(base_tree_stats, .id = "id") %>%
    separate(id, into = c("organism", "replicate"), sep = "_")

  print("Base tree data complete!")

  # Step 1: Get log files
  log_files <- list.files(path = "posterior_simulations/", pattern = "\\.trees$", full.names = TRUE)
  
  # Step 2: Initialize an empty list to store results
  imbalance <- list()
  
  # Step 3: Process each log file individually
  for(file in log_files) {

    print(paste(c("Reading", file)))
    t1 <- Sys.time()

    trees <- ape::read.nexus(file)
    # NB, taking second half of posterior to save space & time (trees[-c(1:(0.5*length(trees)))])
    trees <- trees[-c(1:(0.5 * length(trees)))]
    
    # Step 4: Calculate tree imbalance for each tree
    # NB, taking second half of posterior to save space & time (trees[-c(1:(0.5*length(trees)))])
    tree_stats <- parallel::mclapply(trees, function(tr) {
      data.frame(
        coless_imbalance = treebalance::collessI(tr, method = "corrected"),
        sackin_index = treebalance::sackinI(tr),
        ei_ratio = ei_ratio(tr),
        longest_branch = longest_branch(tr),
        Sample = as.numeric(gsub(names(tr), pattern = "STATE_", replacement = ""))
      )
    })
    
    # Step 5: Bind tree stats with sample step and append to imbalance list
     imbalance[[file]] <- (
        tree_stats %>%
        bind_rows(.id = "Sample") %>%
        summarise(
            coless_imbalance = mean(coless_imbalance),
            sackin_index = mean(sackin_index),
            ei_ratio = mean(ei_ratio),
            longest_branch = mean(longest_branch)
        )
    )

    t2 <- Sys.time()

    print(paste(file, "processed!", t2 - t1))
  }

#   print("Up to 5.5")  
#   # Step 5.5: Test saving earlier to see if tidy is the memory issue
#   save(tree_stats, base_tree_stats_df, file = out_path)

 print("Up to 5.5")  
  names(imbalance) <- gsub(
    gsub(log_files, pattern = "[.]trees", replacement = ""),
    pattern = ".+\\/", replacement = ""
  )
  
  print("Up to 6")  
  # Step 6: Combine all results into a single dataframe
  imbalance_df <- bind_rows(imbalance, .id = "id") %>%
    separate(
        id,
        into = c("organism", "treePrior", "clockodel", "resolution", "replicate"),
        sep = "_"
    )
  print("Up to 7")  
  # Step 7: combine data frames
  tree_stats <- imbalance_df %>%
    left_join(base_tree_stats_df)

  print("Up to 8")   
  # Step 8: Save the result
  t1 <- Sys.time()
  saveRDS(tree_stats, file = out_path)
  t2 <- Sys.time()
  print(paste(out_path, "processed!", t2 - t1))
}
 
ei_ratio <- function(tree) {
    index_external <- which(tree$edge[, 2] %in% seq_along(tree$tip.label))

    sum_external <- sum(tree$edge.length[index_external])
    sum_internal <- sum(sum(tree$edge.length[-index_external]))

    external_internal_ratio <- sum_external / sum_internal

    return(external_internal_ratio)
}

longest_branch <- function(tree) {
    
    longest_branch <- max(tree$edge.length)
    return(longest_branch)
}

# Call the function with output path
get_tree_stats(snakemake@output[[1]])