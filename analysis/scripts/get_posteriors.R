########################################################
########## Combine posteriors to data frame ############
######### Should probably combine with ess.R ###########
########################################################

require(tidyverse)

# Function to calculate and save ESS values
get_posteriors <- function(out_path) {
  # Step 1: Get log files
  log_files <- dir(
    path = "posterior_simulations/",
    pattern = ".+[.]log",
    full.names = TRUE
  )
  
  # Step 2: Read log files
  trace_data <- lapply(
    log_files,
    function(file) { read.table(file, header = TRUE) }
  )

  # Step 3: Add names
  names(trace_data) <- gsub(
    gsub(log_files, pattern = "[.]log", replacement = ""),
        pattern = ".+\\/", replacement = ""
  )

  # Step 4: bind rows
  trace_data <- trace_data %>%
    bind_rows(.id = "id") %>%
    separate_wider_delim(
        id,
        delim = "_",
        names = c("organism", "clockModel", "treePrior", "resolution", "replicate")
    )

  # Step 5: Save ESS values
  saveRDS(trace_data, file = out_path)
}

# Call the function with output path
get_posteriors(snakemake@output[[1]])
