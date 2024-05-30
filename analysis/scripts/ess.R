########################################################
########## Calculate ESS and filter data ###############
########################################################
# Check if tracerer loaded and/or load
if (!requireNamespace("tracerer", quietly = TRUE)) {
  install.packages("tracerer", repos = "https://mirror.aarnet.edu.au/pub/CRAN/")
}
library(tracerer)

# Function to calculate and save ESS values
get_sim_ess <- function(out_path) {
  # Step 1: Get log files
  log_files <- dir(
    path = "posterior_simulations/",
    pattern = ".+[.]log",
    full.names = TRUE
  )
  
  # Step 2: Read non-empty log files
  trace_data <- lapply(
    log_files,
    function(file) {
        read.table(file, header = TRUE)
    }
  )
  
  # Step 3: Test burnin & calculate ESS
  burnin <- c(seq(from = 0.1, to = 0.5, by = 0.05))
  ess_values <- list()

  for (i in seq_along(burnin)) {
    
    df <- lapply(
      trace_data,
      function(data) data[-ceiling(burnin[i] * nrow(data)), ]
    )

    # Calculate ESS
    ess_values[[i]] <- lapply(
      df,
      function(x) calc_esses(x, 50000)
    )

    # Add names
    names(ess_values[[i]]) <- gsub(
      gsub(log_files, pattern = "[.]log", replacement = ""),
      pattern = ".+\\/", replacement = ""
    )

  }

  # Step 6: Save ESS values
  saveRDS(ess_values, file = out_path)
}

# Call the function with output path
get_sim_ess(snakemake@output[[1]])
