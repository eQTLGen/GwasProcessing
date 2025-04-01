#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(arrow)
library(data.table)

# Declare constants

# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Process input
  eqtlgen_index_file <- "1000G-30x_index.parquet"
  eqtlgen_index <- arrow::read_parquet(eqtlgen_index_file)
  eqtlgen_index[, mapping_id := paste0(chromosome, "_", bp, "_", pmin(eff_allele, non_eff_allele), "_", pmax(eff_allele, non_eff_allele))]
  eqtlgen_index_extended <- eqtlgen_index %>% group_by(mapping_id) %>%
    filter(n() == 1) %>% ungroup() %>% as.data.table()

  write_parquet(eqtlgen_index_extended, "1000G-30x_index_with_unique_mapping_ids_50M.parquet")

  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}