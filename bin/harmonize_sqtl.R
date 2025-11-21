#!/usr/bin/env Rscript

# Load libraries
library(data.table)
library(argparse)
library(stringr)
library(tidyverse)
library(dtplyr)
library(arrow)

# Declare constants
SOFTWARE_VERSION <- 1.2

INTERVAL_SQTL_FORMAT_INFO <- list(
  "snp" = "variant_id",
  "chr" = "chr",
  "bp" = "pos_b38",
  "a1" = "effect_allele",
  "a2" = "other_allele",
  "effect" = "slope",
  "se" = "slope_se",
  "p" = "pval_nominal",
  "build" = "hg38")


# None for now

# Declare function definitions


# Define a function to check if required columns exist and are non-NA
check_columns <- function(data, required_columns) {
  missing_columns <- setdiff(required_columns, names(data))  # Find missing columns
  names(missing_columns) <- names(required_columns)[required_columns %in% missing_columns] # Add names to missing columns

  # Check for missing columns
  if (length(missing_columns) > 0) {
    message(sprintf("The following required columns are missing: ,", paste(names(missing_columns), collapse = ", ")))
  } else {
    message("All required columns are present.")
  }

  columns_with_na <- required_columns[sapply(data[, ..required_columns], function(col) any(is.na(col)))]  # Columns containing NA
  columns_with_all_na <- required_columns[sapply(data[, ..required_columns], function(col) all(is.na(col)))]  # Columns containing NA

  # Check for columns with NA values
  if (length(columns_with_na) > 0) {
    message("The following columns contain some NA values: ", paste(columns_with_na, collapse = ", "))
  } else {
    message("All required columns have non-NA values.")
  }

  # Check for columns with NA values
  if (length(columns_with_all_na) > 0) {
    message("The following columns contain just NA values: ", paste(columns_with_all_na, collapse = ", "))
  }

  # If everything is good
  if (length(missing_columns) == 0 && length(columns_with_all_na) == 0) {
    message("All checks passed. Data is ready for further processing.")
    return(TRUE)
  } else {
    warning("Some issues detected with the data. Please review the messages above.")
    print(head(data))
    return(FALSE)
  }
}


# Harmonize each GWAS file
harmonize_gwas_file <- function(gwas_data, format_info, liftover_chain = NULL, output_z_scores = FALSE, output_sample_size = TRUE) {
  build <- format_info$build

  column_types <- list(numeric = c(format_info$p),
                       integer = c(format_info$bp),
                       character = c(format_info$chr, format_info$snp, format_info$a1, format_info$a2))

  # Load the summary statistics file using data.table
  original_length <- nrow(gwas_data)

  gwas_data[, (column_types[["numeric"]]) := lapply(.SD, as.numeric), .SDcols = column_types[["numeric"]]]
  gwas_data[, (column_types[["integer"]]) := lapply(.SD, as.numeric), .SDcols = column_types[["integer"]]]

  required_columns <- c("chr" = format_info$chr, "bp" = format_info$bp, "snp" = format_info$snp,
                        "a1" = format_info$a1, "a2" = format_info$a2, "p" = format_info$p)

  if (format_info$snp == "" | !(format_info$snp %in% names(gwas_data))) {
    gwas_data[, "snp" := paste(get(required_columns['chr']),
                               get(required_columns['bp']),
                               get(required_columns['a1']),
                               get(required_columns['a2']), sep="_")]
    required_columns["snp"] <- "snp"
  }

  gwas_data[, required_columns['chr'] := gsub("^0+", "", get(required_columns['chr']))]

  check_columns(gwas_data, required_columns)

  # Rename columns to a standardized format
  setnames(gwas_data, old = required_columns,
           new = c("CHR", "BP", "SNP", "A1", "A2", "p"))

  # Check and convert p-values
  if (max(gwas_data$p) > 1 & str_detect(tolower(format_info$p), "log")) {
    message("Detected -log10 p-values. Converting...")
    gwas_data[, p := 10^(-gwas_data$p)]
  } else if (max(gwas_data$p) > 1) {
    warning("Max p-value is greater than 1, but column name does not indicate log transformation. Check this!")
  } else if (str_detect(tolower(format_info$p), "log")) {
    warning("Max p-value is not greater than 1, but column name indicates log transformation. Check this!")
  }

  # Check and convert effect estimates
  gwas_data <- convert_effect_to_beta_dynamically(gwas_data, format_info)
  if ("EFFECT" %in% colnames(gwas_data)) {
    gwas_data <- map_standard_error_dynamically(gwas_data, format_info)
  }

  if (output_z_scores) {
    gwas_data <- map_z_score(gwas_data, format_info)
    columns_to_check <- c("CHR", "BP", "SNP", "A1", "A2", "p", "Z")
  } else {
    columns_to_check <- c("CHR", "BP", "SNP", "A1", "A2", "p", "se", "EFFECT")
  }

  if (output_sample_size) {
    gwas_data[, n := as.integer(gwas_data[[format_info$n]])]
    columns_to_check <- c(columns_to_check, "N")
  }

  gwas_data <- gwas_data[rowSums(is.na(gwas_data[, columns_to_check, with = FALSE])) == 0, ]
  removed_na_length <- nrow(gwas_data)

  # Harmonize to build 38 if necessary
  if (build != "hg38" && !is.null(liftover_chain)) {
    gwas_data <- lift_to_build_38(gwas_data, liftover_chain)
  }

  lifted_length <- nrow(gwas_data)

  return(list(gwas_dt = gwas_data,
              original_length = original_length,
              removed_na_length = removed_na_length,
              lifted_length = lifted_length))
}


# Function to map gwas data to z-scores
map_z_score <- function(gwas_data, format_info) {
  z_score_column <- format_info$z_score

  # Case 1: Use Z score column
  if (!is.null(z_score_column) && z_score_column != "" && z_score_column %in% names(gwas_data)) {
    message("Using provided Z-score column.")
    gwas_data[, Z := gwas_data[[z_score_column]]]
    return(gwas_data)
  } else if (all(c("EFFECT", "p") %in% names(gwas_data))) {
    gwas_data[, Z := abs(qnorm(gwas_data[["p"]] / 2, lower.tail=F)) * sign(gwas_data[["EFFECT"]])]
    return(gwas_data)
  } else {
    warning("No Z-score column available, and could not compute Z-scores using effect-size column and P-values.")
  }
  return(gwas_data)
}


# Function to dynamically map standard errors from SE or CI columns
map_standard_error_dynamically <- function(gwas_data, format_info) {
  # Check for the presence of standard error and confidence interval columns
  se_column <- format_info$se
  ci_lower_column <- format_info$ci_lower
  ci_upper_column <- format_info$ci_upper

  # Case 1: Use SE column if present and non-NA
  if (!is.null(se_column) && se_column %in% names(gwas_data) && !all(is.na(gwas_data[[se_column]]))) {
    message("Using provided SE column.")
    gwas_data[, se := gwas_data[[se_column]]]
    return(gwas_data)
  } else if (!is.null(ci_lower_column) && ci_lower_column %in% names(gwas_data) &&
    !is.null(ci_upper_column) && ci_upper_column %in% names(gwas_data)) {
    # Case 2: Use CI columns if SE is not available

    # Check if the CI columns are non-NA
    if (!all(is.na(gwas_data[[ci_lower_column]])) && !all(is.na(gwas_data[[ci_upper_column]]))) {

      message("SE column not found. Calculating SE from confidence intervals.")
      # Calculate the midpoint of CI to determine if we need to log-transform the values
      mean_of_ci_lower_and_upper <- mean(c(gwas_data[[ci_lower_column]], gwas_data[[ci_upper_column]]), na.rm = TRUE)

      likely_beta_like <- abs(mean_of_ci_lower_and_upper - 1) > abs(mean_of_ci_lower_and_upper)

      # Check if CIs are closer to 1 (suggesting ORs), and log-transform if necessary
      if (!likely_beta_like) {
        message("Confidence intervals suggest Odds Ratios (ORs). Applying log transformation.")
        ci_lower_log <- log(gwas_data[[ci_lower_column]])
        ci_upper_log <- log(gwas_data[[ci_upper_column]])
        # Calculate SE from log-transformed CIs
        gwas_data[, se := (ci_upper_log - ci_lower_log) / (2 * qnorm(0.975, lower.tail=F))]
      } else {
        message("Confidence intervals suggest Beta-like values. No log transformation.")
        # Calculate SE directly from CIs
        gwas_data[, se := (gwas_data[[ci_upper_column]] - gwas_data[[ci_lower_column]]) / (2 * 1.96)]
      }
      return(gwas_data)
    }
  }
  message("No valid SE or CI columns found in the data.")
  message("Using p-value to estimate stanard errors.")

  # Step 1: Compute the z-score from the two-sided p-value
  z_score <- abs(qnorm(gwas_data[["p"]] / 2, lower.tail=F))

  # Step 2: Compute the standard error (SE)
  gwas_data[, se := (gwas_data[["EFFECT"]] / z_score)]

  # If none of the above cases are valid, warn and return data without SE
  return(gwas_data)
}


# Dynamically determine if the effect is in Beta or OR and convert to Beta if necessary
convert_effect_to_beta_dynamically <- function(gwas_data, format_info) {
  # Split the effect types on ";" to get available options
  effects <- unlist(strsplit(format_info$effect, ";"))  # Keep original case

  # Initialize a variable to hold the selected effect column
  selected_effect_column <- NULL

  # Check which effect types are present in the summary statistics columns
  for (effect_column in effects) {
    # Check if the respective column exists in the GWAS data
    if (effect_column %in% names(gwas_data)) {
      if (!all(is.na(gwas_data[[effect_column]]))) {
        selected_effect_column <- effect_column  # Found a matching column
        break  # Exit loop after the first match
      }
    }
  }

  # If no valid effect column was found, return the original data
  if (is.null(selected_effect_column)) {
    message("No valid effect column found in summary statistics.")
    return(gwas_data)
  }

  print("Selected effect-size column:")
  print(selected_effect_column)

  # Assign the selected effect values to the EFFECT column
  gwas_data[, EFFECT := gwas_data[[selected_effect_column]]]

  # Calculate the average of the EFFECT column
  mean_effect <- mean(gwas_data$EFFECT, na.rm = TRUE)
  print("Mean effect:")
  print(mean_effect)

  # Is closer to 0 or 1
  likely_beta_like <- abs(mean_effect - 1) > abs(mean_effect)
  print("Likely beta-like:")
  print(likely_beta_like)

  # Determine the likely type of the effect based on the column name
  is_log_transformed <- grepl("log|logodds|logor", selected_effect_column, ignore.case = TRUE)
  is_odds_ratio <- grepl("odds_ratio|or|odds|ratio", selected_effect_column, ignore.case = TRUE)
  is_beta <- grepl("beta|b", selected_effect_column, ignore.case = TRUE)

  # Evaluate the effect type based on average and column name
  if (is_log_transformed) {
    if (likely_beta_like) {
      message("Detected Log-transformed Odds Ratio (OR) values. No conversion needed.")
    } else {
      warning("Log-transformed OR values appear to be very large, possibly odds ratios.")
    }
  } else if (is_odds_ratio) {
    if (!likely_beta_like) {
      message("Detected Odds Ratio (OR) values. Converting OR to Beta.")
      gwas_data[, EFFECT := log(EFFECT)]  # Convert OR to Beta using log transformation
    } else {
      message("Columns indicate odds Ratio (OR) values are used, but average indicates differently")
    }
  } else if (is_beta) {
    if (likely_beta_like) {
      message("Detected Beta values. No conversion needed")
    } else {
      message("Columns indicate Beta values, but average indicates differently")
    }
  } else {
    message("Unknown effect type in the format mapping based on column name and average.")
  }

  return(gwas_data)
}


# Function to lift genome build to hg38 using rtracklayer (LiftOver)
lift_to_build_38 <- function(gwas_data, liftover_chain_file) {
  # Create a GRanges object for liftover
  gwas_data[, index := seq_along(gwas_data$SNP)]

  gr <- GRanges(
    seqnames = Rle(paste0("chr", str_remove(gwas_data$CHR, "^chr"))),
    ranges = IRanges(start = gwas_data$BP, end = gwas_data$BP),
    strand = Rle("*"),
    index = gwas_data$index
  )

  # Perform liftover
  chain <- import.chain(liftover_chain_file)
  lifted_gr <- liftOver(gr, chain)

  # Flatten the result back to the data.table
  new_coords <- unlist(lifted_gr)
  new_coords_dt <- as.data.table(new_coords) %>% group_by(index) %>%
    filter(n() == 1) %>% as.data.table()

  gwas_data_filtered <- gwas_data[new_coords_dt[["index"]], , ]
  gwas_data_filtered[, CHR := sub("chr", "", new_coords_dt[["seqnames"]])]
  gwas_data_filtered[, BP := new_coords_dt[["start"]]]
  gwas_data_filtered[ , index := NULL, ]

  return(gwas_data_filtered)
}

# Read in pQTLs
readpqtl <- function(fn, allele_frequency){
  message(sprintf("Reading file: %s...", fn))
  dt_temp <- fread(fn)
  dt_temp <- dt_temp[A1FREQ > allele_frequency & A1FREQ < 1 - allele_frequency]
  keycols <- c("ID")
  setkeyv(dt_temp, keycols)
  return(dt_temp)
}


# Main function to run the pipeline
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = TRUE)
  }

  # Set up argument parser
  parser <- ArgumentParser(description = "Harmonize GWAS summary statistics files.")

  # Add arguments
  parser$add_argument("--out-dir", type = "character", required = TRUE,
                      help = "Base directory where the processed summary statistics files will be stored.")
  parser$add_argument('--sqtl-folder', type = 'character',
                      help = 'A tab-delimited file with pQTL results.')
  parser$add_argument("--eQTLGen-index", type = "character", required = TRUE, default="/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/processed_data/variants/1000G-30x_index.parquet",
                      help = "eQTLGen parquet index file.")

  # Parse arguments
  args <- parser$parse_args(argv)

  message(sprintf("Starting '%s'\nfile: '%s'", args$pqtl_id, args$sqtl_folder))

  sqtls <- fread(args$sqtl_folder)

  message("pQTLs read in!")

  eqtlgen_index_extended <- arrow::read_parquet(args$eQTLGen_index)
  setkey(eqtlgen_index_extended, mapping_id)

  # Harmonize the file
  harmonized_data <- harmonize_gwas_file(sqtls, INTERVAL_SQTL_FORMAT_INFO, output_z_scores = FALSE, output_sample_size = FALSE)
  harmonized_gwas_data <- harmonized_data$gwas_dt
  original_length <- harmonized_data$original_length
  removed_na_length <- harmonized_data$removed_na_length
  lifted_length <- harmonized_data$lifted_length

  # Save the harmonized data to a new file
  output_file <- file.path(args$out_dir, sprintf("sumstats.harmonized.gz"))
  fwrite(harmonized_gwas_data[, .(SNP, CHR, BP, A1, A2, EFFECT, p, se, phenotype_id)], output_file, sep = "\t", row.names=F, col.names=T, quote=F)
  message(sprintf("Harmonized %s -> %s", args$pqtl_folder, output_file))
  message(sprintf("Original rows: %i\nAfter removing rows where some of the remapped columns have NA's: %i\nAfter liftOver: %i\n", original_length, removed_na_length, lifted_length))
  message("Starting merging...")

  harmonized_gwas_data[, mapping_id := paste0(CHR, "_", BP, "_", pmin(A1, A2), "_", pmax(A1, A2))]
  harmonized_data_subset <- harmonized_gwas_data[, .(mapping_id, A1, A2, EFFECT, p, se, phenotype_id, CHR)]

  setkey(harmonized_data_subset, mapping_id)
  merged_gwas_data <- merge(harmonized_data_subset,eqtlgen_index_extended, all=F)
  merged_length <- nrow(merged_gwas_data)
  message(sprintf("After merging: %i", merged_length))

  merged_gwas_data[, beta := merged_gwas_data[["EFFECT"]]]
  merged_gwas_data[eff_allele != A1,"beta"] <- merged_gwas_data[eff_allele != A1,beta*-1]

  swapped_rows <- nrow(merged_gwas_data[beta != EFFECT])
  message(sprintf("Swapped beta sign: %i / %i (%.3f%%)", swapped_rows, merged_length, swapped_rows/merged_length*100))

  merged_gwas_data[, c("mapping_id", "EFFECT") := NULL, ]

  names(merged_gwas_data)[names(merged_gwas_data) == 'eff_allele'] <- 'str_allele2'
  names(merged_gwas_data)[names(merged_gwas_data) == 'non_eff_allele'] <- 'str_allele1'

  write_dataset(arrange(merged_gwas_data, variant_index), file.path(args$out_dir, "eqtlgen_parquet_dataset"),
                format="parquet", partitioning="CHR", existing_data_behavior="delete_matching")

  log <- tibble(trait_id = merged_gwas_data$CHR[1],
                original_file = args$sqtl_folder,
                format = "INTERVAL_SQTL",
                original_length = original_length,
                removed_na_length = removed_na_length,
                lifted_length = lifted_length,
                merged_length = merged_length,
                swapped_rows = swapped_rows,
                format_info = list(INTERVAL_SQTL_FORMAT_INFO),
                software_version = SOFTWARE_VERSION
  ) %>% unnest_wider(format_info, names_sep="_")

  log_file <- file.path(args$out_dir, sprintf("%s.log.txt", merged_gwas_data$CHR[1]))

  fwrite(log, log_file)
}

# Entry point check
if (sys.nframe() == 0 && !interactive()) {
  main()
}
