#!/usr/bin/env Rscript

# Load libraries
library(data.table)
library(rtracklayer)
library(argparse)
library(stringr)
library(tidyverse)
library(dtplyr)
library(arrow)

# Declare constants
SOFTWARE_VERSION <- 1.2

# None for now

# Declare function definitions

# Function to load the mapping information and summary statistics info
load_mapping_files <- function(gwas_file_path, format_file_path) {
  # Load GWAS file names and formats
  gwas_dt <- fread(gwas_file_path)

  # Load the column mapping for summary statistics formats
  format_mapping_dt <- fread(format_file_path)

  # Create a list of format mappings for quick lookups
  format_mappings <- lapply(seq_len(nrow(format_mapping_dt)), function(i) {
    as.list(format_mapping_dt[i, ])
  })
  names(format_mappings) <- format_mapping_dt$format_name

  return(list(gwas_dt = gwas_dt, format_mappings = format_mappings))
}


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
harmonize_gwas_file <- function(gwas_path, format_info, liftover_chain = NULL, output_z_scores = FALSE) {
  build <- format_info$build

  column_types <- list(numeric = c(format_info$p),
                       integer = c(format_info$bp),
                       character = c(format_info$chr, format_info$snp, format_info$a1, format_info$a2))

  # Load the summary statistics file using data.table
  gwas_data <- fread(gwas_path, skip=format_info$skip_lines)
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


check_if_skip <- function(log_file, format_info, gwas_file) {
  if (file.exists(log_file)) {
    log <- fread(log_file)
    print(log)
    comparison_results <- c(log$software_version == SOFTWARE_VERSION, log$original_file == gwas_file,
                 log$format_info_skip_lines == format_info$skip_lines,
                 log$format_info_chromosome == format_info$chromosome,
                 log$format_info_bp == format_info$bp,
                 log$format_info_a1 == format_info$a1,
                 log$format_info_a2 == format_info$a2,
                 log$format_info_snp == format_info$snp,
                 log$format_info_p == format_info$p,
                 log$format_info_z_score == format_info$z_score,
                 log$format_info_effect == format_info$effect,
                 log$format_info_se == format_info$se)
    if (
      all(c("software_version", "format_info_chromosome", "format_info_bp", "format_info_a1",
      "format_info_a2","format_info_snp","format_info_p","format_info_effect","format_info_se") %in% names(log))
      & all(comparison_results, na.rm=T)) {
      return(1)
    }
    return(-1)
  }
  return(0)
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
  parser$add_argument("--gwas-file", type = "character", required = TRUE,
                      help = "Path to the tabulated file containing GWAS file names and formats (e.g., TSV file).")
  parser$add_argument("--format-file", type = "character", required = TRUE,
                      help = "Path to the tabulated file containing summary statistics format mappings (e.g., TSV file).")
  parser$add_argument("--z", required = FALSE, action="store_true", default=FALSE,
                      help = "Instead of beta's and standard errors, output z-scores")
  parser$add_argument("--liftover-chain", type = "character", required = TRUE,
                      help = "Path to the chain file for liftover (e.g., hg19ToHg38.over.chain.gz).")
  parser$add_argument("--eQTLGen-index", type = "character", required = TRUE, default="/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/processed_data/variants/1000G-30x_index.parquet",
                      help = "eQTLGen parquet index file.")

  # Parse arguments
  args <- parser$parse_args(argv)

  eqtlgen_index_file <- "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/processed_data/variants/1000G-30x_index.parquet"
  gwas_tab_file <- "trait_list.txt"
  gwas_format_file <- "format_list.txt"
  liftover_chain_file <- "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/public_data/liftOver/hg19ToHg38.over.chain"
  out_dir <- "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/public_data/linking_disease_relevant_genes/harmonized"
  output_z_scores <- TRUE

  eqtlgen_index_file <- args$eQTLGen_index
  gwas_tab_file <- args$gwas_file
  gwas_format_file <- args$format_file
  liftover_chain_file <- args$liftover_chain
  out_dir <- args$out_dir
  output_z_scores <- args$z

  # Load eQTLGen index
  eqtlgen_index <- arrow::read_parquet(eqtlgen_index_file)
  eqtlgen_index[, mapping_id := paste0(chromosome, "_", bp, "_", pmin(eff_allele, non_eff_allele), "_", pmax(eff_allele, non_eff_allele))]
  eqtlgen_index_extended <- eqtlgen_index %>% group_by(mapping_id) %>%
    filter(n() == 1) %>% ungroup() %>% as.data.table()
  setkey(eqtlgen_index_extended, mapping_id)

  # Load mapping information
  mappings <- load_mapping_files(gwas_tab_file, gwas_format_file)
  gwas_dt <- mappings$gwas_dt
  format_mappings <- mappings$format_mappings

  # Process each GWAS file
  for (i in seq_len(nrow(gwas_dt))) {
    gwas_file <- gwas_dt[i, gwas_file]
    gwas_dir <- gwas_dt[i, gwas_dir]
    format_name <- gwas_dt[i, format_name]
    trait_id <- gwas_dt[i, trait_identifier]

    message(sprintf("Starting '%s' (%i / %i)\nfile: '%s'\nformat: '%s'", trait_id, i, nrow(gwas_dt), gwas_file, format_name))

    # Get format information for the current GWAS file
    format_info <- format_mappings[[format_name]]

    log_file <- file.path(out_dir, "logs", sprintf("%s.log.txt", trait_id))
    skip <- check_if_skip(log_file, format_info, gwas_file)
    if (skip == 1) {
      message("File already processed. Skipping...")
      next
    } else if (skip == -1) {
      message(sprintf("File already processed, but variables were unexpected. Overwriting"))
    }

    gwas_path <- file.path(gwas_dir, gwas_file)

    if (file.exists(gwas_path)) {

      # Harmonize the file
      harmonized_data <- harmonize_gwas_file(gwas_path, format_info, liftover_chain = liftover_chain_file, output_z_scores = output_z_scores)
      harmonized_gwas_data <- harmonized_data$gwas_dt
      original_length <- harmonized_data$original_length
      removed_na_length <- harmonized_data$removed_na_length
      lifted_length <- harmonized_data$lifted_length

      # Save the harmonized data to a new file
      output_file <- file.path(out_dir, sprintf("%s.sumstats.harmonized.gz", trait_id))
      if (output_z_scores) {
        fwrite(harmonized_gwas_data[, .(SNP, CHR, BP, A1, A2, Z, p)], output_file, sep = "\t", row.names=F, col.names=T, quote=F)
      } else {
        fwrite(harmonized_gwas_data[, .(SNP, CHR, BP, A1, A2, EFFECT, p, se)], output_file, sep = "\t", row.names=F, col.names=T, quote=F)
      }
      message(sprintf("Harmonized %s -> %s", gwas_file, output_file))
      message(sprintf("Original rows: %i\nAfter removing rows where some of the remapped columns have NA's: %i\nAfter liftOver: %i\n", original_length, removed_na_length, lifted_length))
      message("Starting merging...")

      harmonized_gwas_data[, mapping_id := paste0(CHR, "_", BP, "_", pmin(A1, A2), "_", pmax(A1, A2))]
      if (output_z_scores) {
        harmonized_data_subset <- harmonized_gwas_data[, .(mapping_id, A1, A2, Z, p)]
      } else {
        harmonized_data_subset <- harmonized_gwas_data[, .(mapping_id, A1, A2, EFFECT, p, se)]
      }

      setkey(harmonized_data_subset, mapping_id)
      merged_gwas_data <- merge(harmonized_data_subset,eqtlgen_index_extended, all=F)
      merged_length <- nrow(merged_gwas_data)
      message(sprintf("After merging: %i", merged_length))

      if (output_z_scores) {
        merged_gwas_data[, z_score := merged_gwas_data[["Z"]]]
        merged_gwas_data[eff_allele != A1,"z_score"] <- merged_gwas_data[eff_allele != A1,z_score*-1]

        swapped_rows <- nrow(merged_gwas_data[z_score != Z])
        message(sprintf("Swapped beta sign: %i / %i (%.3f%%)", swapped_rows, merged_length, swapped_rows/merged_length*100))

        merged_gwas_data[, c("mapping_id", "Z") := NULL, ]
      } else {
        merged_gwas_data[, beta := merged_gwas_data[["EFFECT"]]]
        merged_gwas_data[eff_allele != A1,"beta"] <- merged_gwas_data[eff_allele != A1,beta*-1]

        swapped_rows <- nrow(merged_gwas_data[beta != EFFECT])
        message(sprintf("Swapped beta sign: %i / %i (%.3f%%)", swapped_rows, merged_length, swapped_rows/merged_length*100))

        merged_gwas_data[, c("mapping_id", "EFFECT") := NULL, ]
      }
      merged_gwas_data[, gwas_id := trait_id]

      names(merged_gwas_data)[names(merged_gwas_data) == 'eff_allele'] <- 'str_allele2'
      names(merged_gwas_data)[names(merged_gwas_data) == 'non_eff_allele'] <- 'str_allele1'

      write_dataset(merged_gwas_data, file.path(out_dir, "eqtlgen_parquet_dataset"),
                    format="parquet", partitioning="gwas_id", existing_data_behavior="delete_matching")

      log <- tibble(trait_id = trait_id,
             original_file = gwas_file,
             format = format_name,
             original_length = original_length,
             removed_na_length = removed_na_length,
             lifted_length = lifted_length,
             merged_length = merged_length,
             swapped_rows = swapped_rows,
             format_info = list(format_info),
             software_version = SOFTWARE_VERSION
      ) %>% unnest_wider(format_info, names_sep="_")

      fwrite(log, log_file)

    } else {
      message(sprintf("File %s not found.", gwas_file))
    }
    message(sprintf("Done! %i of %i", i, nrow(gwas_dt)))
  }

}

# Entry point check
if (sys.nframe() == 0 && !interactive()) {
  main()
}
