#!/usr/bin/env Rscript
# =============================================================================
# convert_epic_samplesheet.R
#
# Author: Jelle R. Dalenberg
# Created: 24-Feb-2026
#
# Description:
#   Convert an Illumina methylation SampleSheet (with [Header]/[Data] sections)
#   into a pipeline-ready CSV with the following columns:
#
#     Sample_ID, Well, Sample_Name, Plate, Slide, Array, Gender, Basename
#
#   The script:
#     - Extracts the [Data] section from an Illumina-style SampleSheet.
#     - Validates the presence of required columns.
#     - Constructs a SHORT 'Basename' in the exact form required by the pipeline:
#           <Slide>/<Slide>_<Array>
#       (e.g., 209548820074/209548820074_R01C01)
#     - Optionally checks that the corresponding *_Grn.idat and *_Red.idat files
#       exist on disk. Checking prepends --basename-prefix to the
#       short Basename, but the output file always contains the short Basename.
#     - Optionally enforces uniqueness of Slide + Array (Sentrix_ID + Position).
#
# Checks performed:
#   Input parsing:
#     --input            → Must exist and contain a [Data] section.
#   Columns in [Data]:
#     Must include: Sample_Well, Sample_Name, Sample_Plate,
#                   Sentrix_ID, Sentrix_Position, Gender
#   Duplicate positions:
#     If --stop-on-duplicate is provided, stop when any duplicate (Slide, Array)
#     is found; otherwise, emit a warning.
#   File existence (optional):
#     If --check-files is provided, verify that for every row the following exist:
#       <prefix>/<Slide>/<Slide>_<Array>_Grn.idat
#       <prefix>/<Slide>/<Slide>_<Array>_Red.idat
#     where <prefix> is --basename-prefix (may be empty).
#
# Usage:
#   Rscript convert_epic_samplesheet.R \
#       --input S2025-71_EPIC_SampleSheet.csv \
#       --output output_samplesheet.csv \
#       [--basename-prefix /abs/or/rel/path/to/idats] \
#       [--check-files] \
#       [--stop-on-duplicate]
#
# Output columns:
#   Sample_ID, Well, Sample_Name, Plate, Slide, Array, Gender, Basename
#
# Notes:
#   - This script does not include 'Sample_Group' or any other custom columns.
#   - Intended for EPIC/EPICv2-style IDAT stems:
#         <Slide>/<Slide>_<Array>_Grn.idat
#         <Slide>/<Slide>_<Array>_Red.idat
# =============================================================================

suppressWarnings(suppressMessages({
  requireNamespace("readr", quietly = TRUE) || stop("Package 'readr' is required.")
  requireNamespace("dplyr", quietly = TRUE) || stop("Package 'dplyr' is required.")
  requireNamespace("stringr", quietly = TRUE) || stop("Package 'stringr' is required.")
}))

# -----------------------------
# Minimal command-line flag parser
# -----------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  get_val <- function(flag, default = NULL, takes_value = TRUE) {
    i <- which(args == flag)
    if (length(i) == 0) {
      if (isTRUE(takes_value)) return(default) else return(FALSE)
    }
    if (isFALSE(takes_value)) return(TRUE)
    if (i == length(args)) stop("Flag ", flag, " requires a value.")
    return(args[i + 1])
  }

  list(
    input             = get_val("--input"),
    output            = get_val("--output", "output_samplesheet.csv"),
    basename_prefix   = get_val("--basename-prefix", ""),
    check_files       = get_val("--check-files", takes_value = FALSE),
    stop_on_duplicate = get_val("--stop-on-duplicate", takes_value = FALSE)
  )
}

args <- parse_args()

# Ensure input file exists
if (is.null(args$input) || !file.exists(args$input)) {
  stop("Please provide --input <path> to an existing SampleSheet file.")
}

# -----------------------------
# 1) Read file and extract [Data] section
# -----------------------------
lines <- readr::read_lines(args$input, progress = FALSE)
data_start <- which(stringr::str_trim(lines) == "[Data]")
if (length(data_start) == 0) stop("Could not find a [Data] section in the input file.")

after_data <- lines[(data_start + 1):length(lines)]
non_empty_rel <- which(stringr::str_trim(after_data) != "")[1]
if (is.na(non_empty_rel)) stop("No header row found after [Data].")
header_idx <- data_start + non_empty_rel

# Rebuild a CSV text from header row to end of file
csv_text <- paste(lines[header_idx:length(lines)], collapse = "\n")

# -----------------------------
# 2) Parse the [Data] table
# -----------------------------
df <- readr::read_csv(file = I(csv_text), show_col_types = FALSE, progress = FALSE)

required_cols <- c("Sample_Well","Sample_Name","Sample_Plate",
                   "Sentrix_ID","Sentrix_Position","Gender")
missing <- setdiff(required_cols, names(df))
if (length(missing) > 0) {
  stop("Missing required column(s) in [Data]: ", paste(missing, collapse = ", "))
}

# Clean whitespace and ensure character types
df <- df |>
  dplyr::mutate(
    Sample_Well      = stringr::str_trim(as.character(.data$Sample_Well)),
    Sample_Name      = stringr::str_trim(as.character(.data$Sample_Name)),
    Sample_Plate     = stringr::str_trim(as.character(.data$Sample_Plate)),
    Sentrix_ID       = stringr::str_trim(as.character(.data$Sentrix_ID)),
    Sentrix_Position = stringr::str_trim(as.character(.data$Sentrix_Position)),
    Gender           = stringr::str_trim(as.character(.data$Gender))
  )

# -----------------------------
# 3) Build Basename
#    <Slide>/<Slide>_<Array>
# -----------------------------
make_basename_short <- function(slide, pos) {
  file.path(slide, paste0(slide, "_", pos))
}

out <- df |>
  dplyr::transmute(
    Sample_ID     = .data$Sample_Name,   # Adjust if you prefer a different ID
    Well          = .data$Sample_Well,
    Sample_Name   = .data$Sample_Name,
    Plate         = .data$Sample_Plate,
    Slide         = .data$Sentrix_ID,
    Array         = .data$Sentrix_Position,
    Gender        = .data$Gender,
    Basename      = make_basename_short(.data$Sentrix_ID, .data$Sentrix_Position)
  )

# -----------------------------
# 4) Detect duplicate (Slide, Array) pairs
# -----------------------------
dup_key <- paste(out$Slide, out$Array, sep = "__")
dup_tab <- table(dup_key)
if (any(dup_tab > 1)) {
  dup_vals <- names(dup_tab)[dup_tab > 1]
  msg <- paste0(
    "Duplicate Slide+Array entries detected:\n  ",
    paste(gsub("__", " / ", dup_vals), collapse = "\n  "),
    "\nEach Sentrix_ID (Slide) + Sentrix_Position (Array) must be unique."
  )
  if (isTRUE(args$stop_on_duplicate)) stop(msg) else warning(msg, call. = FALSE)
}

# -----------------------------
# 5) Optional IDAT file existence check
# -----------------------------
if (isTRUE(args$check_files)) {

  idat_ok <- vapply(seq_len(nrow(out)), function(i) {

    # Short basename (as written to output CSV)
    short_base <- out$Basename[i]

    # Full path for checking only: <prefix>/<short_base>
    long_base <- if (nzchar(args$basename_prefix)) {
      file.path(args$basename_prefix, short_base)
    } else {
      short_base
    }

    grn <- paste0(long_base, "_Grn.idat")
    red <- paste0(long_base, "_Red.idat")

    file.exists(grn) && file.exists(red)
  }, logical(1))

  if (!all(idat_ok)) {
    missing_rows <- which(!idat_ok)
    warn_list <- paste0("Row ", missing_rows, ": ", out$Basename[!idat_ok])
    warning(
      "Some IDAT pairs were not found on disk (check working dir or --basename-prefix):\n  ",
      paste(warn_list, collapse = "\n  "),
      call. = FALSE
    )
  }
}

# -----------------------------
# 6) Write output CSV
# -----------------------------
readr::write_csv(out, args$output)
message("Wrote ", nrow(out), " rows to: ", args$output)