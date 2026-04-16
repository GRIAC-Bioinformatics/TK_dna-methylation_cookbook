#!/usr/bin/env Rscript
# =============================================================================
# merge_excluded_probes.R
#
# Author: Jelle R. Dalenberg
# Created: 13-Apr-2026
#
# Description:
#   Merge multiple probe exclusion lists into one final pipeline-ready file:
#
#     excluded_probes.csv
#
#   The script:
#     - Reads one or more probe exclusion files from a directory.
#     - Supports common one-column formats produced by the QC pipeline, such as:
#         * detectionP.csv
#         * high_intensity.csv
#         * low_beadcount.csv
#         * snp_containing_SBE_CpG.csv
#         * EPICV2_probes_950K_CrossHybridization.csv
#     - Handles different header names, including:
#         * FailedProbes
#         * FlaggedProbes
#         * Flagged_Probes
#         * EPICV2_950K
#         * V1 / Probe / ProbeID / IlmnID
#     - Removes empty placeholder content such as:
#         * x
#         * X
#     - Strips quotes and whitespace from values.
#     - Concatenates all probe IDs into one combined list.
#     - Removes overlap / duplicates while preserving first-seen order.
#     - Writes a final file with:
#         * no header
#         * one probe ID per line
#
# Intended use:
#   This helper is intended for the manual probe-review checkpoint before
#   running 06_filter_flagged_probes.R, where the user decides which flagged
#   probes should be excluded from downstream processing.
#
# Checks performed:
#   Input handling:
#     --input_dir        → Directory containing probe exclusion files.
#     --files            → Optional comma-separated list of filenames to merge.
#     --output           → Path to the final merged exclusion file.
#
#   Per-file parsing:
#     - If a file does not exist, it is skipped with a warning.
#     - If a file is empty or contains only placeholder content (x / X),
#       it is treated as empty.
#     - If a file contains a recognized header and one column of probe IDs,
#       that column is extracted.
#     - If a file is not a standard CSV, raw line-based fallback parsing is used.
#
#   Output validation:
#     - Final output contains unique probe identifiers only.
#     - Final output has no header and no quotes.
#
# Usage:
#   Rscript merge_excluded_probes.R \
#       --input_dir flag/probe \
#       --output flag/excluded_probes.csv
#
#   Optional:
#       --files detectionP.csv,high_intensity.csv,low_beadcount.csv,snp_containing_SBE_CpG.csv,EPICV2_probes_950K_CrossHybridization.csv
#
# Examples:
#   Rscript merge_excluded_probes.R \
#       --input_dir flag/probe \
#       --output flag/excluded_probes.csv
#
#   Rscript merge_excluded_probes.R \
#       --input_dir . \
#       --output excluded_probes.csv \
#       --files detectionP.csv,high_intensity.csv
#
# Output:
#   A plain text / CSV-style file with one unique probe ID per line, e.g.:
#
#     cg01499815_TC11
#     cg03683899_BC11
#     cg12941754_TC11
#     cg21985559_TC11
#
# Notes:
#   - The script is robust to mixed file conventions across QC outputs.
#   - Duplicate probe IDs across files are removed automatically.
#   - Placeholder files containing only 'x' or 'X' are ignored.
#   - Designed to create the final exclusion list required by the probe
#     filtering stage of the methylation QC workflow.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
  out <- list(
    input_dir = ".",
    output = "excluded_probes.csv",
    files = c(
      "detectionP.csv",
      "high_intensity.csv",
      "low_beadcount.csv",
      "snp_containing_SBE_CpG.csv",
      "EPICV2_probes_950K_CrossHybridization.csv"
    )
  )

  i <- 1
  while (i <= length(x)) {
    key <- x[i]
    if (key == "--input_dir") {
      out$input_dir <- x[i + 1]
      i <- i + 2
    } else if (key == "--output") {
      out$output <- x[i + 1]
      i <- i + 2
    } else if (key == "--files") {
      out$files <- strsplit(x[i + 1], ",", fixed = TRUE)[[1]]
      out$files <- trimws(out$files)
      out$files <- out$files[nzchar(out$files)]
      i <- i + 2
    } else if (key %in% c("-h", "--help")) {
      cat(
"Usage:\n",
"  Rscript merge_excluded_probes.R --input_dir <dir> --output <file>\n\n",
"Optional:\n",
"  --files detectionP.csv,high_intensity.csv,low_beadcount.csv,snp_containing_SBE_CpG.csv,EPICV2_probes_950K_CrossHybridization.csv\n\n",
"Examples:\n",
"  Rscript merge_excluded_probes.R --input_dir flag/probe --output flag/excluded_probes.csv\n",
"  Rscript merge_excluded_probes.R --input_dir . --output excluded_probes.csv --files detectionP.csv,high_intensity.csv\n",
sep = "")
      quit(save = "no", status = 0)
    } else {
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    }
  }
  out
}

cfg <- parse_args(args)

clean_values <- function(x) {
  if (length(x) == 0) return(character(0))

  x <- enc2utf8(as.character(x))
  x <- gsub("^\\ufeff", "", x)          # remove BOM
  x <- gsub("\\r", "", x)               # CRLF safety
  x <- trimws(x)
  x <- gsub('^"|"$', "", x)              # strip surrounding quotes
  x <- gsub("^'|'$", "", x)
  x <- trimws(x)

  x <- x[!is.na(x)]
  x <- x[nzchar(x)]

  # Empty-placeholder files may contain only x or X
  x <- x[!tolower(x) %in% c("x")]

  # Remove common header names if they slipped into the values
  header_tokens <- c(
    "failedprobes",
    "flaggedprobes",
    "flagged_probes",
    "epicv2_950k",
    "v1",
    "probe",
    "probes",
    "probeid",
    "probe_id",
    "ilmnid",
    "targetid"
  )
  x <- x[!tolower(gsub('[[:space:]]+', '', x)) %in% header_tokens]

  x
}

looks_like_probe <- function(x) {
  grepl("^(cg|ch|rs)[A-Za-z0-9._-]+$", x, perl = TRUE)
}

read_probe_file <- function(path) {
  if (!file.exists(path)) {
    warning(sprintf("File not found, skipping: %s", path))
    return(character(0))
  }

  # Read raw lines for fallback logic
  raw_lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  raw_lines <- clean_values(raw_lines)

  # If raw file is empty or contains only placeholder/header tokens
  if (length(raw_lines) == 0) {
    message(sprintf("Skipping empty placeholder file: %s", basename(path)))
    return(character(0))
  }

  # Try structured CSV read first
  df <- tryCatch(
    read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )

  if (!is.null(df)) {
    # Header-only file, e.g. file with just x or just a header token
    if (nrow(df) == 0) {
      # Fallback: if raw line itself is a probe, keep it; otherwise empty.
      only_probe_like <- raw_lines[looks_like_probe(raw_lines)]
      if (length(only_probe_like) > 0) {
        return(unique(only_probe_like))
      }
      message(sprintf("Skipping empty/header-only file: %s", basename(path)))
      return(character(0))
    }

    preferred_cols <- c("FailedProbes", "FlaggedProbes", "Flagged_Probes", "EPICV2_950K",
                        "V1", "Probe", "Probes", "ProbeID", "Probe_ID", "IlmnID", "TargetID")

    chosen_col <- NULL
    for (nm in preferred_cols) {
      if (nm %in% names(df)) {
        chosen_col <- nm
        break
      }
    }
    if (is.null(chosen_col)) {
      chosen_col <- names(df)[1]
    }

    vals <- clean_values(df[[chosen_col]])

    # Keep only probe-like values when possible; otherwise keep cleaned non-header values.
    probe_like <- vals[looks_like_probe(vals)]
    if (length(probe_like) > 0) {
      return(unique(probe_like))
    }

    if (length(vals) > 0) {
      return(unique(vals))
    }
  }

  # Fallback for strange one-item-per-line files
  raw_probe_like <- raw_lines[looks_like_probe(raw_lines)]
  if (length(raw_probe_like) > 0) {
    return(unique(raw_probe_like))
  }

  character(0)
}

# Build full paths
input_paths <- file.path(cfg$input_dir, cfg$files)

# Read and concatenate while preserving first-seen order
all_probes <- character(0)
per_file_counts <- integer(length(input_paths))
names(per_file_counts) <- basename(input_paths)

for (i in seq_along(input_paths)) {
  p <- input_paths[i]
  probes <- read_probe_file(p)
  per_file_counts[i] <- length(probes)
  if (length(probes) > 0) {
    all_probes <- c(all_probes, probes)
  }
}

final_probes <- unique(all_probes)

# Ensure output directory exists
out_dir <- dirname(cfg$output)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

write.table(
  final_probes,
  file = cfg$output,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

message("----------------------------------------")
message("Probe merge summary")
message("----------------------------------------")

name_width  <- max(nchar(names(per_file_counts)))
count_width <- max(nchar(as.character(per_file_counts)), nchar(as.character(length(final_probes))))

for (nm in names(per_file_counts)) {
  message(sprintf(
    paste0("%-", name_width, "s  %", count_width, "d"),
    nm,
    per_file_counts[[nm]]
  ))
}

message("----------------------------------------")
message(sprintf(
  paste0("%-", name_width, "s  %", count_width, "d"),
  "Unique probes written:",
  length(final_probes)
))
message(sprintf("Output file: %s", cfg$output))
message("----------------------------------------")
