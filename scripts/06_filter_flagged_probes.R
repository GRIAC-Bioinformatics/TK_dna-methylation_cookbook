#!/usr/bin/env Rscript
# =============================================================================
# 06_filter_flagged_probes.R
# Author: Vartika Bisht, Jelle R. Dalenberg
# Updated: 14-Apr-2026
#
# Description:
#   Filter flagged probes from one or more minfi RData objects and generate a
#   compact PDF summary showing the effect of probe removal.
#
#   Features:
#     - supports --plot_only to skip filtering when filtered objects already
#       exist, and only (re)generate the PDF summary.
#     - PDF pages are intentionally compact:
#         1) retained vs removed probes per object
#         2) summary table
#         3) per-sample median beta before vs after filtering (grSet preferred,
#            MSet fallback; only if available)
# =============================================================================

suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(optparse)
  }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

option_list <- list(
  make_option(c("-s", "--flagged"), type = "character", default = NULL,
              help = "Path to a headerless file or CSV containing probes to be excluded."),
  make_option(c("-m", "--mset"), type = "character", default = NULL,
              help = "Path to the MSet.RData file (optional)."),
  make_option(c("-g", "--grset"), type = "character", default = NULL,
              help = "Path to the grSet.RData file (optional)."),
  make_option(c("-R", "--rgset"), type = "character", default = NULL,
              help = "Path to the RGset.RData file (optional)."),
  make_option(c("-x", "--rgsetext"), type = "character", default = NULL,
              help = "Path to the RGsetEXT.RData file (optional)."),
  make_option(c("-b", "--base_suffix"), type = "character", default = "_probe_filtered",
              help = "Suffix to add to output file names. Default: '_probe_filtered'."),
  make_option(c("-p", "--pdf_output"), type = "character", default = NULL,
              help = "Path to the output PDF report (optional but recommended)."),
  make_option(c("--beta_plot_max_loci"), type = "integer", default = 100000,
              help = "Maximum number of loci to use for the median-beta scatter plot. Default: 100000."),
  make_option(c("--plot_only"), action = "store_true", default = FALSE,
              help = "Skip filtering and only (re)generate the PDF using already filtered output files inferred from --base_suffix.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
if (is.null(opt$flagged)) {
  stop("Error: Required argument --flagged is missing. Use --help for more information.", call. = FALSE)
}

if (all(sapply(list(opt$mset, opt$grset, opt$rgset, opt$rgsetext), is.null))) {
  stop("Error: At least one of --mset, --grset, --rgset, or --rgsetext must be provided.", call. = FALSE)
}

if (!file.exists(opt$flagged)) {
  stop("The specified flagged probe file does not exist: ", opt$flagged, call. = FALSE)
}

input_rdata_paths <- list(
  MSet = opt$mset,
  grSet = opt$grset,
  RGset = opt$rgset,
  RGsetEXT = opt$rgsetext
)
input_rdata_paths <- input_rdata_paths[!sapply(input_rdata_paths, is.null)]

for (file_path in input_rdata_paths) {
  if (!file.exists(file_path)) {
    stop("The specified RData file does not exist: ", file_path, call. = FALSE)
  }
}

if (opt$plot_only && is.null(opt$pdf_output)) {
  stop("Error: --plot_only requires --pdf_output because the purpose is to re-generate the PDF.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Read and clean failed probes
# -----------------------------------------------------------------------------
read_failed_probes <- function(path) {
  raw <- read.csv(path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)[[1]]
  raw <- enc2utf8(as.character(raw))
  raw <- gsub("^\\ufeff", "", raw)
  raw <- trimws(raw)
  raw <- gsub('^"|"$', "", raw)
  raw <- gsub("^'|'$", "", raw)
  raw <- trimws(raw)
  raw <- raw[!is.na(raw) & nzchar(raw)]
  unique(raw)
}

failed.probes <- read_failed_probes(opt$flagged)
message(length(failed.probes), " unique probes in the exclusion list.")
if (opt$plot_only) {
  message("Running in --plot_only mode: filtering will be skipped and the PDF will be regenerated from existing filtered objects.")
}

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
get_output_file_name <- function(file_path, base_suffix) {
  base_file_name <- tools::file_path_sans_ext(basename(file_path))
  file.path(dirname(file_path), paste0(base_file_name, base_suffix, ".RData"))
}

load_rdata_object <- function(file_path, object_name) {
  e <- new.env(parent = emptyenv())
  load(file_path, envir = e)
  if (!exists(object_name, envir = e, inherits = FALSE)) {
    stop("Expected object '", object_name, "' was not found in ", file_path, call. = FALSE)
  }
  get(object_name, envir = e)
}

save_rdata_object <- function(file_path, object_name, object_value) {
  e <- new.env(parent = emptyenv())
  assign(object_name, object_value, envir = e)
  save(list = object_name, file = file_path, envir = e)
}

sample_loci_for_beta <- function(obj, max_loci = 100000) {
  n <- nrow(obj)
  if (n <= max_loci) {
    rownames(obj)
  } else {
    set.seed(1)
    sample(rownames(obj), size = max_loci, replace = FALSE)
  }
}

extract_beta_scatter_data <- function(before_obj, after_obj, object_name, max_loci = 100000) {
  loci_before <- sample_loci_for_beta(before_obj, max_loci = max_loci)
  loci_after <- intersect(loci_before, rownames(after_obj))

  if (length(loci_after) < 1000) {
    return(NULL)
  }

  beta_before <- getBeta(before_obj[loci_before, ])
  beta_after  <- getBeta(after_obj[loci_after, ])

  common_samples <- intersect(colnames(beta_before), colnames(beta_after))
  beta_before <- beta_before[, common_samples, drop = FALSE]
  beta_after  <- beta_after[, common_samples, drop = FALSE]

  beta_before_common <- beta_before[rownames(beta_before) %in% loci_after, , drop = FALSE]
  if (nrow(beta_before_common) == 0 || nrow(beta_after) == 0) {
    return(NULL)
  }

  med_before <- apply(beta_before_common, 2, median, na.rm = TRUE)
  med_after  <- apply(beta_after, 2, median, na.rm = TRUE)

  list(
    object_name = object_name,
    loci_before = length(loci_before),
    loci_after = length(loci_after),
    med_before = med_before,
    med_after = med_after
  )
}

process_minfi_object <- function(file_path, object_name, failed.probes, base_suffix,
                                 beta_capture = FALSE, beta_plot_max_loci = 100000,
                                 plot_only = FALSE) {
  if (is.null(file_path)) {
    message("No file path provided for ", object_name, ". Skipping.")
    return(list(summary = NULL, beta_plot = NULL))
  }
  if (!file.exists(file_path)) {
    message("Warning: File not found for ", object_name, ": ", file_path, ". Skipping.")
    return(list(summary = NULL, beta_plot = NULL))
  }

  output_file_name <- get_output_file_name(file_path, base_suffix)
  current_object <- load_rdata_object(file_path, object_name)
  initial_probe_count <- nrow(current_object)

  if (plot_only) {
    if (!file.exists(output_file_name)) {
      stop("In --plot_only mode, the filtered file does not exist for ", object_name, ": ", output_file_name, call. = FALSE)
    }
    message("Loading existing filtered ", object_name, " from: ", output_file_name)
    current_object_filtered <- load_rdata_object(output_file_name, object_name)
  } else {
    message("Loading ", object_name, " from: ", file_path)
    probes_to_remove_indices <- which(rownames(current_object) %in% failed.probes)
    current_object_filtered <- current_object[-probes_to_remove_indices, ]
    save_rdata_object(output_file_name, object_name, current_object_filtered)
    message("Filtered and saved: ", output_file_name)
  }

  removed_probe_count <- initial_probe_count - nrow(current_object_filtered)
  percent_removed <- if (initial_probe_count > 0) round(100 * removed_probe_count / initial_probe_count, 2) else 0

  message("Compared ", object_name, ": Removed ", removed_probe_count,
          " probes out of ", initial_probe_count,
          sprintf(" (%.2f%%)", percent_removed))

  beta_plot <- NULL
  if (beta_capture) {
    beta_plot <- tryCatch(
      extract_beta_scatter_data(current_object, current_object_filtered, object_name,
                                max_loci = beta_plot_max_loci),
      error = function(err) {
        message("Warning: Median-beta scatter data could not be generated for ", object_name, ": ", err$message)
        NULL
      }
    )
  }

  list(
    summary = data.frame(
      object_name = object_name,
      initial_probe_count = initial_probe_count,
      final_probe_count = nrow(current_object_filtered),
      removed_probe_count = removed_probe_count,
      percent_removed = percent_removed,
      output_file = output_file_name,
      stringsAsFactors = FALSE
    ),
    beta_plot = beta_plot
  )
}

process_minfi_RGobject <- function(file_path, object_name, failed.probes, base_suffix,
                                   plot_only = FALSE) {
  if (is.null(file_path)) {
    message("No file path provided for ", object_name, ". Skipping.")
    return(list(summary = NULL))
  }
  if (!file.exists(file_path)) {
    message("Warning: File not found for ", object_name, ": ", file_path, ". Skipping.")
    return(list(summary = NULL))
  }

  output_file_name <- get_output_file_name(file_path, base_suffix)
  current_object <- load_rdata_object(file_path, object_name)
  initial_probe_count <- nrow(current_object)

  if (plot_only) {
    if (!file.exists(output_file_name)) {
      stop("In --plot_only mode, the filtered file does not exist for ", object_name, ": ", output_file_name, call. = FALSE)
    }
    message("Loading existing filtered ", object_name, " from: ", output_file_name)
    current_object_filtered <- load_rdata_object(output_file_name, object_name)
  } else {
    message("Loading ", object_name, " from: ", file_path)
    current_object_filtered <- subsetByLoci(
      current_object,
      includeLoci = NULL,
      excludeLoci = failed.probes,
      keepControls = TRUE,
      keepSnps = TRUE
    )
    save_rdata_object(output_file_name, object_name, current_object_filtered)
    message("Filtered and saved: ", output_file_name)
  }

  removed_probe_count <- initial_probe_count - nrow(current_object_filtered)
  percent_removed <- if (initial_probe_count > 0) round(100 * removed_probe_count / initial_probe_count, 2) else 0

  message("Compared ", object_name, ": Removed ", removed_probe_count,
          " probes out of ", initial_probe_count,
          sprintf(" (%.2f%%)", percent_removed))

  list(
    summary = data.frame(
      object_name = object_name,
      initial_probe_count = initial_probe_count,
      final_probe_count = nrow(current_object_filtered),
      removed_probe_count = removed_probe_count,
      percent_removed = percent_removed,
      output_file = output_file_name,
      stringsAsFactors = FALSE
    )
  )
}

# -----------------------------------------------------------------------------
# Process objects
# -----------------------------------------------------------------------------
summary_list <- list()
beta_plot_info <- NULL

res_grset <- process_minfi_object(opt$grset, "grSet", failed.probes, opt$base_suffix,
                                  beta_capture = TRUE,
                                  beta_plot_max_loci = opt$beta_plot_max_loci,
                                  plot_only = opt$plot_only)
summary_list[[length(summary_list) + 1]] <- res_grset$summary
if (!is.null(res_grset$beta_plot)) {
  beta_plot_info <- res_grset$beta_plot
}

res_mset <- process_minfi_object(opt$mset, "MSet", failed.probes, opt$base_suffix,
                                 beta_capture = is.null(beta_plot_info),
                                 beta_plot_max_loci = opt$beta_plot_max_loci,
                                 plot_only = opt$plot_only)
summary_list[[length(summary_list) + 1]] <- res_mset$summary
if (is.null(beta_plot_info) && !is.null(res_mset$beta_plot)) {
  beta_plot_info <- res_mset$beta_plot
}

res_rgset <- process_minfi_RGobject(opt$rgset, "RGset", failed.probes, opt$base_suffix,
                                    plot_only = opt$plot_only)
summary_list[[length(summary_list) + 1]] <- res_rgset$summary

res_rgsetext <- process_minfi_RGobject(opt$rgsetext, "RGsetEXT", failed.probes, opt$base_suffix,
                                       plot_only = opt$plot_only)
summary_list[[length(summary_list) + 1]] <- res_rgsetext$summary

summary_df <- do.call(rbind, Filter(Negate(is.null), summary_list))

if (is.null(summary_df) || nrow(summary_df) == 0) {
  stop("No objects were processed successfully.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# PDF report
# -----------------------------------------------------------------------------
if (!is.null(opt$pdf_output)) {
  pdf(opt$pdf_output, width = 11, height = 8.5)

  # Page 1: retained vs removed probes
  par(mfrow = c(1, 1), mar = c(5, 8, 4, 2) + 0.1)
  plot_matrix <- rbind(
    Retained = summary_df$final_probe_count,
    Removed  = summary_df$removed_probe_count
  )
  colnames(plot_matrix) <- summary_df$object_name

  barplot(
    plot_matrix,
    horiz = TRUE,
    las = 1,
    col = c("steelblue", "firebrick"),
    main = if (opt$plot_only) "Effect of probe filtering per object (plot-only mode)" else "Effect of probe filtering per object",
    xlab = "Number of probes",
    legend.text = TRUE,
    args.legend = list(x = "bottomright", bty = "n")
  )
  mtext(sprintf("Input exclusion list: %d unique probes", length(failed.probes)), side = 3, line = 0.5, cex = 0.9)

  # Page 2: summary table only
  par(mfrow = c(1, 1), mar = c(1, 1, 2, 1))
  plot.new()
  title("Probe filtering summary", line = -1)
  summary_print <- summary_df[, c("object_name", "initial_probe_count", "final_probe_count", "removed_probe_count")]
  summary_print$initial_probe_count <- format(summary_print$initial_probe_count, big.mark = ",")
  summary_print$final_probe_count   <- format(summary_print$final_probe_count, big.mark = ",")
  summary_print$removed_probe_count <- format(summary_print$removed_probe_count, big.mark = ",")
  txt <- capture.output(print(summary_print, row.names = FALSE, right = TRUE))
  text(x = 0.02, y = 0.95, labels = paste(txt, collapse = "\n"), adj = c(0, 1), family = "mono", cex = 1.0)
  mtext(paste("Filtered objects suffix:", opt$base_suffix), side = 1, line = -1, cex = 0.9)

  # Page 3: per-sample median beta before vs after, if possible
  if (!is.null(beta_plot_info)) {
    par(mfrow = c(1, 1), mar = c(5, 5, 4, 2) + 0.1)
    x <- beta_plot_info$med_before
    y <- beta_plot_info$med_after
    lims <- range(c(x, y), na.rm = TRUE)
    plot(
      x, y,
      pch = 19,
      col = "navy",
      xlab = "Median beta before filtering",
      ylab = "Median beta after filtering",
      main = paste0(beta_plot_info$object_name, ": Sample median beta before vs after"),
      xlim = lims,
      ylim = lims
    )
    abline(0, 1, col = "firebrick", lty = 2, lwd = 2)
    mtext("Points near the diagonal indicate little sample-level distribution shift.", side = 3, line = 0.5, cex = 0.85)
    if (!is.null(names(x)) && length(names(x)) == length(x) && length(x) <= 30) {
      text(x, y, labels = names(x), pos = 3, cex = 0.7)
    }
  }

  dev.off()
  message("PDF summary written to: ", opt$pdf_output)
}

if (opt$plot_only) {
  message("\nPlot-only run completed successfully.")
} else {
  message("\nProbe filtering process completed for all specified R data files.")
}
