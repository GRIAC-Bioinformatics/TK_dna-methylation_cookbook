# filter flagged probe 
# Rscript 05_filter_flagged_probes.R --probe_qc_file probe_qc_overview.csv --min_flag_overlap 2 --mset "MSet_sample_filtered.RData" --grset "grSet_sample_filtered.RData" --ratioset "ratioSet_sample_filtered.RData" --rgset "RGset_sample_filtered.RData" --rgsetext "RGsetEXT_sample_filtered.RData" --base_suffix "_probe_filtered"

suppressPackageStartupMessages({
  tryCatch({
    library(minfi)
    library(optparse)
    }, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

# Parse command line arguments
option_list <- list(
  make_option(c("-s", "--probe_qc_file"), type = "character", default = NULL,
              help = "Path to the probe QC CSV file"),
  make_option(c("-m", "--mset"), type = "character", default = NULL,
              help = "Path to the MSet.RData file"),
  make_option(c("-g", "--grset"), type = "character", default = NULL,
              help = "Path to the grSet.RData file"),
  make_option(c("-r", "--ratioset"), type = "character", default = NULL,
              help = "Path to the ratioSet.RData file"),
  make_option(c("-R", "--rgset"), type = "character", default = NULL,
              help = "Path to the RGset.RData file"),
  make_option(c("-E", "--rgsetext"), type = "character", default = NULL,
              help = "Path to the RGsetEXT.RData file"),
  make_option(c("-b", "--base_suffix"), type = "character", default = "_probe_filtered",
              help = "Suffix to add to output file names (e.g., '_filtered'). Default is '_probe_filtered'.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# CSV file with flagged probes
failed.probes <- read.csv(opt$probe_qc_file, stringsAsFactors= FALSE, header = TRUE)
failed.probes <- failed.probes$X
message(paste("Number of failed probes identified:", length(failed.probes)))

# --- Function to load, filter, and save an RData object ---
process_minfi_object <- function(file_path, object_name, failed.probes, base_suffix) {
  if (!is.null(file_path) && file.exists(file_path)) {
    message(paste("Loading", object_name, "from:", file_path))
    # Load the R data object. It will create a variable with 'object_name' (e.g., 'MSet').
    load(file_path)

    # Get the loaded object dynamically by its name.
    current_object <- get(object_name)

    # Filter the object by removing failed probes
    initial_probe_count <- nrow(current_object)
    print(head(rownames(current_object)))
    print(head(failed.probes))
    probes_to_remove_indices <- which(rownames(current_object) %in% failed.probes)
    current_object_filtered <- current_object[-probes_to_remove_indices,]
    
    message(paste("Filtered", object_name, ": Removed", 
                  initial_probe_count - nrow(current_object_filtered), 
                  "probes out of", initial_probe_count))

    # Subset the object to keep only the desired samples.
    # Assign it directly back to the original object_name variable.
    assign(object_name, current_object_filtered)

    # Construct the output file name.
    base_file_name <- tools::file_path_sans_ext(basename(file_path))
    output_file_name <- file.path(dirname(file_path), paste0(base_file_name, base_suffix, ".RData"))

    # Save the *modified* object (which has the original 'object_name' variable name).
    save(list = object_name, file = output_file_name)
    message(paste("  Filtered and saved:", output_file_name))

    message(paste("Filtered and saved:", output_file_name))
  } else if (!is.null(file_path)) {
    message(paste("Warning: File not found for", object_name, ":", file_path, ". Skipping filtering for this object."))
  } else {
    message(paste("No file path provided for", object_name, ". Skipping filtering for this object."))
  }
}


# --- Function to load, filter, and save an RData object ( RGset or RGsetEXT) ---
process_minfi_RGobject <- function(file_path, object_name, failed.probes, base_suffix) {
  if (!is.null(file_path) && file.exists(file_path)) {
    message(paste("Loading", object_name, "from:", file_path))
    # Load the R data object. It will create a variable with 'object_name' (e.g., 'MSet').
    load(file_path)

    # Get the loaded object dynamically by its name.
    current_object <- get(object_name)

    # Filter the object by removing failed probes
    initial_probe_count <- nrow(current_object)
    current_object_filtered = subsetByLoci(current_object, includeLoci = NULL, excludeLoci = failed.probes,
             keepControls = TRUE, keepSnps = TRUE)
    
    message(paste("Filtered", object_name, ": Removed", 
                  initial_probe_count - nrow(current_object_filtered), 
                  "probes out of", initial_probe_count))

    # Subset the object to keep only the desired samples.
    # Assign it directly back to the original object_name variable.
    assign(object_name, current_object_filtered)

    # Construct the output file name.
    base_file_name <- tools::file_path_sans_ext(basename(file_path))
    output_file_name <- file.path(dirname(file_path), paste0(base_file_name, base_suffix, ".RData"))

    # Save the *modified* object (which has the original 'object_name' variable name).
    save(list = object_name, file = output_file_name)
    message(paste("  Filtered and saved:", output_file_name))

    message(paste("Filtered and saved:", output_file_name))
  } else if (!is.null(file_path)) {
    message(paste("Warning: File not found for", object_name, ":", file_path, ". Skipping filtering for this object."))
  } else {
    message(paste("No file path provided for", object_name, ". Skipping filtering for this object."))
  }
}


# --- Process each specified RData file ---

# Process MSet
if (!is.null(opt$mset)) {
  message("MSet file provided: ", opt$mset)
  process_minfi_object(opt$mset, "MSet", failed.probes, opt$base_suffix)
} else {
  message("No MSet file provided. Skipping MSet processing.")
}

# Process grSet
if (!is.null(opt$grset)) {
  message("grSet file provided: ", opt$grset)
  process_minfi_object(opt$grset, "grSet", failed.probes, opt$base_suffix)
} else {
  message("No grSet file provided. Skipping grSet processing.")
}

# Process ratioSet
if (!is.null(opt$ratioset)) {
  message("ratioSet file provided: ", opt$ratioset)
  process_minfi_object(opt$ratioset, "ratioSet", failed.probes, opt$base_suffix)
} else {
  message("No ratioSet file provided. Skipping ratioSet processing.")
}

# For RGset and RGsetEXT, you cannot simply subset using subsetByLoci
# Process RGset
if (!is.null(opt$rgset)) {
  message("RGset file provided: ", opt$rgset)
  process_minfi_RGobject(opt$rgset, "RGset", failed.probes, opt$base_suffix)
} else {
  message("No RGset file provided. Skipping RGset processing.")
}

# Process RGsetEXT
if (!is.null(opt$rgsetext)) {
  message("RGsetEXT file provided: ", opt$rgsetext)
  process_minfi_RGobject(opt$rgsetext, "RGsetEXT", failed.probes, opt$base_suffix)
} else {
  message("No RGsetEXT file provided. Skipping RGsetEXT processing.")
}

message("Script execution finished.")
