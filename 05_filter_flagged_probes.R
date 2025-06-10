# filter flagged probe 
# Rscript 05_filter_flagged_probes.R --probe_qc_file probe_qc_overview.csv --min_flag_overlap 2 --mset "MSet_sample_filtered.RData" --grset "grSet_sample_filtered.RData" --ratioset "ratioSet_sample_filtered.RData" --rgset "RGset_sample_filtered.RData" --rgsetext "RGsetEXT_sample_filtered.RData" --base_suffix "_probe_filtered"

library(optparse)
library(minfi) 

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
              help = "Suffix to add to output file names (e.g., '_filtered'). Default is '_probe_filtered'."),
  make_option(c("-f", "--min_flag_overlap"), type = "character", default = NULL,
              help = "Minimum number of overlapping flags required to filter a probe (required).")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# CSV file with flagged probes
failed.probes <- read.csv(opt$probe_qc_file, stringsAsFactors= FALSE, header = TRUE)
failed.probes <- failed.probes[failed.probes$Total_Occurrence,]$X
message(paste("Number of failed probes identified:", length(failed.probes)))

# --- Function to load, filter, and save an RData object ---
process_minfi_object <- function(file_path, object_name, failed_probes_list, suffix) {
  if (!is.null(file_path) && file.exists(file_path)) {
    message(paste("Loading", object_name, "from:", file_path))
    # Load the object. It will be named 'filtered_object' due to how it was likely saved previously.
    load(file_path)
    
    # Assign the loaded object to its intended name
    current_object <- get("filtered_object")
    rm(filtered_object) # Clean up the temporary loaded name

    # Filter the object by removing failed probes
    initial_probe_count <- nrow(current_object)
    probes_to_remove_indices <- which(rownames(current_object) %in% failed_probes_list)
    current_object_filtered <- current_object[-probes_to_remove_indices,]
    
    message(paste("Filtered", object_name, ": Removed", 
                  initial_probe_count - nrow(current_object_filtered), 
                  "probes out of", initial_probe_count))

    # Save the filtered object to the new R data file
    base_file_name <- tools::file_path_sans_ext(basename(file_path))
    output_file_name <- file.path(dirname(file_path), paste0(base_file_name, suffix, ".RData"))
    
    # Use assign to save with the correct variable name
    assign(object_name, current_object_filtered)
    save(list = object_name, file = output_file_name)
    message(paste("Filtered and saved:", output_file_name))
  } else if (!is.null(file_path)) {
    message(paste("Warning: File not found for", object_name, ":", file_path, ". Skipping filtering for this object."))
  } else {
    message(paste("No file path provided for", object_name, ". Skipping filtering for this object."))
  }
}

# --- Process each specified RData file ---
# Process grSet
process_minfi_object(opt$mset, "MSet_probe_filtered", failed.probes, opt$base_suffix)

# Process grSet
process_minfi_object(opt$grset, "grSet_probe_filtered", failed.probes, opt$base_suffix)

# Process ratioSet
process_minfi_object(opt$ratioset, "ratioSet_probe_filtered", failed.probes, opt$base_suffix)

# Process RGset
process_minfi_object(opt$rgset, "RGset_probe_filtered", failed.probes, opt$base_suffix)

# Process RGsetEXT
process_minfi_object(opt$rgsetext, "RGsetEXT_probe_filtered", failed.probes, opt$base_suffix)

message("Script execution finished.")
