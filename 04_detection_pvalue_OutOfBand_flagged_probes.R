# After filtering the samples, we are going to look at pOOBAH to remove probes
# Rscript 04_detection_pvalue_OutOfBand_flagged_probes.R --sigsetlist SigSetList.RData --cutoff 0.01 --threshold 0.99 --flagged flagged_probes.csv --pdf detectionP_OOB_flagged_probes.pdf


suppressPackageStartupMessages({
  tryCatch({
  library(sesame)
  library(optparse)
  library(BiocParallel)
  library(ggplot2)
  library(matrixStats)
}, error = function(e) {
    stop("Package loading failed: ", e$message, call. = FALSE)
  })
})

# Define command-line options
option_list = list(
  make_option(c("-s", "--sigsetlist"), type="character", default=NULL, 
              help="base directory containing IDAT files (e.g., /groups/umcg-griac/tmp02/rawdata/medall/IDAT)", metavar="character"),
  make_option(c("-o", "--cutoff"), type="numeric", default=0.01,
                help="Cutoff for detection p-value (default: 0.05)", metavar="NUM"),
  make_option(c("-t", "--threshold"), type="numeric", default=0.99,
                help="Fraction of samples in which the probe should pass the p value cut off, for a probe to pass detectionP filter", metavar="NUM"),
  make_option(c("-f", "--flagged"), type="character", default=NULL,
                help="CSV file to save flagged probes which have high detection p-values", metavar="FILE"),
  make_option(c("-p", "--pdf"), type="character", default=NULL,
                help="Output PDF file path", metavar="FILE")
)

# Parse command-line arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Check if required arguments are provided
if (is.null(opt$sigsetlist)){
  print_help(opt_parser)
  stop("SigSetList Rdata file must be supplied (-s/--sigsetlist)", call.=FALSE)
}
if (is.null(opt$flagged)){
  print_help(opt_parser)
  stop("Output flagged CSV file argument must be supplied (-f/--flagged)", call.=FALSE)
}
if (is.null(opt$pdf)){
  print_help(opt_parser)
  stop("Output PDF file argument must be supplied (-p/--pdf)", call.=FALSE)
}
tryCatch({

    cutoff <- as.numeric(opt$cutoff)
    threshold <- as.numeric(opt$threshold)

    message("Loading SigSetList object...")
    load(opt$sigsetlist)

    message("SigSetList loaded successfully. Number of samples in SigSets: ", length(SigSetList))

    # Initialize an empty list to store the results
    pSigSetList <- vector("list", length(SigSetList))
    names(pSigSetList) <- names(SigSetList)

    # Loop through each element of SigSetList
    for (i in names(SigSetList)) {
        pfx <- SigSetList[[i]] # Get the current sample
        message("Processing sample: ", i)
        # Use a tryCatch block to gracefully handle errors for individual samples
        tryCatch({
            # Call pOOBAH for the current sample
            pSigSetList[[i]] <- pOOBAH(pfx,return.pval = TRUE,combine.neg = TRUE,pval.threshold = cutoff,verbose = TRUE) 

        }, error = function(e) {
            # If an error occurs, store the error message in the list
            message("Error processing sample ", i, ": ", e$message)
        })
    }

    pSigSetList <- data.frame(pSigSetList, check.names = FALSE)
    PassProbeDF <- data.frame("Fraction" = rowSums(pSigSetList < cutoff)/ncol(pSigSetList))
    failed.probes <- data.frame("FailedProbes" = names(which(rowSums(pSigSetList < cutoff)/ncol(pSigSetList) < threshold)))

    # Open PDF for plotting
    pdf(opt$pdf, width = 10, height = 7)
    message("Generating detection p-value plots...")
    print(ggplot(PassProbeDF, aes(x = Fraction)) +
        geom_histogram(binwidth = 0.001, fill = "blue", color = "black") +
        labs(title = "Distribution of number of samples in which each probes passes the detection P-value cutoff",
            subtitle = paste0("Cutoff: ", cutoff, " (red line indicates threshold)\n",
                                "Probes passing the cutoff in at least ", threshold * 100, "% of samples are retained"),
            x = "Fraction of Samples Passing Cutoff",
            y = "Number of Probes") +
        theme_minimal() +
        scale_x_log10() + # Added log scale for x-axis
        geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
        annotate("text", x = threshold, y = max(table(PassProbeDF$Fraction)),
                label = paste0("Cutoff ",threshold), color = "red", vjust = -1))

    dev.off()

    message("Filtering probes based on detection p-values...")
    
    # Create a data frame for flagged samples and save to CSV
    if (dim(failed.probes)[1] > 0) {
        write.csv(failed.probes, file = opt$flagged, row.names = FALSE)
        message(paste0("Flagged ", dim(failed.probes)[1], " probes with p-value < ",cutoff," in ",threshold*100,"% samples saved to ", opt$flagged))
    } else {
        message("No probes flagged with p-value < ", cutoff, " in ", threshold * 100, "% of samples.")
        # Ensure an empty file is still created if no samples are flagged
        write.csv(data.frame("FailedProbes" = character()), file = opt$flagged, row.names = FALSE)
    }
    # Close the PDF device
    message("Detection p-value plots generated successfully and saved to ", opt$pdf)

        }, error = function(e) {    
    # Handle errors and stop execution with an informative message
    stop("Failed to process IDAT files or save output file: ", e$message)
    })
