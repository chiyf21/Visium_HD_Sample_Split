#' Merge Visium HD Sample Split Results into a Seurat Object
#' ============================================================
#'
#' Usage:
#'   Rscript merge_to_seurat.R --data-dir binned_outputs/square_016um \
#'                             --assignment spot_sample_assignment.csv \
#'                             --out-dir ./seurat_output
#'
#' Or interactively:
#'   source("merge_to_seurat.R")
#'   obj <- merge_sample_assignment(
#'       data.dir  = "binned_outputs/square_016um",
#'       csv.path  = "spot_sample_assignment.csv"
#'   )

library(Seurat)

#' Read CSV assignment and merge into a Seurat object
#'
#' @param data.dir  Path to Space Ranger binned output folder
#' @param csv.path  Path to spot_sample_assignment.csv
#' @param bin.size  Bin size in microns (default: 16)
#' @return A Seurat object with sample_id / sample_name / in_tissue metadata
merge_sample_assignment <- function(
    data.dir,
    csv.path,
    bin.size = 16
) {
    # --- Load original Visium HD data ---
    message("Loading Visium HD data from: ", data.dir)
    obj <- Load10X_Spatial(
        data.dir  = data.dir,
        bin.size  = bin.size
    )

    # --- Read assignment ---
    message("Reading assignment: ", csv.path)
    assignment <- read.csv(csv.path, row.names = 1, stringsAsFactors = FALSE)

    # --- Check barcode overlap ---
    common <- intersect(Cells(obj), rownames(assignment))
    message(sprintf("Barcodes in Seurat:  %d", ncol(obj)))
    message(sprintf("Barcodes in CSV:     %d", nrow(assignment)))
    message(sprintf("Common barcodes:     %d", length(common)))

    if (length(common) == 0) {
        stop("No barcode overlap. Check that the CSV matches the Visium HD dataset.")
    }

    # --- Merge ---
    obj <- AddMetaData(obj, metadata = assignment)

    # --- Report ---
    cat("\nSample distribution:\n")
    print(table(obj$sample_name))

    invisible(obj)
}


#' Split a Seurat object by sample_id into a list
split_by_sample <- function(obj) {
    samples <- setdiff(unique(obj$sample_id), 0)
    result <- list()
    for (sid in sort(samples)) {
        nm <- sprintf("sample_%03d", sid)
        result[[nm]] <- subset(obj, subset = sample_id == sid)
    }
    result
}


# ============================================================
# CLI entry point
# ============================================================
if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)

    # Defaults
    data.dir   <- "binned_outputs/square_016um"
    csv.path   <- "spot_sample_assignment.csv"
    out.dir    <- "seurat_output"
    bin.size   <- 16

    # Parse
    i <- 1
    while (i <= length(args)) {
        switch(args[i],
            "--data-dir"   = { data.dir <- args[i+1]; i <- i+2 },
            "--assignment" = { csv.path <- args[i+1]; i <- i+2 },
            "--out-dir"    = { out.dir  <- args[i+1]; i <- i+2 },
            "--bin-size"   = { bin.size <- as.numeric(args[i+1]); i <- i+2 },
            { stop("Unknown argument: ", args[i]) }
        )
    }

    message("===== Visium HD Sample Split → Seurat =====")
    message("  data.dir:   ", data.dir)
    message("  csv.path:   ", csv.path)
    message("  bin.size:   ", bin.size)
    message("  out.dir:    ", out.dir)

    # Run
    obj <- merge_sample_assignment(data.dir, csv.path, bin.size)

    # Save
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

    # Full object
    saveRDS(obj, file.path(out.dir, "visium_hd_sample_split.rds"))
    message("Saved: ", file.path(out.dir, "visium_hd_sample_split.rds"))

    # Per-sample objects
    sample_list <- split_by_sample(obj)
    for (nm in names(sample_list)) {
        saveRDS(sample_list[[nm]], file.path(out.dir, paste0(nm, ".rds")))
    }
    message("Saved ", length(sample_list), " per-sample RDS files.")

    # Metadata CSV
    write.csv(obj@meta.data, file.path(out.dir, "metadata_with_sample.csv"))
    message("Saved: metadata_with_sample.csv")
    message("===== Done =====")
}
