#!/usr/bin/env Rscript

rm(list = ls())

# Optional: set your working directory
setwd("/Volumes/Public/SNP_array_data/horse_data/ref/manifest_files")

# ----------------------------- #
# Helper: wrap FASTA sequences
# ----------------------------- #
wrap_fasta <- function(seq, width = 60) {
  seq <- gsub("\\s+", "", seq)
  n <- nchar(seq)
  if (is.na(seq) || n == 0) return("")
  starts <- seq(1, n, by = width)
  ends <- pmin(starts + width - 1, n)
  paste(mapply(function(s, e) substr(seq, s, e), starts, ends), collapse = "\n")
}

# ----------------------------------------------------------- #
# Read Illumina manifest and return only the SNP assay table
# ----------------------------------------------------------- #
read_illumina_manifest_assay <- function(input_file) {
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }
  
  lines <- readLines(input_file, warn = FALSE)
  
  # Trim only for section detection; keep original lines for CSV parsing
  trimmed <- trimws(lines)
  
  # Locate [Assay]
  assay_idx <- which(trimmed == "[Assay]")
  
  if (length(assay_idx) == 0) {
    # Fallback: find the first line containing the actual column header
    header_idx <- which(grepl("(^|,)AlleleA_ProbeSeq(,|$)", lines))
    if (length(header_idx) == 0) {
      stop("Could not find [Assay] section or a header line containing 'AlleleA_ProbeSeq'.")
    }
    start_idx <- header_idx[1]
  } else {
    start_idx <- assay_idx[1] + 1
  }
  
  # Find the next section header after [Assay], e.g. [Controls]
  # Important: only lines that START with [SectionName]
  next_section_idx <- which(grepl("^\\[[^]]+\\]$", trimmed) & seq_along(trimmed) > start_idx)
  
  if (length(next_section_idx) > 0) {
    end_idx <- next_section_idx[1] - 1
  } else {
    end_idx <- length(lines)
  }
  
  assay_lines <- lines[start_idx:end_idx]
  assay_lines <- assay_lines[nzchar(trimws(assay_lines))]
  
  if (length(assay_lines) < 2) {
    stop("Assay section was found, but no tabular data could be read.")
  }
  
  con <- textConnection(assay_lines)
  on.exit(close(con), add = TRUE)
  
  dat <- read.csv(
    con,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE,
    quote = "\"",
    comment.char = ""
  )
  
  return(dat)
}

# ---------------------------------------------------------------- #
# Extract SNP IDs + AlleleA probe sequences and write FASTA output
# ---------------------------------------------------------------- #
extract_manifest_probes_to_fasta <- function(input_file, output_fasta = NULL) {
  dat <- read_illumina_manifest_assay(input_file)
  
  # Candidate ID columns in order of preference
  id_candidates <- c("Name", "SNP Name", "SNP_Name", "Locus ID", "Locus_ID", "IlmnID")
  id_col <- id_candidates[id_candidates %in% colnames(dat)][1]
  
  if (is.na(id_col) || length(id_col) == 0) {
    stop(
      "Could not find a SNP identifier column. Tried: ",
      paste(id_candidates, collapse = ", ")
    )
  }
  
  if (!"AlleleA_ProbeSeq" %in% colnames(dat)) {
    stop("Column 'AlleleA_ProbeSeq' was not found in the assay table.")
  }
  
  snp_id <- trimws(dat[[id_col]])
  probe_seq <- toupper(gsub("\\s+", "", trimws(dat[["AlleleA_ProbeSeq"]])))
  
  keep <- !is.na(snp_id) & snp_id != "" & !is.na(probe_seq) & probe_seq != ""
  out <- data.frame(
    snp_id = snp_id[keep],
    probe_seq = probe_seq[keep],
    stringsAsFactors = FALSE
  )
  
  if (nrow(out) == 0) {
    stop("No valid SNP ID / probe sequence pairs were found.")
  }
  
  # Ensure FASTA headers are unique if duplicates exist
  if (anyDuplicated(out$snp_id)) {
    warning("Duplicate SNP IDs detected. Making FASTA headers unique with '__dup' suffix.")
    out$snp_id <- make.unique(out$snp_id, sep = "__dup")
  }
  
  if (is.null(output_fasta)) {
    base_name <- tools::file_path_sans_ext(basename(input_file))
    output_fasta <- file.path(dirname(input_file), paste0(base_name, "_AlleleA_ProbeSeq.fasta"))
  }
  
  con <- file(output_fasta, open = "w")
  on.exit(close(con), add = TRUE)
  
  for (i in seq_len(nrow(out))) {
    writeLines(paste0(">", out$snp_id[i]), con)
    writeLines(wrap_fasta(out$probe_seq[i], width = 60), con)
  }
  
  message("Input file:   ", input_file)
  message("ID column:    ", id_col)
  message("Sequences:    ", nrow(out))
  message("Output FASTA: ", output_fasta)
  
  invisible(out)
}

# ----------------------------- #
# Command-line interface
# ----------------------------- #
args <- commandArgs(trailingOnly = TRUE)

if (interactive()) {
  # Example runs for your two manifest files
  extract_manifest_probes_to_fasta("GGP65_manifest.csv")
  extract_manifest_probes_to_fasta("GGPPlus_manifest.csv")
} else {
  if (length(args) < 1 || length(args) > 2) {
    stop(
      "Usage:\n",
      "  Rscript extract_manifest_probes.R <input_manifest.csv> [output.fasta]\n\n",
      "Examples:\n",
      "  Rscript extract_manifest_probes.R GGP65_manifest.csv\n",
      "  Rscript extract_manifest_probes.R GGPPlus_manifest.csv GGPPlus_probes.fasta\n"
    )
  }
  
  input_file <- args[1]
  output_fasta <- if (length(args) == 2) args[2] else NULL
  
  extract_manifest_probes_to_fasta(input_file, output_fasta)
}
