#!/usr/bin/env Rscript

rm(list = ls())

input_file <- "/Volumes/Public/SNP_array_data/horse_data/ref/manifest_files/Affy_annotation_file.csv"
output_vcf <- "/Volumes/Public/SNP_array_data/horse_data/ref/manifest_files/Affy_annotation_file.vcf"
add_chr_prefix <- FALSE

clean_missing <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "---", "NA", "NaN", "NULL")] <- NA
  x
}

clean_base <- function(x) {
  toupper(clean_missing(x))
}

make_id <- function(cust_id, probe_set_id, affy_snp_id) {
  ifelse(!is.na(cust_id), cust_id,
         ifelse(!is.na(probe_set_id), probe_set_id,
                ifelse(!is.na(affy_snp_id), affy_snp_id, ".")))
}

escape_info <- function(x) {
  x <- clean_missing(x)
  x[is.na(x)] <- "."
  x <- gsub(";", ",", x, fixed = TRUE)
  x <- gsub("\\s+", "_", x)
  x
}

chrom_rank <- function(chr_vec) {
  chr_vec2 <- toupper(as.character(chr_vec))
  rank <- suppressWarnings(as.integer(chr_vec2))
  extra <- rep(999L, length(chr_vec2))
  extra[chr_vec2 == "X"] <- 1000L
  extra[chr_vec2 == "Y"] <- 1001L
  extra[chr_vec2 %in% c("MT", "M", "MITO", "MITOCHONDRIA")] <- 1002L
  ifelse(!is.na(rank), rank, extra)
}

# --------------------------------------------------
# Locate the actual table header:
# first non-empty line that does not start with '#'
# --------------------------------------------------
all_lines <- readLines(input_file, warn = FALSE, encoding = "UTF-8")
all_lines <- sub("^\ufeff", "", all_lines)

non_hash_idx <- which(
  nzchar(trimws(all_lines)) &
  !grepl("^#", trimws(all_lines))
)

if (length(non_hash_idx) == 0) {
  stop("Could not find any non-comment table lines in the file.")
}

header_idx <- non_hash_idx[1]
header_line <- all_lines[header_idx]

cat("Detected header line number:", header_idx, "\n")
cat("Header line content:\n", header_line, "\n\n", sep = "")

# detect separator
tab_count   <- lengths(regmatches(header_line, gregexpr("\t", header_line)))
comma_count <- lengths(regmatches(header_line, gregexpr(",", header_line)))
semi_count  <- lengths(regmatches(header_line, gregexpr(";", header_line)))

sep <- c("\t", ",", ";")[which.max(c(tab_count, comma_count, semi_count))]

cat("Detected separator:",
    if (sep == "\t") "TAB" else if (sep == ",") "COMMA" else "SEMICOLON",
    "\n")

ann <- read.table(
  file = input_file,
  sep = sep,
  header = TRUE,
  skip = header_idx - 1,
  quote = "\"",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fill = TRUE
)

colnames(ann)[1] <- sub("^\ufeff", "", colnames(ann)[1])

cat("Rows read:", nrow(ann), "\n")
cat("Columns read:", ncol(ann), "\n")
cat("First 20 column names:\n")
print(colnames(ann)[1:min(20, ncol(ann))])

required_cols <- c(
  "Probe Set ID",
  "Affy SNP ID",
  "Chromosome",
  "Physical Position",
  "Ref Allele",
  "Alt Allele",
  "cust_id"
)

missing_cols <- setdiff(required_cols, colnames(ann))
if (length(missing_cols) > 0) {
  cat("\nAll detected column names:\n")
  print(colnames(ann))
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- data.frame(
  CHROM        = clean_missing(ann[["Chromosome"]]),
  POS          = suppressWarnings(as.integer(clean_missing(ann[["Physical Position"]]))),
  REF          = clean_base(ann[["Ref Allele"]]),
  ALT          = clean_base(ann[["Alt Allele"]]),
  probe_set_id = clean_missing(ann[["Probe Set ID"]]),
  affy_snp_id  = clean_missing(ann[["Affy SNP ID"]]),
  cust_id      = clean_missing(ann[["cust_id"]]),
  stringsAsFactors = FALSE
)

df$ID <- make_id(df$cust_id, df$probe_set_id, df$affy_snp_id)

if (add_chr_prefix) {
  df$CHROM <- ifelse(is.na(df$CHROM), NA, paste0("chr", df$CHROM))
}

valid_base <- c("A", "C", "G", "T", "N")

keep <- !is.na(df$CHROM) &
        !is.na(df$POS) &
        !is.na(df$REF) &
        !is.na(df$ALT) &
        df$REF %in% valid_base &
        df$ALT %in% valid_base &
        df$REF != df$ALT

cat("Rows removed due to invalid/missing CHROM, POS, REF, or ALT:", sum(!keep), "\n")

df <- df[keep, , drop = FALSE]

if (nrow(df) == 0) {
  stop("No valid variants remained after filtering.")
}

df$chrom_order <- chrom_rank(df$CHROM)
df <- df[order(df$chrom_order, df$CHROM, df$POS, df$ID), , drop = FALSE]

df$INFO <- paste0(
  "PROBESET=", escape_info(df$probe_set_id),
  ";AFFYSNPID=", escape_info(df$affy_snp_id),
  ";CUSTID=", escape_info(df$cust_id)
)

vcf_df <- data.frame(
  "#CHROM" = df$CHROM,
  POS      = df$POS,
  ID       = df$ID,
  REF      = df$REF,
  ALT      = df$ALT,
  QUAL     = ".",
  FILTER   = "PASS",
  INFO     = df$INFO,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

con <- file(output_vcf, open = "wt")
writeLines("##fileformat=VCFv4.2", con)
writeLines("##source=Affy_annotation_to_VCF_R", con)
writeLines("##reference=EquCab3.0", con)
writeLines('##INFO=<ID=PROBESET,Number=1,Type=String,Description="Affymetrix Probe Set ID">', con)
writeLines('##INFO=<ID=AFFYSNPID,Number=1,Type=String,Description="Affymetrix SNP ID">', con)
writeLines('##INFO=<ID=CUSTID,Number=1,Type=String,Description="Custom marker ID from annotation file">', con)

write.table(
  vcf_df,
  file = con,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

close(con)

cat("VCF written to:", output_vcf, "\n")
cat("Final number of variants:", nrow(vcf_df), "\n")
