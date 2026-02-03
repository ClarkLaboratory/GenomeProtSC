library(optparse)
library(rtracklayer)
library(dplyr)        # For %>%
library(data.table)   # For fread()
library(tibble)
library(Matrix)       # For readMM()
library(BiocGenerics) # For a specific version of as.data.frame

set.seed(42)

option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL,
              help="FLAMES directory", metavar="character"),
  make_option(c("-m", "--mode"), type="character", default=NULL,
              help="Sample mode ('multi' or 'single')", metavar="character"),
  make_option(c("-i", "--single_sample_id"), type="character", default=NULL,
              help="Sample ID to use for single-sample mode (ignored if in multi-sample mode)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

flames_dir <- opt$directory
sample_mode <- opt$mode
single_sample_id <- opt$single_sample_id

# Rename oarfish output files with the sample ID when in single-sample mode
if (sample_mode == "single") {
  file.rename(from = file.path(flames_dir, "oarfish.barcodes.txt"), to = file.path(flames_dir, paste0(single_sample_id, ".barcodes.txt")))
  file.rename(from = file.path(flames_dir,    "oarfish.count.mtx"), to = file.path(flames_dir, paste0(single_sample_id,    ".count.mtx")))
  file.rename(from = file.path(flames_dir, "oarfish.features.txt"), to = file.path(flames_dir, paste0(single_sample_id, ".features.txt")))
}

# Rename the gene count file with the sample ID when in single-sample mode
if (sample_mode == "single") {
  file.rename(from = file.path(flames_dir, "gene_count.csv"), to = file.path(flames_dir, paste0(single_sample_id, "_gene_count.csv")))
}

# Function to make csv naming resource 
make_isoform_gene_symbol_dict <- function(FLAMES_gtf, reference_gtf, output_file) {
  # Import the first GTF file (transcripts GTF)
  gtf1 <- rtracklayer::readGFF(FLAMES_gtf)
  gtf1_df <- as.data.frame(gtf1) # Unnecessary? gtf1 should already be a data frame using the above method

  # Select relevant columns from the first GTF
  selected_columns1 <- gtf1_df[, c("transcript_id", "gene_id")]
  unique_selected_cols <- unique(selected_columns1)

  # Import the second GTF file (reference GTF with gene symbols)
  gtf2 <- import(reference_gtf)
  gtf2_df <- as.data.frame(gtf2)

  # Select relevant columns from the second GTF
  selected_columns2 <- gtf2_df[, c("gene_name", "gene_id")]
  unique_gene_symbol <- unique(selected_columns2)

  # Merge the two data frames on 'gene_id'
  combined_data <- merge(unique_selected_cols, unique_gene_symbol, by = "gene_id", all.x = TRUE)

  # If 'gene_name' is missing, replace it with 'gene_id'
  combined_data$gene_symbol <- ifelse(is.na(combined_data$gene_name), combined_data$gene_id, combined_data$gene_name)

  # Select relevant columns
  final_combined_data <- combined_data[, c("transcript_id", "gene_id", "gene_symbol")]

  # Write to a CSV file
  write.csv(final_combined_data, file = output_file, row.names = FALSE)

  return(final_combined_data)
}

# The FLAMES ref can be found in your selected output folder after running the Flames pipeline. 
FLAMES_gtf_file <- file.path(flames_dir, "isoform_annotated.gtf")
reference_gtf_file <- "/volstorage/refs/gencode.v47.annotation.gtf" # ensure file is unzipped
output_file <- file.path(flames_dir, "isoform_gene_dict.csv")

# Call the helper function defined in code block above to create a dictionary containing corresponding gene information for each isoform
# This may take a few minutes 
isoform_gene_dict <- make_isoform_gene_symbol_dict(FLAMES_gtf_file, reference_gtf_file, output_file)

convert_ENSGID_to_geneSymbol <- function(gene_count_matrix_path, id_symbol_df = isoform_gene_dict, output_file, return_df = FALSE) {
  # Load the reference dictionary we made earlier - select gene-level cols
  id_symbol_df <- as_tibble(id_symbol_df) %>%
    dplyr::select(gene_id, gene_symbol)

  # Load the data object with ENSGID row names
  gene_count_matrix <- fread(gene_count_matrix_path, header = TRUE)
  colnames(gene_count_matrix)[1] <- "gene_id"

  # Turn any missing gene count into 0
  for (i in names(gene_count_matrix)) {
    gene_count_matrix[is.na(get(i)), (i):=0]
  }

  # Replace ENSGIDs with gene symbols in original flames gene-level count matrix
  formatted_gene_count_matrix <- gene_count_matrix %>%
    merge(id_symbol_df, by.x = "gene_id", by.y = "gene_id") %>%     # Add gene symbol information
    distinct(gene_symbol, .keep_all = TRUE) %>%                     # Remove duplicates based on gene symbol
    dplyr::select(-gene_id) %>%                                     # Remove the ENSGID column
    column_to_rownames(var = "gene_symbol")                         # use the gene symbols we added as rownames

  # Write out the processed data frame
  fwrite(formatted_gene_count_matrix, output_file, row.names = TRUE)

  # Return the processed count matrix for further use if needed
  if (return_df) {
    return(formatted_gene_count_matrix)
  }
}

get_oarfish_df <- function(mtx_file, feature_file, barcode_file) {
  sparse_matrix <- t(readMM(mtx_file))
  rownames(sparse_matrix) <- read.table(feature_file)$V1
  colnames(sparse_matrix) <- read.table(barcode_file)$V1
  sparse_matrix <- as.data.frame(as.matrix(sparse_matrix))
  return(sparse_matrix)
}

annotate_oarfish_df_with_geneSymbols <- function(counts_df, id_symbol_df = isoform_gene_dict, output_file, return_df = FALSE) {
  # Add transcript_id as the first column
  counts_df$transcript_id <- rownames(counts_df)
  counts_df <- counts_df[, c(ncol(counts_df), 1:(ncol(counts_df) - 1))]

  # Merge with the resource table to add gene symbols
  df_genesymbol <- counts_df %>%
    left_join(id_symbol_df, by = "transcript_id")

  # Remove the gene_id column and reorder the columns
  df_genesymbol$gene_id <- NULL
  df_genesymbol <- df_genesymbol[, c(ncol(df_genesymbol), 1:(ncol(df_genesymbol) - 1))]

  # Update row names to include gene symbol instead of transcript_id
  rownames(df_genesymbol) <- paste0(df_genesymbol$transcript_id, "_", df_genesymbol$gene_symbol)
  df_genesymbol$transcript_id <- NULL
  df_genesymbol$gene_symbol <- NULL

  # Write the output to a CSV file
  fwrite(df_genesymbol, output_file, row.names = TRUE)
}

gene_count_files <- list.files(flames_dir, pattern = "_gene_count\\.csv$", full.names = TRUE)
for (i in seq_along(gene_count_files)) {
  # convert Gene_id to gene symbol for gene counts
  convert_ENSGID_to_geneSymbol(gene_count_matrix_path = gene_count_files[i], output_file = gene_count_files[i])
}

barcode_files <- list.files(flames_dir, pattern = "\\.barcodes\\.txt$", full.names = TRUE)
sample_ids <- gsub("\\.barcodes\\.txt$", "", basename(barcode_files))
sample_files <- lapply(sample_ids, function(sample_id) {
  list(
    barcode_file = file.path(flames_dir, paste0(sample_id, ".barcodes.txt")),
        mtx_file = file.path(flames_dir, paste0(sample_id,    ".count.mtx")),
    feature_file = file.path(flames_dir, paste0(sample_id, ".features.txt"))
  )
})
names(sample_files) <- sample_ids

for (sample_id in names(sample_files)) {
  files <- sample_files[[sample_id]]

  # get a data frame of the transcript counts matrix from Oarfish files
  oarfish_df <- get_oarfish_df(mtx_file = files$mtx_file, feature_file = files$feature_file, barcode_file = files$barcode_file)

  # annotate transcript IDs in the data frame with gene symbols (or gene IDs if there's no alternative)
  output_file <- file.path(flames_dir, paste0(sample_id, "_transcript_count.csv"))
  annotate_oarfish_df_with_geneSymbols(counts_df = oarfish_df, output_file = output_file)
}