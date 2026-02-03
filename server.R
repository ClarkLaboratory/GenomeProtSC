library(shiny)
library(shinyjs)

### USER INPUT VALIDATION FUNCTIONS: Return FALSE if user input doesn't pass the requirements, otherwise return TRUE ###

# FLAMES module
is_flames_input_valid <- function(input) {
  error_msg <- ""

  # Single-cell FASTQs: At least one file must be uploaded
  user_fastq_files <- input$user_fastq_files
  if (is.null(user_fastq_files)) {
    error_msg <- "Error: At least one single-cell FASTQ file must be uploaded."
    return(error_msg)
  }

  # Single-cell FASTQs: Must have a file extension of either '.fastq.gz' or '.fastq'
  num_files <- nrow(user_fastq_files)
  filenames <- user_fastq_files$name
  sample_ids <- gsub("\\.fastq(?:\\.gz)?$", "", filenames)
  for (i in 1:num_files) {
    invalid_sample_ids <- c("", filenames[i])
    if (is.element(sample_ids[i], invalid_sample_ids)) {
      error_msg <- "Error: Please ensure that all uploaded single-cell FASTQ files have a file extension of either '.fastq.gz' or '.fastq'."
      return(error_msg)
    }
  }

  # Expected cell count: Integer > 0
  expected_cell_count <- floor(input$expected_cell_count)
  if (!(is.finite(expected_cell_count) && (expected_cell_count > 0))) {
    error_msg <- paste0("Error: Expected cell count must be an integer larger than 0. Entered value: ", expected_cell_count)
    return(error_msg)
  }

  # Organism: "human" or "mouse"
  organism <- input$organism
  allowed_organisms <- c("human", "mouse")
  if (!is.element(organism, allowed_organisms)) {
    error_msg <- paste0("Error: Organism must be 'human' or 'mouse'. Entered value: ", organism)
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

# Database module
is_prot_db_input_valid <- function(input) {
  error_msg <- ""

  # Isoform annotation file: One file must be uploaded
  prot_db_gtf_file <- input$prot_db_gtf_file
  if (is.null(prot_db_gtf_file)) {
    error_msg <- "Error: An isoform annotation file must be uploaded."
    return(error_msg)
  }

  # Single-cell FASTQs: Must have a file extension of either '.fastq.gz' or '.fastq'
  filename <- prot_db_gtf_file$name
  if (!endsWith(filename, ".gtf")) {
    error_msg <- "Error: The isoform annotation file must have a file extension of '.gtf'."
    return(error_msg)
  }

  # Organism: "human" or "mouse"
  organism <- input$prot_db_organism
  allowed_organisms <- c("human", "mouse")
  if (!is.element(organism, allowed_organisms)) {
    error_msg <- paste0("Error: Organism must be 'human' or 'mouse'. Entered value: ", organism)
    return(error_msg)
  }

  # Minimum ORF length: Integer > 0
  min_orf_length <- floor(input$min_orf_length)
  if (!(is.finite(min_orf_length) && (min_orf_length > 0))) {
    error_msg <- paste0("Error: Minimum ORF length must be an integer larger than 0. Entered value: ", min_orf_length)
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

# The 'view statistics' part of the Seurat module
is_view_stats_input_valid <- function(input) {
  error_msg <- ""

  gene_count_files <- input$gene_count_files
  if (is.null(gene_count_files)) {
    error_msg <- "Error: At least one gene count file must be uploaded."
    return(error_msg)
  }

  num_files <- nrow(gene_count_files)
  filenames <- gene_count_files$name
  sample_ids <- gsub("_gene_count\\.csv$", "", filenames)
  for (i in 1:num_files) {
    invalid_sample_ids <- c("", filenames[i])
    if (is.element(sample_ids[i], invalid_sample_ids)) {
      error_msg <- "Error: Please ensure that all uploaded gene count files begin with a sample ID and end with '_gene_count.csv'."
      return(error_msg)
    }
  }

  transcript_count_files <- input$transcript_count_files
  if (is.null(transcript_count_files)) {
    error_msg <- "Error: At least one transcript count file must be uploaded."
    return(error_msg)
  }

  num_transcript_files <- nrow(transcript_count_files)
  if (num_transcript_files != num_files) {
    error_msg <- "Error: The number of transcript count files uploaded must be the same as the number of gene count files uploaded."
    return(error_msg)
  }

  filenames <- transcript_count_files$name
  transcript_sample_ids <- gsub("_transcript_count\\.csv$", "", filenames)
  for (i in 1:num_files) {
    invalid_sample_ids <- c("", filenames[i])
    if (is.element(transcript_sample_ids[i], invalid_sample_ids)) {
      error_msg <- "Error: Please ensure that all uploaded transcript count files begin with a sample ID and end with '_transcript_count.csv'."
      return(error_msg)
    }
  }

  if (!all(sort(sample_ids) == sort(transcript_sample_ids))) {
    error_msg <- "Error: The sample IDs used for gene count files and transcript count files must be identical."
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

# The 'view statistics after QC' part of the Seurat module
is_qc_view_stats_input_valid <- function(input) {
  error_msg <- ""

  gene_count_files <- input$gene_count_files
  if (is.null(gene_count_files)) {
    error_msg <- "Error: At least one gene count file must be uploaded."
    return(error_msg)
  }

  num_files <- nrow(gene_count_files)
  filenames <- gene_count_files$name
  sample_ids <- gsub("_gene_count\\.csv$", "", filenames)
  for (i in 1:num_files) {
    invalid_sample_ids <- c("", filenames[i])
    if (is.element(sample_ids[i], invalid_sample_ids)) {
      error_msg <- "Error: Please ensure that all uploaded gene count files begin with a sample ID and end with '_gene_count.csv'."
      return(error_msg)
    }
  }

  transcript_count_files <- input$transcript_count_files
  if (is.null(transcript_count_files)) {
    error_msg <- "Error: At least one transcript count file must be uploaded."
    return(error_msg)
  }

  num_transcript_files <- nrow(transcript_count_files)
  if (num_transcript_files != num_files) {
    error_msg <- "Error: The number of transcript count files uploaded must be the same as the number of gene count files uploaded."
    return(error_msg)
  }

  filenames <- transcript_count_files$name
  transcript_sample_ids <- gsub("_transcript_count\\.csv$", "", filenames)
  for (i in 1:num_files) {
    invalid_sample_ids <- c("", filenames[i])
    if (is.element(transcript_sample_ids[i], invalid_sample_ids)) {
      error_msg <- "Error: Please ensure that all uploaded transcript count files begin with a sample ID and end with '_transcript_count.csv'."
      return(error_msg)
    }
  }

  if (!all(sort(sample_ids) == sort(transcript_sample_ids))) {
    error_msg <- "Error: The sample IDs used for gene count files and transcript count files must be identical."
    return(error_msg)
  }

  # min_nCount_RNA: Integer >= 0
  min_nCount_RNA <- floor(input$min_nCount_RNA)
  if (!(is.finite(min_nCount_RNA) && (min_nCount_RNA >= 0))) {
    error_msg <- paste0("Error: Minimum number of molecules detected per cell must be a non-negative integer. Entered value: ", min_nCount_RNA)
    return(error_msg)
  }

  # max_nCount_RNA: Integer >= 0
  max_nCount_RNA <- floor(input$max_nCount_RNA)
  if (!(is.finite(max_nCount_RNA) && (max_nCount_RNA >= 0))) {
    error_msg <- paste0("Error: Maximum number of molecules detected per cell must be a non-negative integer. Entered value: ", max_nCount_RNA)
    return(error_msg)
  }

  # max_nCount_RNA: Must be larger than min_nCount_RNA
  if (max_nCount_RNA <= min_nCount_RNA) {
    error_msg <- paste0("Error: Maximum number of molecules detected per cell must be larger than the minimum number of molecules detected per cell. Entered value: ", max_nCount_RNA)
    return(error_msg)
  }

  # min_nFeature_RNA: Integer >= 0
  min_nFeature_RNA <- floor(input$min_nFeature_RNA)
  if (!(is.finite(min_nFeature_RNA) && (min_nFeature_RNA >= 0))) {
    error_msg <- paste0("Error: Minimum number of genes detected per cell must be a non-negative integer. Entered value: ", min_nFeature_RNA)
    return(error_msg)
  }

  # max_nFeature_RNA: Integer >= 0
  max_nFeature_RNA <- floor(input$max_nFeature_RNA)
  if (!(is.finite(max_nFeature_RNA) && (max_nFeature_RNA >= 0))) {
    error_msg <- paste0("Error: Maximum number of genes detected per cell must be a non-negative integer. Entered value: ", max_nFeature_RNA)
    return(error_msg)
  }

  # max_nFeature_RNA: Must be larger than min_nFeature_RNA
  if (max_nFeature_RNA <= min_nFeature_RNA) {
    error_msg <- paste0("Error: Maximum number of genes detected per cell must be larger than the minimum number of genes detected per cell. Entered value: ", max_nFeature_RNA)
    return(error_msg)
  }

  # max_percent_mt: Number >= 0
  max_percent_mt <- input$max_percent_mt
  if (!(is.finite(max_percent_mt) && (max_percent_mt >= 0))) {
    error_msg <- paste0("Error: Maximum percentage of molecules from mitochondrial genes per cell must be a non-negative number. Entered value: ", max_percent_mt)
    return(error_msg)
  }

  # npc_input: Integer > 1
  npc_input <- floor(input$npc_input)
  if (!(is.finite(npc_input) && (npc_input > 1))) {
    error_msg <- paste0("Error: Number of principal components must be an integer larger than 1. Entered value: ", npc_input)
    return(error_msg)
  }

  # k_input: Integer > 1
  k_input <- floor(input$k_input)
  if (!(is.finite(k_input) && (k_input > 1))) {
    error_msg <- paste0("Error: Number of neighbours must be an integer larger than 1. Entered value: ", k_input)
    return(error_msg)
  }

  # cluster_res_input: Number > 0
  cluster_res_input <- input$cluster_res_input
  if (!(is.finite(cluster_res_input) && (cluster_res_input > 0))) {
    error_msg <- paste0("Error: Cluster resolution must be a positive number. Entered value: ", cluster_res_input)
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

# The 'perform sample integration' part of the Seurat module
is_perform_integ_input_valid <- function(input) {
  error_msg <- ""

  # To perform sample integration, there must be multiple samples
  gene_count_files <- input$gene_count_files
  if (is.null(gene_count_files) || (nrow(gene_count_files) <= 1)) {
    error_msg <- "Error: At least two gene count files must be uploaded to perform sample integration."
    return(error_msg)
  }

  num_files <- nrow(gene_count_files)
  filenames <- gene_count_files$name
  sample_ids <- gsub("_gene_count\\.csv$", "", filenames)
  for (i in 1:num_files) {
    invalid_sample_ids <- c("", filenames[i])
    if (is.element(sample_ids[i], invalid_sample_ids)) {
      error_msg <- "Error: Please ensure that all uploaded gene count files begin with a sample ID and end with '_gene_count.csv'."
      return(error_msg)
    }
  }

  transcript_count_files <- input$transcript_count_files
  if (is.null(transcript_count_files)) {
    error_msg <- "Error: At least one transcript count file must be uploaded."
    return(error_msg)
  }

  num_transcript_files <- nrow(transcript_count_files)
  if (num_transcript_files != num_files) {
    error_msg <- "Error: The number of transcript count files uploaded must be the same as the number of gene count files uploaded."
    return(error_msg)
  }

  filenames <- transcript_count_files$name
  transcript_sample_ids <- gsub("_transcript_count\\.csv$", "", filenames)
  for (i in 1:num_files) {
    invalid_sample_ids <- c("", filenames[i])
    if (is.element(transcript_sample_ids[i], invalid_sample_ids)) {
      error_msg <- "Error: Please ensure that all uploaded transcript count files begin with a sample ID and end with '_transcript_count.csv'."
      return(error_msg)
    }
  }

  if (!all(sort(sample_ids) == sort(transcript_sample_ids))) {
    error_msg <- "Error: The sample IDs used for gene count files and transcript count files must be identical."
    return(error_msg)
  }

  # min_nCount_RNA: Integer >= 0
  min_nCount_RNA <- floor(input$min_nCount_RNA)
  if (!(is.finite(min_nCount_RNA) && (min_nCount_RNA >= 0))) {
    error_msg <- paste0("Error: Minimum number of molecules detected per cell must be a non-negative integer. Entered value: ", min_nCount_RNA)
    return(error_msg)
  }

  # max_nCount_RNA: Integer >= 0
  max_nCount_RNA <- floor(input$max_nCount_RNA)
  if (!(is.finite(max_nCount_RNA) && (max_nCount_RNA >= 0))) {
    error_msg <- paste0("Error: Maximum number of molecules detected per cell must be a non-negative integer. Entered value: ", max_nCount_RNA)
    return(error_msg)
  }

  # max_nCount_RNA: Must be larger than min_nCount_RNA
  if (max_nCount_RNA <= min_nCount_RNA) {
    error_msg <- paste0("Error: Maximum number of molecules detected per cell must be larger than the minimum number of molecules detected per cell. Entered value: ", max_nCount_RNA)
    return(error_msg)
  }

  # min_nFeature_RNA: Integer >= 0
  min_nFeature_RNA <- floor(input$min_nFeature_RNA)
  if (!(is.finite(min_nFeature_RNA) && (min_nFeature_RNA >= 0))) {
    error_msg <- paste0("Error: Minimum number of genes detected per cell must be a non-negative integer. Entered value: ", min_nFeature_RNA)
    return(error_msg)
  }

  # max_nFeature_RNA: Integer >= 0
  max_nFeature_RNA <- floor(input$max_nFeature_RNA)
  if (!(is.finite(max_nFeature_RNA) && (max_nFeature_RNA >= 0))) {
    error_msg <- paste0("Error: Maximum number of genes detected per cell must be a non-negative integer. Entered value: ", max_nFeature_RNA)
    return(error_msg)
  }

  # max_nFeature_RNA: Must be larger than min_nFeature_RNA
  if (max_nFeature_RNA <= min_nFeature_RNA) {
    error_msg <- paste0("Error: Maximum number of genes detected per cell must be larger than the minimum number of genes detected per cell. Entered value: ", max_nFeature_RNA)
    return(error_msg)
  }

  # max_percent_mt: Number >= 0
  max_percent_mt <- input$max_percent_mt
  if (!(is.finite(max_percent_mt) && (max_percent_mt >= 0))) {
    error_msg <- paste0("Error: Maximum percentage of molecules from mitochondrial genes per cell must be a non-negative number. Entered value: ", max_percent_mt)
    return(error_msg)
  }

  # npc_input: Integer > 1
  npc_input <- floor(input$npc_input)
  if (!(is.finite(npc_input) && (npc_input > 1))) {
    error_msg <- paste0("Error: Number of principal components must be an integer larger than 1. Entered value: ", npc_input)
    return(error_msg)
  }

  # k_input: Integer > 1
  k_input <- floor(input$k_input)
  if (!(is.finite(k_input) && (k_input > 1))) {
    error_msg <- paste0("Error: Number of neighbours must be an integer larger than 1. Entered value: ", k_input)
    return(error_msg)
  }

  # cluster_res_input: Number > 0
  cluster_res_input <- input$cluster_res_input
  if (!(is.finite(cluster_res_input) && (cluster_res_input > 0))) {
    error_msg <- paste0("Error: Cluster resolution must be a positive number. Entered value: ", cluster_res_input)
    return(error_msg)
  }

  # npcs: Integer > 1
  npcs <- floor(input$npcs)
  if (!(is.finite(npcs) && (npcs > 1))) {
    error_msg <- paste0("Error: npcs must be an integer larger than 1. Entered value: ", npcs)
    return(error_msg)
  }

  # theta: Number > 0
  theta <- input$theta
  if (!(is.finite(theta) && (theta > 0))) {
    error_msg <- paste0("Error: theta must be a positive number. Entered value: ", theta)
    return(error_msg)
  }

  # lambda: Number > 0
  lambda <- input$lambda
  if (!(is.finite(lambda) && (lambda > 0))) {
    error_msg <- paste0("Error: lambda must be a positive number. Entered value: ", lambda)
    return(error_msg)
  }

  # sigma: Number > 0
  sigma <- input$sigma
  if (!(is.finite(sigma) && (sigma > 0))) {
    error_msg <- paste0("Error: sigma must be a positive number. Entered value: ", sigma)
    return(error_msg)
  }

  # nclust: Either -1 or a positive integer
  nclust <- floor(input$nclust)
  if (!(is.finite(nclust) && ((nclust >= 1) || (nclust == -1)))) {
    error_msg <- paste0("Error: nclust must be either -1 (default value) or a positive integer. Entered value: ", nclust)
    return(error_msg)
  }

  # tau: Integer > 1
  tau <- floor(input$tau)
  if (!(is.finite(tau) && (tau >= 0))) {
    error_msg <- paste0("Error: tau must be a non-negative integer. Entered value: ", tau)
    return(error_msg)
  }

  # block.size: Number between 0 and 1
  block.size <- input$block.size
  if (!(is.finite(block.size) && (block.size >= 0) && (block.size <= 1))) {
    error_msg <- paste0("Error: block.size must be a number between 0 and 1. Entered value: ", block.size)
    return(error_msg)
  }

  # max.iter.harmony: Positive integer
  max.iter.harmony <- floor(input$max.iter.harmony)
  if (!(is.finite(max.iter.harmony) && (max.iter.harmony > 0))) {
    error_msg <- paste0("Error: max.iter.harmony must be a positive integer. Entered value: ", max.iter.harmony)
    return(error_msg)
  }

  # max.iter.cluster: Positive integer
  max.iter.cluster <- floor(input$max.iter.cluster)
  if (!(is.finite(max.iter.cluster) && (max.iter.cluster > 0))) {
    error_msg <- paste0("Error: max.iter.cluster must be a positive integer. Entered value: ", max.iter.cluster)
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

# scRNA-seq + proteomics integration module
is_integ_input_valid <- function(input) {
  error_msg <- ""

  # Proteomics results: A file must be uploaded
  user_proteomics_file <- input$user_proteomics_file
  if (is.null(user_proteomics_file)) {
    error_msg <- "Error: Proteomics results must be uploaded."
    return(error_msg)
  }

  # 'proteome_database.fasta': A file must be uploaded
  user_fasta_file <- input$user_fasta_file
  if (is.null(user_fasta_file)) {
    error_msg <- "Error: 'proteome_database.fasta' must be uploaded."
    return(error_msg)
  }

  # 'proteome_database_metadata.txt': A file must be uploaded
  user_metadata_file <- input$user_metadata_file
  if (is.null(user_metadata_file)) {
    error_msg <- "Error: 'proteome_database_metadata.txt' must be uploaded."
    return(error_msg)
  }

  # 'proteome_database_transcripts.gtf': A file must be uploaded
  user_post_gtf_file <- input$user_post_gtf_file
  if (is.null(user_post_gtf_file)) {
    error_msg <- "Error: 'proteome_database_transcripts.gtf' must be uploaded."
    return(error_msg)
  }

  # Seurat object RDS file: A file must be uploaded
  user_rds_file <- input$user_rds_file
  if (is.null(user_rds_file)) {
    error_msg <- "Error: A Seurat object RDS file must be uploaded."
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

is_integ_marker_input_valid <- function(input) {
  error_msg <- ""

  # 'combined_annotations.gtf': A file must be uploaded
  user_gtf_file <- input$user_gtf_file
  if (is.null(user_gtf_file)) {
    error_msg <- "Error: 'combined_annotations.gtf' must be uploaded."
    return(error_msg)
  }

  # 'marker_genes.csv': A file must be uploaded
  user_marker_file <- input$user_marker_file
  if (is.null(user_marker_file)) {
    error_msg <- "Error: 'marker_genes.csv' from module 1b must be uploaded."
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

# Visualization module
is_vis_isovis_input_valid <- function(input) {
  error_msg <- ""

  transcript_info_file <- input$transcript_info_file_isovis
  seurat_rds_file <- input$seurat_rds_file_isovis

  if (is.null(transcript_info_file) && is.null(seurat_rds_file)) {
    error_msg <- "Error: Please upload a Seurat RDS file and optionally 'transcript_expression_info.csv'."
    return(error_msg)
  }

  if (!is.null(transcript_info_file) && is.null(seurat_rds_file)) {
    error_msg <- "Error: A Seurat object needs to be uploaded if 'transcript_expression_info.csv' is uploaded."
    return(error_msg)
  }

  # All checks pass, so the user input is valid
  return(error_msg)
}

# Function for exporting VSN peptide counts
normalize_peptide_counts <- function(input, session_id) {
  error_msg <- ""

  # Peptide counts file: One file must be uploaded
  peptide_counts_file <- input$peptide_counts_file
  if (is.null(peptide_counts_file)) {
    error_msg <- "Error: A peptide counts file must be uploaded."
    return(error_msg)
  }

  # The peptide counts file must be a CSV or tab-separated text file
  filename <- peptide_counts_file$name
  if (!endsWith(filename, ".csv") && !endsWith(filename, ".tsv") && !endsWith(filename, ".txt")) {
    error_msg <- "Error: The peptide counts file must have a file extension of '.csv', '.tsv', or '.txt' (tab-separated text)."
    return(error_msg)
  }

  countsp <- fread(peptide_counts_file$datapath)

  # store peptide sequences in the 'Peptide' column if it doesn't already exist
  if (!("Peptide" %in% colnames(countsp))) {
    if ("Stripped.Sequence" %in% colnames(countsp)) {
      countsp$Peptide <- countsp$Stripped.Sequence
      countsp$Stripped.Sequence <- NULL
    } else if ("peptide" %in% colnames(countsp)) {
      countsp$Peptide <- countsp$peptide
      countsp$peptide <- NULL
    } else if ("peptide_sequence" %in% colnames(countsp)) {
      countsp$Peptide <- countsp$peptide_sequence
      countsp$peptide_sequence <- NULL
    }
  }

  if (!("Peptide" %in% colnames(countsp))) {
    error_msg <- "Error: No column named 'Peptide', 'peptide', 'Stripped.Sequence' or 'peptide_sequence' found. Ensure this column contains the amino acid sequences of peptides and is present in the uploaded file."
    return(error_msg)
  }

  # filter for sample names
  sample_names <- c()
  forbidden_column_names <- c("1/k0", "1/k0", "all mapped proteins", "all mapped genes", "average.missed.tryptic.cleavages", "average.peptide.charge", "average.peptide.length",
                              "best.fr.mz", "best.fr.mz.delta", "channel", "channel.evidence", "channel.l", "channel.q.value", "corr", "cscore", "decoy", "decoy.cscore", "decoy.evidence",
                              "delta", "description", "empirical.quality", "evidence", "exclude.from.quant", "excludefromassay", "file.name", "filename", "first.protein.description",
                              "fragment.charge", "fragment.correlations", "fragment.info", "fragment.loss.type", "fragment.quant.corrected", "fragment.quant.raw", "fragment.series.number",
                              "fragment.sum", "fragment.type", "fragmentcharge", "fragmentlosstype", "fragmentseriesnumber", "fragmenttype", "fullunimodpeptidename", "fwhm", "fwhm.rt",
                              "fwhm.scans", "gene", "gene.names", "genes", "genes.maxlfq", "genes.maxlfq.quality", "genes.maxlfq.unique", "genes.maxlfq.unique.quality", "genes.normalised",
                              "genes.quantity", "genes.topn", "gg.q.value", "global.peptidoform.q.value", "global.pg.q.value", "global.q.value", "intensities", "ion.mobility", "irt",
                              "label.ratio", "lib.peptidoform.q.value", "lib.pg.q.value", "lib.ptm.site.confidence", "lib.q.value", "libraryintensity", "m/z", "mass.evidence",
                              "median.mass.acc.ms1", "median.mass.acc.ms1.corrected", "median.mass.acc.ms2", "median.mass.acc.ms2.corrected", "median.rt.prediction.acc", "modification",
                              "modified.sequence", "modifiedpeptide", "ms.level", "ms1.apex.area", "ms1.apex.mz.delta", "ms1.area", "ms1.normalised", "ms1.profile.corr", "ms1.signal",
                              "ms1.total.signal.after", "ms1.total.signal.before", "ms2.scan", "ms2.scan", "ms2.signal", "n.proteotypic.sequences", "n.sequences", "normalisation.factor",
                              "normalisation.instability", "normalisation.noise", "peptide", "peptidegrouplabel", "peptidesequence", "peptidoform.q.value", "pg.maxlfq", "pg.maxlfq.quality",
                              "pg.normalised", "pg.pep", "pg.q.value", "pg.quantity", "pg.topn", "pgqvalue", "precursor.charge", "precursor.id", "precursor.lib.index", "precursor.mz",
                              "precursor.normalised", "precursor.quantity", "precursorcharge", "precursormz", "precursors.identified", "predicted.iim", "predicted.im", "predicted.irt",
                              "predicted.rt", "probability", "product.mz", "productmz", "protein", "protein.group", "protein.id", "protein.ids", "protein.index.in.group", "protein.name",
                              "protein.names", "protein.q.value", "protein.sites", "proteingroup", "proteinname", "proteins.identified", "proteotypic", "ptm.site.confidence", "q.value",
                              "quantity.quality", "qvalue", "relative.intensity", "residue", "retention.times", "rt.start", "rt.stop", "run", "run.index", "sequence", "site",
                              "site.occupancy.probabilities", "theoretical.mz", "total.quantity", "translated.q.value", "uniprotid", "window.high", "window.low")
  column_names <- colnames(countsp)
  for (i in seq_along(column_names)) {
    column_name <- tolower(column_names[i])
    if (!(column_name %in% forbidden_column_names)) {
      sample_names <- c(sample_names, column_names[i])
    }
  }

  # order sample names
  sample_names <- sample_names[order(match(sample_names, colnames(countsp)))]

  # rarely, a peptide is in the data twice, so take max count value and get unique peptide IDs
  countsp <- countsp %>% 
    dplyr::select(tidyselect::all_of(c("Peptide", sample_names))) %>% 
    dplyr::mutate(sum = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>% 
    dplyr::group_by(Peptide) %>% 
    dplyr::slice_max(sum) %>% dplyr::ungroup() %>% dplyr::select(-sum)

  # convert to matrix
  countsp_matrix <- as.matrix(countsp[,-1])

  # run VSN for normalisation of peptide intensities
  # TODO: justvsn() errors out when there are less rows than data points per stratum (allow users to modify the justvsn parameters?)
  min_data_points_per_stratum <- 42
  if (nrow(countsp_matrix) < min_data_points_per_stratum) {
    min_data_points_per_stratum <- nrow(countsp_matrix)
  }

  if (min_data_points_per_stratum == 0) {
    error_msg <- "Error: No data rows found from the peptide counts. Please ensure there is at least one data row in the uploaded file."
    return(error_msg)
  }

  vsnp <- as.data.frame(justvsn(countsp_matrix, minDataPointsPerStratum = min_data_points_per_stratum))

  # prepare and export the normalized counts
  peptide <- countsp$Peptide
  vsnp <- cbind(peptide, vsnp)

  write.csv(vsnp, file.path(session_id, "vsn_peptide_counts.csv"), row.names = FALSE)

  # no errors occurred; exit
  return(error_msg)
}

### Helper functions and variables for hiding and showing different parts of the server ###

show_elements <- function(element_ids) {
  for (i in seq_along(element_ids)) {
    showElement(element_ids[i])
  }
}

hide_elements <- function(element_ids) {
  for (i in seq_along(element_ids)) {
    hide(element_ids[i])
  }
}

seurat_module_outputs <- c("show_statistics_output_part", "hide_statistics_output_part", "statistics_heading",
                           "scatter_plot_heading", "scatter_plot", "vln_plot_heading", "vln_plot",
                           "elbow_plot_heading", "elbow_plot", "umap_plot_heading", "umap_plot",
                           "marker_gene_plot_heading", "marker_gene_plot")

seurat_module_qc_params_parts <- c("show_qc_part", "hide_qc_part", "qc_heading", "show_adv_qc_params_button", "hide_adv_qc_params_button",
                                   "min_nCount_RNA", "max_nCount_RNA", "min_nFeature_RNA", "max_nFeature_RNA", "max_percent_mt",
                                   "npc_input", "k_input", "cluster_res_input",
                                   "view_qc_stats_button")

seurat_module_harmony_integ_params_parts <- c("show_integ_part", "hide_integ_part", "harmony_integ_heading", "show_adv_harmony_integ_params_button", "hide_adv_harmony_integ_params_button",
                                              "tau", "theta",
                                              "npcs", "lambda", "sigma", "nclust", "block.size", "max.iter.harmony", "max.iter.cluster", "epsilon.cluster", "epsilon.harmony",
                                              "perform_integration_button",
                                              "harmony_integ_footer_1", "harmony_integ_footer_2")

seurat_module_rename_clusters_parts <- c("show_rename_clusters_part", "hide_rename_clusters_part", "rename_clusters_heading", "rename_clusters_info",
                                         "cluster_list", "rename_clusters_button")

seurat_module_find_marker_genes_parts <- c("show_marker_genes_part", "hide_marker_genes_part", "find_marker_genes_heading", "find_marker_genes_info",
                                           "ident.1", "ident.2", "p_adj.thresh", "logfc.threshold", "min.pct", "min.diff.pct", "only.pos", "find_marker_genes_button")

seurat_module_find_all_marker_genes_parts <- c("show_all_marker_genes_part", "hide_all_marker_genes_part", "find_all_marker_genes_heading",
                                               "p_adj.thresh_all", "logfc.threshold_all", "min.pct_all", "min.diff.pct_all", "only.pos_all", "find_all_marker_genes_button")

do_show_qc_part_action <- function() {
  hide("show_qc_part")
  show_elements(seurat_module_qc_params_parts[c(2:4, 6:10, 14)])
  hide_elements(seurat_module_qc_params_parts[c(5, 11:13)])
}

do_hide_qc_part_action <- function() {
  showElement("show_qc_part")
  hide_elements(seurat_module_qc_params_parts[c(2:length(seurat_module_qc_params_parts))])
}

do_show_integ_part_action <- function() {
  hide("show_integ_part")
  show_elements(seurat_module_harmony_integ_params_parts[c(2:4, 6:7, 17:19)])
  hide_elements(seurat_module_harmony_integ_params_parts[c(1, 5, 8:16)])
}

do_hide_integ_part_action <- function() {
  showElement("show_integ_part")
  hide_elements(seurat_module_harmony_integ_params_parts[c(2:length(seurat_module_harmony_integ_params_parts))])
}

do_show_rename_clusters_part_action <- function() {
  hide("show_rename_clusters_part")
  show_elements(seurat_module_rename_clusters_parts[2:length(seurat_module_rename_clusters_parts)])
}

do_hide_rename_clusters_part_action <- function() {
  showElement("show_rename_clusters_part")
  hide_elements(seurat_module_rename_clusters_parts[2:length(seurat_module_rename_clusters_parts)])
}

do_show_marker_genes_part_action <- function() {
  hide("show_marker_genes_part")
  show_elements(seurat_module_find_marker_genes_parts[2:length(seurat_module_find_marker_genes_parts)])
}

do_hide_marker_genes_part_action <- function() {
  showElement("show_marker_genes_part")
  hide_elements(seurat_module_find_marker_genes_parts[2:length(seurat_module_find_marker_genes_parts)])
}

do_show_all_marker_genes_part_action <- function() {
  hide("show_all_marker_genes_part")
  show_elements(seurat_module_find_all_marker_genes_parts[2:length(seurat_module_find_all_marker_genes_parts)])
}

do_hide_all_marker_genes_part_action <- function() {
  showElement("show_all_marker_genes_part")
  hide_elements(seurat_module_find_all_marker_genes_parts[2:length(seurat_module_find_all_marker_genes_parts)])
}

do_show_adv_qc_params_action <- function() {
  hide("show_adv_qc_params_button")
  show_elements(seurat_module_qc_params_parts[c(5, 11:13)])
}

do_hide_adv_qc_params_action <- function() {
  showElement("show_adv_qc_params_button")
  hide_elements(seurat_module_qc_params_parts[c(5, 11:13)])
}

do_show_adv_harmony_integ_params_action <- function() {
  hide("show_adv_harmony_integ_params_button")
  show_elements(seurat_module_harmony_integ_params_parts[c(5, 8:16)])
}

do_hide_adv_harmony_integ_params_action <- function() {
  showElement("show_adv_harmony_integ_params_button")
  hide_elements(seurat_module_harmony_integ_params_parts[c(5, 8:16)])
}

do_show_statistics_output_part_action <- function() {
  hide("show_statistics_output_part")
  showElement("seurat_statistics_msg_container")
  show_elements(seurat_module_outputs[c(2:9)])
}

do_hide_statistics_output_part_action <- function() {
  showElement("show_statistics_output_part")
  hide("seurat_statistics_msg_container")
  hide_elements(seurat_module_outputs[c(2:9)])
}

hide_post_qc_parts_seurat_module <- function() {
  hide_elements(seurat_module_rename_clusters_parts)
  hide_elements(seurat_module_find_marker_genes_parts)
  hide_elements(seurat_module_find_all_marker_genes_parts)
}

show_default_seurat_module <- function() {
  # Hide plots, the rename clusters menu, the find marker genes menu, and the find all marker genes menu
  hide_elements(seurat_module_outputs)
  hide_post_qc_parts_seurat_module()

  # Hide the quality control menu behind a button
  showElement(seurat_module_qc_params_parts[1])
  hide_elements(seurat_module_qc_params_parts[2:length(seurat_module_qc_params_parts)])

  # Hide the Harmony sample integration menu behind a button
  showElement(seurat_module_harmony_integ_params_parts[1])
  hide_elements(seurat_module_harmony_integ_params_parts[2:length(seurat_module_harmony_integ_params_parts)])
}

show_basic_outputs_seurat_module <- function() {
  showElement("seurat_statistics_msg_container")
  show_elements(seurat_module_outputs[3:7])
}

show_more_outputs_seurat_module <- function() {
  showElement("seurat_statistics_msg_container")
  do_hide_rename_clusters_part_action()
  show_elements(seurat_module_outputs[2:11])
  show_elements(c("show_marker_genes_part", "show_all_marker_genes_part"))
}

show_marker_gene_outputs_seurat_module <- function() {
  show_elements(c("marker_gene_plot_heading", "marker_gene_plot"))
}

hide_marker_gene_outputs_seurat_module <- function() {
  hide_elements(c("marker_gene_plot_heading", "marker_gene_plot"))
}

show_rename_clusters_part_and_umap_outputs_seurat_module <- function() {
  show_elements(c("show_rename_clusters_part", "umap_plot_heading", "umap_plot"))
}

### INTERNAL SERVER FUNCTIONS ###

flames_server <- function(input, is_multi_sample_mode, single_sample_id, session) {
  # store session ID
  session_id <- session$token

  # set output dir
  outdir_flames <- file.path(getwd(), session_id, "single_cell_output")

  # create output dir
  dir.create(outdir_flames, recursive = TRUE)
  dir.create(file.path(outdir_flames, "flames"), recursive = TRUE) # The directory the config file will be stored in also needs to be created

  old_wd <- getwd()
  setwd(outdir_flames)

  # set genome
  if (input$organism == "human") {
    genome_file <- "/volstorage/refs/human.fasta"
  } else if (input$organism == "mouse") {
    genome_file <- "/volstorage/refs/mouse.fasta"
  }

  # process uploaded fastq files
  user_fastq_files <- input$user_fastq_files
  renamed_files <- file.path(dirname(user_fastq_files$datapath), user_fastq_files$name)

  # rename to original file names
  mapply(file.rename, from = user_fastq_files$datapath, to = renamed_files)
  message("Renamed FASTQ files: ", paste(renamed_files, collapse = ", "))

  # get location of uploaded fastqs
  fastq_dir <- dirname(input$user_fastq_files$datapath[1])

  bambu_ndr <- input$bambu_ndr
  if (!(is.numeric(bambu_ndr) && is.finite(bambu_ndr) && ((bambu_ndr >= 0.0) && (bambu_ndr <= 1.0)))) {
    bambu_ndr <- -1
  }

  # Run FLAMES on the provided FASTQs
  command_run_flames <- paste0(
    "Rscript ../../bin/database_module/flames-run.R",
    " --genome ", shQuote(genome_file),
    " --annotation ", "/volstorage/refs/gencode.v47.annotation.gtf",
    " --cell_count ", shQuote(floor(input$expected_cell_count)),
    " --ndr ", shQuote(bambu_ndr),
    " --input ", shQuote(fastq_dir),
    " --output ", shQuote(file.path(outdir_flames, "flames"))
  )

  print(command_run_flames)
  system(command_run_flames)

  # export cell counts as files
  mode <- ifelse(is_multi_sample_mode, "multi", "single")
  command_flames_counts <- paste0(
    "Rscript ../../bin/database_module/flames-export-counts.R",
    " --directory ", shQuote(file.path(outdir_flames, "flames")),
    " --mode ", shQuote(mode),
    " --single_sample_id ", shQuote(single_sample_id)
  )

  print(command_flames_counts)
  system(command_flames_counts)

  # Export the files
  top_level_dir <- getwd()
  output_files_dir <- file.path(outdir_flames, "flames")
  setwd(output_files_dir)

  gene_count_files <- list.files(pattern = "_gene_count\\.csv$")
  transcript_count_files <- list.files(pattern = "_transcript_count\\.csv$")

  files_to_zip_db <- c("isoform_annotated.gtf", gene_count_files, transcript_count_files)
  zipfile_path_db <- file.path(old_wd, session_id, "database_results.zip")
  zip(zipfile = zipfile_path_db, files = files_to_zip_db)

  setwd(top_level_dir)
  setwd(old_wd)
}

database_server <- function(input, session) {
  # store session ID
  session_id <- session$token

  # set output dir
  outdir_db <- file.path(session_id, "database_output")

  # create output dir
  dir.create(outdir_db, recursive = TRUE)

  # Get the filename of the Bambu transcript annotations
  db_gtf_file <- input$prot_db_gtf_file$datapath

  # construct the command
  command_generate_proteome <- paste0(
    "Rscript bin/database_module/generate_proteome.R",
    " --gtf ", shQuote(db_gtf_file),
    " --reference ", "/volstorage/refs/gencode.v47.annotation.gtf",
    " --organism ", shQuote(input$organism),
    " --length ", shQuote(floor(input$min_orf_length)),
    " --uorfs ", shQuote(input$user_find_utr_5_orfs),
    " --dorfs ", shQuote(input$user_find_utr_3_orfs),
    " --savepath ", shQuote(outdir_db)
  )

  print(command_generate_proteome)
  system(command_generate_proteome)

  print("Generated ORFs")

  # set reference protein database per organism 
  if (input$organism == "human") {
    ref_proteome <- "data/openprot_uniprotDb_hs.txt"
  } else if (input$organism == "mouse") {
    ref_proteome <- "data/openprot_uniprotDb_mm.txt"
  }

  # run python script to create proteome fasta
  command_annotate_proteome <- paste0(
    ". /miniconda3/etc/profile.d/conda.sh; ",
    "conda activate; ",
    "python bin/database_module/annotate_proteome.py ",
    "/volstorage/refs/gencode.v47.annotation.gtf ", 
    shQuote(ref_proteome), " ",
    shQuote(file.path(outdir_db, "ORFome_aa.txt")), " ",
    shQuote(file.path(outdir_db, "proteome_database_transcripts.gtf")), " ",
    shQuote(outdir_db), " ",
    "all ",
    shQuote(floor(input$min_orf_length)), " ",
    "None"
  )

  print(command_annotate_proteome)
  system(command_annotate_proteome)

  # Export the files
  top_level_dir <- getwd()
  setwd(outdir_db)

  files_to_zip_db <- c("proteome_database.fasta", "proteome_database_metadata.txt", "proteome_database_transcripts.gtf")
  zipfile_path_db <- file.path(top_level_dir, session_id, "protein_database_results.zip")
  zip(zipfile = zipfile_path_db, files = files_to_zip_db)

  setwd(top_level_dir)
}

integration_server <- function(input, output, session) {
  # req(input$user_proteomics_file, input$user_post_gtf_file, input$user_fasta_file, input$user_metadata_file)  # GTF is required

  # store session ID
  session_id <- session$token

  # set output dir
  outdir_integ <- file.path(session_id, "integ_output")

  # create output dir
  dir.create(outdir_integ, recursive = TRUE)

  # run Rscript
  command_map_peptides <- paste0(
    "Rscript bin/integration_module/map_peptides_generate_outputs.R",
    " --proteomics ", shQuote(input$user_proteomics_file$datapath),
    " --fasta ", shQuote(input$user_fasta_file$datapath),
    " --metadata ", shQuote(input$user_metadata_file$datapath),
    " --gtf ", shQuote(input$user_post_gtf_file$datapath),
    " --rds ", shQuote(input$user_rds_file$datapath),
    " --savepath ", shQuote(outdir_integ)
  )

  print(command_map_peptides)
  system(command_map_peptides)

  # get the top level dir
  top_level_dir <- getwd()
  dir.create(file.path(outdir_integ, "report_images"))

  # create report
  rmarkdown::render(input = file.path(top_level_dir, "bin", "integration_module", "integration_summary_report.Rmd"),
                    output_file = file.path(top_level_dir, outdir_integ, "summary_report.html"),
                    output_format = "html_document",
                    params = list(
                      directory = file.path(top_level_dir, outdir_integ),
                      file = "peptide_info.csv"
                    ))

  # zip all the result files
  if (file.exists(file.path(outdir_integ, "peptide_info.csv")) && file.exists(file.path(outdir_integ, "summary_report.html"))) {
    files_to_zip_int <- c("summary_report.html", "peptide_info.csv", "report_images/",
                          "combined_annotations.gtf", "transcripts_and_ORFs_for_isovis.gtf",
                          "peptides.bed12", "ORFs.bed12", "transcripts.bed12", "transcript_expression_info.csv")
    zipfile_path_int <- file.path("..", "integration_results.zip")
    setwd(outdir_integ)
    zip(zipfile = zipfile_path_int, files = files_to_zip_int)
    setwd(top_level_dir)
  }
}

integration_marker_server <- function(input, session) {
  # store session ID
  session_id <- session$token

  marker_genes_csv_filename <- input$user_marker_file$datapath
  marker_genes_csv <- read.csv(marker_genes_csv_filename)

  gtf_annotation_filename <- input$user_gtf_file$datapath
  new_gtf_annotation_filename <- file.path(session_id, "combined_annotations_with_marker_genes.gtf")

  gtf_annotation <- file(gtf_annotation_filename, 'r')
  new_gtf_annotation <- file(new_gtf_annotation_filename, 'w')

  while (TRUE) {
    line <- readLines(gtf_annotation, n = 1)
    if (length(line) == 0) {
      break
    }

    cols <- base::strsplit(line, '\t')[[1]]
    if ((length(cols) != 9) || (cols[2] != "Bambu") || (cols[3] != "exon"))
    {
      writeLines(line, new_gtf_annotation, useBytes = TRUE)
      next
    }

    start_index <- gregexpr('gene_name "', cols[9], fixed = TRUE)[[1]][1]
    if (start_index == -1)
    {
      writeLines(line, new_gtf_annotation, useBytes = TRUE)
      next
    }

    start_index <- start_index + 11
    gene_name <- substr(cols[9], start_index, nchar(cols[9]))

    end_index <- gregexpr('"', gene_name, fixed = TRUE)[[1]][1]
    if (end_index == -1)
    {
      writeLines(line, new_gtf_annotation, useBytes = TRUE)
      next
    }

    gene_name <- substr(gene_name, 1, end_index - 1)

    row_indices <- which(marker_genes_csv$gene == gene_name)
    if (length(row_indices) == 0)
    {
      writeLines(line, new_gtf_annotation, useBytes = TRUE)
      next
    }

    marker_gene_rows <- marker_genes_csv[row_indices, ]
    marker_genes_csv <- marker_genes_csv[-row_indices, ]

    p_val_adjs <- c()
    avg_log2fcs <- c()
    pct.1s <- c()
    pct.2s <- c()
    clusters <- c()
    againsts <- c()
    for (i in 1:nrow(marker_gene_rows)) {
      marker_gene_row <- marker_gene_rows[i,]
      p_val_adj <- marker_gene_row$p_val_adj
      avg_log2fc <- marker_gene_row$avg_log2FC
      pct.1 <- marker_gene_row$pct.1
      pct.2 <- marker_gene_row$pct.2
      cluster <- marker_gene_row$cluster
      against <- marker_gene_row$against
      if (is.null(p_val_adj) || is.null(avg_log2fc) || is.null(pct.1) || is.null(pct.2) || is.null(cluster)) {
        next
      }

      p_val_adjs <- c(p_val_adjs, p_val_adj)
      avg_log2fcs <- c(avg_log2fcs, avg_log2fc)
      pct.1s <- c(pct.1s, pct.1)
      pct.2s <- c(pct.2s, pct.2)
      clusters <- c(clusters, cluster)
      if (!is.null(against)) {
        againsts <- c(againsts, against)
        break
      }
    }

    if (is.null(p_val_adjs) || is.null(avg_log2fcs) || is.null(pct.1s) || is.null(pct.2s) || is.null(clusters)) {
      writeLines(line, new_gtf_annotation, useBytes = TRUE)
      next
    }

    p_val_adjs <- paste(p_val_adjs, collapse = ',')
    avg_log2fcs <- paste(avg_log2fcs, collapse = ',')
    pct.1s <- paste(pct.1s, collapse = ',')
    pct.2s <- paste(pct.2s, collapse = ',')
    clusters <- paste(clusters, collapse = ',')

    cols[9] <- trimws(cols[9])
    if (!endsWith(cols[9], ';')) {
      cols[9] <- paste0(cols[9], ';')
    }

    cols[9] <- paste0(cols[9], ' ',
                      'marker_gene "TRUE"; ',
                      'p_val_adj "', p_val_adjs, '"; ',
                      'avg_log2fc "', avg_log2fcs, '"; ',
                      'pct.1 "', pct.1s, '"; ',
                      'pct.2 "', pct.2s, '"; ',
                      'cluster "', clusters, '"')

    if (!is.null(againsts)) {
      cols[9] <- paste0(cols[9], "; ",
                        'against "', againsts, '"')
    }

    new_line <- paste(cols, collapse = '\t')
    writeLines(new_line, new_gtf_annotation, useBytes = TRUE)
  }

  close(gtf_annotation)
  close(new_gtf_annotation)

  files_to_zip <- "combined_annotations_with_marker_genes.gtf"

  setwd(session_id)
  zipfile_path <- "updated_combined_annotations.zip"

  if (file.exists(zipfile_path)) {
    file.remove(zipfile_path)
  }
  zip(zipfile = zipfile_path, files = files_to_zip)

  setwd("..")
}

# main shiny app server
server <- function(input, output, session) {

  # store session ID
  # create session id tmp directory each time app is run
  session_id <- session$token
  print(paste0("Session: ", session_id))

  # create the dir
  dir.create(session_id)

  # FLAMES MODULE

  # create reactive value for the database zip
  file_available_db <- reactiveVal(FALSE)

  # run database function when submit is pressed
  observeEvent(input$db_submit_button, {
    message("Validating FLAMES module input...")
    error_msg <- is_flames_input_valid(input)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "db-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "db-status-msg-container"))
    message("Validated!")

    # ensure download button remains greyed out (if submit is re-pressed)
    shinyjs::disable("db_download_button")
    shinyjs::runjs("document.getElementById('db_download_button').style.backgroundColor = '#d3d3d3';")
    # disable submit button after it is pressed
    session$sendCustomMessage("disableButton", list(id = "db_submit_button", spinnerId = "db-loading-container"))

    # get the number of files uploaded
    num_files <- nrow(input$user_fastq_files)
    is_multi_sample_mode <- (num_files > 1)
    single_sample_id <- ifelse(is_multi_sample_mode, "not_applicable",
                               gsub("\\.fastq(?:\\.gz)?$", "", input$user_fastq_files$name[1]))

    # run flames
    flames_server(input, is_multi_sample_mode, single_sample_id, session)

    # check if the zip file is created
    if (file.exists(file.path(session_id, "database_results.zip"))) {
      file_available_db(TRUE)
    }
  })

  # enable download once files are available
  observe({
    if (file_available_db()) {
      shinyjs::enable("db_download_button")
      shinyjs::runjs("document.getElementById('db_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "db_submit_button", spinnerId = "db-loading-container")) # re-enable submit button
    }
  })

  # download handler for the database results.zip file
  output$db_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_database_results.zip")
    },
    content = function(file) {
      file.copy(file.path(session_id, "database_results.zip"), file)
    }
  )

  # END FLAMES MODULE

  # PROTEOME DATABASE MODULE

  # create reactive value for the protein database zip
  file_available_prot_db <- reactiveVal(FALSE)

  observeEvent(input$prot_db_submit_button, {
    message("Validating proteome DB module input...")
    error_msg <- is_prot_db_input_valid(input)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "prot-db-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "prot-db-status-msg-container"))
    message("Validated!")

    # ensure download button remains greyed out (if submit is re-pressed)
    shinyjs::disable("prot_db_download_button")
    shinyjs::runjs("document.getElementById('prot_db_download_button').style.backgroundColor = '#d3d3d3';")
    # disable submit button after it is pressed
    session$sendCustomMessage("disableButton", list(id = "prot_db_submit_button", spinnerId = "prot-db-loading-container"))

    # identify ORFs and generate the proteome database
    database_server(input, session)

    # # check if the zip file is created
    if (file.exists(file.path(session_id, "protein_database_results.zip"))) {
      file_available_prot_db(TRUE)
    }
  })

  # enable download once files are available
  observe({
    if (file_available_prot_db()) {
      shinyjs::enable("prot_db_download_button")
      shinyjs::runjs("document.getElementById('prot_db_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "prot_db_submit_button", spinnerId = "prot-db-loading-container")) # re-enable submit button
    }
  })

  # download handler 
  output$prot_db_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_protein_database_results.zip")
    },
    content = function(file) {
      file.copy(file.path(session_id, "protein_database_results.zip"), file)
    }
  )

  # SEURAT FILTERING AND INTEGRATION MODULE

  file_available_seurat <- reactiveVal(FALSE)
  total_cluster_num <- -1
  global_npc_input <- -1
  global_seurat_obj <- FALSE

  remove_cluster_name_inputs <- function(session) {
    session$sendCustomMessage("clearClusterList", list())
  }

  # set output dir
  outdir_seurat <- file.path(session_id, "seurat_integ_output")

  # create output dir
  dir.create(outdir_seurat, recursive = TRUE)

  # enable download once files are available
  observe({
    if (file_available_seurat()) {
      shinyjs::enable("seurat_download_button")
      shinyjs::runjs("document.getElementById('seurat_download_button').style.backgroundColor = '#4CAF50';")
    } else {
      shinyjs::disable("seurat_download_button")
      shinyjs::runjs("document.getElementById('seurat_download_button').style.backgroundColor = '#FFF';")
    }
  })

  seurat_file_to_download <- function() {
    mode <- seurat_data_storage$mode
    if (mode == "view") {
      message("seurat_view_results.zip")
      return("seurat_view_results.zip")
    } else if (mode == "qc_view") {
      message("seurat_qc_view_results.zip")
      return("seurat_qc_view_results.zip")
    }
    message("seurat_integration_results.zip")
    return("seurat_integration_results.zip")
  }

  output$seurat_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_", seurat_file_to_download())
    },
    content = function(file) {
      file.copy(file.path(session_id, seurat_file_to_download()), file)
    }
  )

  show_default_seurat_module()

  # Create some temporary variables for the filtering and integration module
  seurat_data_storage <- reactiveValues()
  seurat_data_storage$mode <- ""

  # Function for updating the single-cell sample selector
  update_sample_selector <- function(session, samples_available) {
    updateSelectInput(session, "sample_selector", choices = samples_available)
  }

  view_stats_function <- function(session, gene_counts_file, transcript_count_file, output) {
    file_available_seurat(FALSE)
    session$sendCustomMessage("clearStatusMessage", list(container = "seurat-status-msg-container"))
    hide("seurat_statistics_msg_container")
    output$seurat_statistics_msg_container <- renderTable(data.frame(c()), rownames = FALSE, colnames = FALSE)
    session$sendCustomMessage("disableButtonOnly", list(id = "view_qc_stats_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "perform_integration_button"))
    session$sendCustomMessage("disableButton", list(id = "view_stats_button", spinnerId = "seurat-loading-container"))
    global_npc_input <<- -1
    global_seurat_obj <<- FALSE
    total_cluster_num <<- -1

    output$scatter_plot <- NULL
    output$vln_plot <- NULL
    output$elbow_plot <- NULL
    output$umap_plot <- NULL
    output$marker_gene_plot <- NULL

    remove_cluster_name_inputs(session)
    hide_elements(seurat_module_outputs)
    hide_post_qc_parts_seurat_module()

    seurat_obj <- get_seurat_obj_from_file(gene_counts_file)

    nCount_RNA <- seurat_obj$nCount_RNA
    quartile5 <- as.integer(floor(quantile(nCount_RNA, probs = c(0.05))))
    quartile95 <- as.integer(ceiling(quantile(nCount_RNA, probs = c(0.95))))
    updateNumericInput(session, "min_nCount_RNA", value = quartile5)
    updateNumericInput(session, "max_nCount_RNA", value = quartile95)

    nFeature_RNA <- seurat_obj$nFeature_RNA
    quartile5 <- as.integer(floor(quantile(nFeature_RNA, probs = c(0.05))))
    quartile95 <- as.integer(ceiling(quantile(nFeature_RNA, probs = c(0.95))))
    updateNumericInput(session, "min_nFeature_RNA", value = quartile5)
    updateNumericInput(session, "max_nFeature_RNA", value = quartile95)

    stats_df <- get_stats_df(seurat_obj)
    output$seurat_statistics_msg_container <- renderTable(stats_df, rownames = TRUE, colnames = TRUE)

    scatter_plot <- get_scatter_plot(seurat_obj)
    output$scatter_plot <- renderPlot({
      suppressWarnings(print(scatter_plot))
    })

    vln_plot <- get_vln_plot(seurat_obj)
    output$vln_plot <- renderPlot({
      suppressWarnings(print(vln_plot))
    })

    top_level_dir <- getwd()
    output_files_dir <- outdir_seurat
    setwd(output_files_dir)

    stats_msg_file <- write_stats_msg(seurat_obj)
    scatter_plot_file <- write_scatter_plot(seurat_obj)
    vln_plot_file <- write_vln_plot(seurat_obj)

    seurat_file <- "seurat_object.rds"
    saveRDS(seurat_obj, file = seurat_file)

    files_to_zip <- c(stats_msg_file, scatter_plot_file, vln_plot_file, seurat_file)
    zipfile_path <- file.path("..", seurat_file_to_download())

    if (file.exists(zipfile_path)) {
      file.remove(zipfile_path)
    }
    zip(zipfile = zipfile_path, files = files_to_zip)
    setwd(top_level_dir)

    show_basic_outputs_seurat_module()

    session$sendCustomMessage("enableButtonOnly", list(id = "view_qc_stats_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "perform_integration_button"))
    session$sendCustomMessage("enableButton", list(id = "view_stats_button", spinnerId = "seurat-loading-container"))

    if (file.exists(file.path(session_id, seurat_file_to_download()))) {
      file_available_seurat(TRUE)
    }
  }

  qc_stats_function <- function(session, gene_counts_file, transcript_count_file, input, output, session_id) {
    file_available_seurat(FALSE)
    session$sendCustomMessage("clearStatusMessage", list(container = "seurat-status-msg-container"))
    hide("seurat_statistics_msg_container")
    output$seurat_statistics_msg_container <- renderTable(data.frame(c()), rownames = FALSE, colnames = FALSE)
    session$sendCustomMessage("disableButtonOnly", list(id = "view_stats_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "perform_integration_button"))
    session$sendCustomMessage("disableButton", list(id = "view_qc_stats_button", spinnerId = "seurat-loading-container"))
    global_npc_input <<- -1
    global_seurat_obj <<- FALSE
    total_cluster_num <<- -1

    output$scatter_plot <- NULL
    output$vln_plot <- NULL
    output$elbow_plot <- NULL
    output$umap_plot <- NULL
    output$marker_gene_plot <- NULL

    remove_cluster_name_inputs(session)
    hide_elements(seurat_module_outputs)
    hide_post_qc_parts_seurat_module()

    min_nCount_RNA <- floor(input$min_nCount_RNA)
    max_nCount_RNA <- floor(input$max_nCount_RNA)
    min_nFeature_RNA <- floor(input$min_nFeature_RNA)
    max_nFeature_RNA <- floor(input$max_nFeature_RNA)
    max_percent_mt <- input$max_percent_mt
    npc_input <- floor(input$npc_input)
    k_input <- floor(input$k_input)
    cluster_res_input <- input$cluster_res_input

    qc_seurat_obj <- perform_qc_from_file(gene_counts_file, min_nCount_RNA, max_nCount_RNA, min_nFeature_RNA, max_nFeature_RNA, max_percent_mt, npc_input, k_input, cluster_res_input)
    if (isFALSE(qc_seurat_obj)) {
      msg <- paste0("Error: Less than 2 cells remained after filtering '", basename(gene_counts_file), "'. Quality control is not possible. Please use less stringent quality control parameters.")
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      session$sendCustomMessage("enableButtonOnly", list(id = "view_stats_button"))
      session$sendCustomMessage("enableButtonOnly", list(id = "perform_integration_button"))
      session$sendCustomMessage("enableButton", list(id = "view_qc_stats_button", spinnerId = "seurat-loading-container"))
      return()
    }
    global_npc_input <<- npc_input

    # Add transcript count information into the Seurat object
    tx_counts <- fread(file = transcript_count_file)

    # Turn any missing transcript count into 0
    for (i in names(tx_counts)) {
      tx_counts[is.na(get(i)), (i):=0]
    }

    tx_counts <- as.data.frame(tx_counts)
    rownames(tx_counts) <- tx_counts$V1
    tx_counts <- subset(tx_counts, select = -c(V1))
    tx_counts <- tx_counts[rowSums(tx_counts) != 0, ]

    iso_seurat_obj <- CreateSeuratObject(counts = tx_counts)

    # Make sure both the gene counts and transcript counts Seurat objects have the same cells
    iso_seurat_obj <- subset(iso_seurat_obj, cells = Cells(qc_seurat_obj))
    iso_seurat_obj <- JoinLayers(iso_seurat_obj)

    counts_table_iso <- iso_seurat_obj[["RNA"]]$counts
    df_iso <- as.data.frame(counts_table_iso)

    # Remove rows where the sum is 0
    df_iso <- df_iso[rowSums(df_iso) != 0, ]

    qc_seurat_obj[["iso"]] <- CreateAssay5Object(counts = df_iso)

    qc_seurat_obj <- NormalizeData(qc_seurat_obj, assay = "iso")
    qc_seurat_obj <- FindVariableFeatures(qc_seurat_obj, assay = "iso")
    qc_seurat_obj <- ScaleData(qc_seurat_obj, assay = "iso")

    global_seurat_obj <<- qc_seurat_obj

    stats_df <- get_stats_df(qc_seurat_obj)
    output$seurat_statistics_msg_container <- renderTable(stats_df, rownames = TRUE, colnames = TRUE)

    scatter_plot <- get_scatter_plot(qc_seurat_obj)
    output$scatter_plot <- renderPlot({
      suppressWarnings(print(scatter_plot))
    })

    vln_plot <- get_vln_plot(qc_seurat_obj)
    output$vln_plot <- renderPlot({
      suppressWarnings(print(vln_plot))
    })

    elbow_plot <- get_elbow_plot(qc_seurat_obj, npc_input)
    output$elbow_plot <- renderPlot({
      suppressWarnings(print(elbow_plot))
    })

    umap_plot <- get_umap_plot(qc_seurat_obj)
    output$umap_plot <- renderPlot({
      suppressWarnings(print(umap_plot))
    })

    top_level_dir <- getwd()
    output_files_dir <- outdir_seurat
    setwd(output_files_dir)

    cluster_names <- levels(qc_seurat_obj)
    for (i in seq_along(cluster_names)) {
      cluster_name <- as.character(cluster_names[i])
      insertUI(
        selector = '#cluster_list',
        ui = textInput(inputId = paste0("cluster_name_", i), label = paste0("Rename cluster name '", cluster_name, "' to..."), placeholder = "Insert new cluster name here")
      )
    }
    total_cluster_num <<- length(cluster_names)

    # extract cell clusters
    metadata_seurat <- qc_seurat_obj@meta.data
    metadata_seurat$sample <- as.character(metadata_seurat$orig.ident)
    metadata_seurat$barcode <- row.names(metadata_seurat)
    metadata_seurat$cluster <- as.character(Idents(qc_seurat_obj))

    metadata_seurat <- metadata_seurat %>% dplyr::select(sample, barcode, cluster)
    row.names(metadata_seurat) <- NULL

    # export
    sample_barcode_cluster_file <- "sample_barcode_cluster.txt"
    write.table(metadata_seurat, file = sample_barcode_cluster_file, sep = "\t", quote = FALSE)

    stats_msg_file <- write_stats_msg(qc_seurat_obj)
    scatter_plot_file <- write_scatter_plot(qc_seurat_obj)
    vln_plot_file <- write_vln_plot(qc_seurat_obj)
    elbow_plot_file <- write_elbow_plot(qc_seurat_obj, npc_input)
    umap_plot_file <- write_umap_plot(qc_seurat_obj)

    seurat_file <- "seurat_object.rds"
    saveRDS(qc_seurat_obj, file = seurat_file)

    files_to_zip <- c(stats_msg_file, scatter_plot_file, vln_plot_file, elbow_plot_file, umap_plot_file, seurat_file, sample_barcode_cluster_file)
    zipfile_path <- file.path("..", seurat_file_to_download())

    if (file.exists(zipfile_path)) {
      file.remove(zipfile_path)
    }
    zip(zipfile = zipfile_path, files = files_to_zip)
    setwd(top_level_dir)

    show_more_outputs_seurat_module()

    session$sendCustomMessage("enableButtonOnly", list(id = "view_stats_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "perform_integration_button"))
    session$sendCustomMessage("enableButton", list(id = "view_qc_stats_button", spinnerId = "seurat-loading-container"))

    if (file.exists(file.path(session_id, seurat_file_to_download()))) {
      file_available_seurat(TRUE)
    }
  }

  integ_samples_function <- function(session, gene_counts_files, transcript_counts_files, input, output, session_id) {
    file_available_seurat(FALSE)
    session$sendCustomMessage("clearStatusMessage", list(container = "seurat-status-msg-container"))
    hide("seurat_statistics_msg_container")
    output$seurat_statistics_msg_container <- renderTable(data.frame(c()), rownames = FALSE, colnames = FALSE)
    session$sendCustomMessage("disableButtonOnly", list(id = "view_stats_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "view_qc_stats_button"))
    session$sendCustomMessage("disableButton", list(id = "perform_integration_button", spinnerId = "seurat-loading-container"))
    global_npc_input <<- -1
    global_seurat_obj <<- FALSE
    total_cluster_num <<- -1

    output$scatter_plot <- NULL
    output$vln_plot <- NULL
    output$elbow_plot <- NULL
    output$umap_plot <- NULL
    output$marker_gene_plot <- NULL

    remove_cluster_name_inputs(session)
    hide_elements(seurat_module_outputs)
    hide_post_qc_parts_seurat_module()

    min_nCount_RNA <- floor(input$min_nCount_RNA)
    max_nCount_RNA <- floor(input$max_nCount_RNA)
    min_nFeature_RNA <- floor(input$min_nFeature_RNA)
    max_nFeature_RNA <- floor(input$max_nFeature_RNA)
    max_percent_mt <- input$max_percent_mt
    npc_input <- floor(input$npc_input)
    k_input <- floor(input$k_input)
    cluster_res_input <- input$cluster_res_input

    # Perform quality control on each sample
    sample_names <- c()

    gene_counts_file <- gene_counts_files[1]
    qc_seurat_obj <- perform_qc_from_file(gene_counts_file, min_nCount_RNA, max_nCount_RNA, min_nFeature_RNA, max_nFeature_RNA, max_percent_mt, npc_input, k_input, cluster_res_input)
    if (isFALSE(qc_seurat_obj)) {
      msg <- paste0("Error: Less than 2 cells remained after filtering '",  basename(gene_counts_file), "'. Quality control is not possible. Please use less stringent quality control parameters.")
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      session$sendCustomMessage("enableButtonOnly", list(id = "view_stats_button"))
      session$sendCustomMessage("enableButtonOnly", list(id = "view_qc_stats_button"))
      session$sendCustomMessage("enableButton", list(id = "perform_integration_button", spinnerId = "seurat-loading-container"))
      return()
    }

    sample_name <- gsub("_gene_count\\.csv$", "", basename(gene_counts_file))
    message(paste0(sample_name, " is processed!"))
    first_umap_object <- qc_seurat_obj
    sample_names <- c(sample_names, sample_name)

    umap_objects <- list()
    for (i in seq_along(gene_counts_files)) {
      if (i == 1) {
        next
      }
      gene_counts_file <- gene_counts_files[i]
      qc_seurat_obj <- perform_qc_from_file(gene_counts_file, min_nCount_RNA, max_nCount_RNA, min_nFeature_RNA, max_nFeature_RNA, max_percent_mt, npc_input, k_input, cluster_res_input)
      if (isFALSE(qc_seurat_obj)) {
        msg <- paste0("Error: Less than 2 cells remained after filtering '",  basename(gene_counts_file), "'. Quality control is not possible. Please use less stringent quality control parameters.")
        session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
        session$sendCustomMessage("enableButtonOnly", list(id = "view_stats_button"))
        session$sendCustomMessage("enableButtonOnly", list(id = "view_qc_stats_button"))
        session$sendCustomMessage("enableButton", list(id = "perform_integration_button", spinnerId = "seurat-loading-container"))
        return()
      }

      sample_name <- gsub("_gene_count\\.csv$", "", basename(gene_counts_file))
      message(paste0(sample_name, " is processed!"))
      umap_objects <- base::append(umap_objects, qc_seurat_obj)
      sample_names <- c(sample_names, sample_name)
    }

    # Merge the quality-controlled samples together
    merged_seurat <- merge(first_umap_object, y = umap_objects, add.cell.ids = sample_names, project = "Integrated")

    merged_seurat <- NormalizeData(object = merged_seurat)
    merged_seurat <- FindVariableFeatures(object = merged_seurat)
    merged_seurat <- ScaleData(object = merged_seurat)
    merged_seurat <- RunPCA(object = merged_seurat)
    merged_seurat <- FindNeighbors(object = merged_seurat, dims = 1:16)
    merged_seurat <- FindClusters(object = merged_seurat, resolution = 0.6)
    merged_seurat <- RunUMAP(object = merged_seurat, dims = 1:30, seed.use = 42)

    merged_seurat <- JoinLayers(object = merged_seurat)
    merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]], f = merged_seurat$orig.ident)

    # Perform integration with Harmony

    npcs <- floor(input$npcs)
    theta <- input$theta
    lambda <- input$lambda
    sigma <- input$sigma
    nclust <- floor(input$nclust)
    if (nclust == -1) {
      nclust <- NULL
    }
    tau <- floor(input$tau)
    block.size <- input$block.size
    max.iter.harmony <- floor(input$max.iter.harmony)
    max.iter.cluster <- floor(input$max.iter.cluster)
    epsilon.cluster <- ifelse(is.finite(input$epsilon.cluster) && (input$epsilon.cluster > 0), input$epsilon.cluster, NULL)
    epsilon.harmony <- ifelse(is.finite(input$epsilon.harmony) && (input$epsilon.harmony > 0), input$epsilon.harmony, NULL)

    obj <- IntegrateLayers(object = merged_seurat,
                           method = HarmonyIntegration,
                           orig.reduction = "pca",
                           new.reduction = 'integrated.harm',
                           npcs = npcs, theta = theta, lambda = lambda, sigma = sigma, nclust = nclust, tau = tau, block.size = block.size, max.iter.harmony = max.iter.harmony,
                           max.iter.cluster = max.iter.cluster, epsilon.cluster = epsilon.cluster, epsilon.harmony = epsilon.harmony,
                           verbose = TRUE)
    obj <- FindNeighbors(obj, reduction="integrated.harm")
    obj <- FindClusters(obj, resolution=0.4, cluster.name="harm_cluster")
    obj <- RunUMAP(obj, reduction="integrated.harm", dims=1:20, reduction.name = "umap.harm", seed.use = 42)

    global_npc_input <<- npc_input

    # Perform quality control on each sample
    sample_names <- c()

    transcript_counts_file <- transcript_counts_files[1]

    tx_counts <- fread(file = transcript_counts_file)

    # Turn any missing transcript count into 0
    for (i in names(tx_counts)) {
      tx_counts[is.na(get(i)), (i):=0]
    }

    tx_counts <- as.data.frame(tx_counts)
    rownames(tx_counts) <- tx_counts$V1
    tx_counts <- subset(tx_counts, select = -c(V1))
    iso_seurat_obj <- CreateSeuratObject(counts = tx_counts)

    sample_name <- gsub("_transcript_count\\.csv$", "", basename(transcript_counts_file))
    message(paste0(sample_name, " is processed!"))
    first_iso_object <- iso_seurat_obj
    sample_names <- c(sample_names, sample_name)

    iso_objects <- list()
    for (i in seq_along(transcript_counts_files)) {
      if (i == 1) {
        next
      }
      transcript_counts_file <- transcript_counts_files[i]

      tx_counts <- fread(file = transcript_counts_file)

      # Turn any missing transcript count into 0
      for (i in names(tx_counts)) {
        tx_counts[is.na(get(i)), (i):=0]
      }

      tx_counts <- as.data.frame(tx_counts)
      rownames(tx_counts) <- tx_counts$V1
      tx_counts <- subset(tx_counts, select = -c(V1))
      iso_seurat_obj <- CreateSeuratObject(counts = tx_counts)

      sample_name <- gsub("_transcript_count\\.csv$", "", basename(transcript_counts_file))
      message(paste0(sample_name, " is processed!"))
      iso_objects <- base::append(iso_objects, iso_seurat_obj)
      sample_names <- c(sample_names, sample_name)
    }

    merged_iso_obj <- merge(first_iso_object, y = iso_objects, add.cell.ids = sample_names, project = "Integrated")
    merged_iso_obj <- subset(merged_iso_obj, cells = obj@graphs[["RNA_nn"]]@Dimnames[[1]])
    merged_iso_obj <- JoinLayers(merged_iso_obj)
    counts_table_iso <- merged_iso_obj[["RNA"]]$counts
    df_iso <- as.data.frame(counts_table_iso)

    # Remove rows where the sum is 0
    df_iso <- df_iso[rowSums(df_iso) != 0, ]

    obj[["iso"]] <- CreateAssay5Object(counts = df_iso)

    obj <- NormalizeData(obj, assay = "iso")
    obj <- FindVariableFeatures(obj, assay = "iso")
    obj <- ScaleData(obj, assay = "iso")

    global_seurat_obj <<- obj

    stats_df <- get_stats_df(obj)
    output$seurat_statistics_msg_container <- renderTable(stats_df, rownames = TRUE, colnames = TRUE)

    scatter_plot <- get_scatter_plot(obj)
    output$scatter_plot <- renderPlot({
      suppressWarnings(print(scatter_plot))
    })

    vln_plot <- get_vln_plot(obj)
    output$vln_plot <- renderPlot({
      suppressWarnings(print(vln_plot))
    })

    elbow_plot <- get_elbow_plot(qc_seurat_obj, npc_input)
    output$elbow_plot <- renderPlot({
      suppressWarnings(print(elbow_plot))
    })

    umap_plot <- get_umap_harm_plot(obj)
    output$umap_plot <- renderPlot({
      suppressWarnings(print(umap_plot))
    })

    top_level_dir <- getwd()
    output_files_dir <- outdir_seurat
    setwd(output_files_dir)

    cluster_names <- levels(obj@meta.data$harm_cluster)
    message(cluster_names)
    for (i in seq_along(cluster_names)) {
      cluster_name <- as.character(cluster_names[i])
      insertUI(
        selector = '#cluster_list',
        ui = textInput(inputId = paste0("cluster_name_", i), label = paste0("Rename cluster name '", cluster_name, "' to..."), placeholder = "Insert new cluster name here")
      )
    }
    total_cluster_num <<- length(cluster_names)

    # extract cell clusters
    metadata_seurat <- obj@meta.data
    metadata_seurat$sample <- as.character(metadata_seurat$orig.ident)
    metadata_seurat$barcode <- row.names(metadata_seurat)
    metadata_seurat$cluster <- metadata_seurat$harm_cluster

    metadata_seurat <- metadata_seurat %>% dplyr::select(sample, barcode, cluster)
    row.names(metadata_seurat) <- NULL

    # export
    sample_barcode_cluster_file <- "sample_barcode_cluster.txt"
    write.table(metadata_seurat, file = sample_barcode_cluster_file, sep = "\t", quote = FALSE)

    stats_msg_file <- write_stats_msg(obj)
    scatter_plot_file <- write_scatter_plot(obj)
    vln_plot_file <- write_vln_plot(obj)
    elbow_plot_file <- write_elbow_plot(qc_seurat_obj, npc_input)
    umap_harm_plot_file <- write_umap_harm_plot(obj)

    seurat_file <- "seurat_object.rds"
    saveRDS(obj, file = seurat_file)

    files_to_zip <- c(stats_msg_file, scatter_plot_file, vln_plot_file, elbow_plot_file, umap_harm_plot_file, seurat_file, sample_barcode_cluster_file)
    zipfile_path <- file.path("..", seurat_file_to_download())

    if (file.exists(zipfile_path)) {
      file.remove(zipfile_path)
    }
    zip(zipfile = zipfile_path, files = files_to_zip)
    setwd(top_level_dir)

    show_more_outputs_seurat_module()

    session$sendCustomMessage("enableButtonOnly", list(id = "view_stats_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "view_qc_stats_button"))
    session$sendCustomMessage("enableButton", list(id = "perform_integration_button", spinnerId = "seurat-loading-container"))

    if (file.exists(file.path(session_id, seurat_file_to_download()))) {
      file_available_seurat(TRUE)
    }
  }

  observeEvent(input$show_qc_part, {do_show_qc_part_action()})
  observeEvent(input$hide_qc_part, {do_hide_qc_part_action()})
  observeEvent(input$show_integ_part, {do_show_integ_part_action()})
  observeEvent(input$hide_integ_part, {do_hide_integ_part_action()})
  observeEvent(input$show_rename_clusters_part, {do_show_rename_clusters_part_action()})
  observeEvent(input$hide_rename_clusters_part, {do_hide_rename_clusters_part_action()})
  observeEvent(input$show_marker_genes_part, {do_show_marker_genes_part_action()})
  observeEvent(input$hide_marker_genes_part, {do_hide_marker_genes_part_action()})
  observeEvent(input$show_all_marker_genes_part, {do_show_all_marker_genes_part_action()})
  observeEvent(input$hide_all_marker_genes_part, {do_hide_all_marker_genes_part_action()})
  observeEvent(input$show_adv_qc_params_button, {do_show_adv_qc_params_action()})
  observeEvent(input$hide_adv_qc_params_button, {do_hide_adv_qc_params_action()})
  observeEvent(input$show_adv_harmony_integ_params_button, {do_show_adv_harmony_integ_params_action()})
  observeEvent(input$hide_adv_harmony_integ_params_button, {do_hide_adv_harmony_integ_params_action()})
  observeEvent(input$show_statistics_output_part, {do_show_statistics_output_part_action()})
  observeEvent(input$hide_statistics_output_part, {do_hide_statistics_output_part_action()})

  # sample selection in dropdown menu
  observeEvent(input$sample_selector, {
    sample_id <- input$sample_selector
    mode <- seurat_data_storage$mode
    req(sample_id, mode)

    # Get the path to the sample's gene count CSV
    gene_count_file <- file.path(session_id, "sample_files", paste0(sample_id, "_gene_count.csv"))

    # If the file doesn't exist, something's wrong
    if (!file.exists(gene_count_file)) {
      error_msg <- paste0("Gene count file for sample '", sample_id, "' does not exist! Please contact the GenomeProtSC team to report this bug.")
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    # Get the path to the sample's transcript count CSV
    transcript_count_file <- file.path(session_id, "sample_files", paste0(sample_id, "_transcript_count.csv"))

    # If the file doesn't exist, something's wrong
    if (!file.exists(transcript_count_file)) {
      error_msg <- paste0("Transcript count file for sample '", sample_id, "' does not exist! Please contact the GenomeProtSC team to report this bug.")
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    if (mode == "view") {
      view_stats_function(session, gene_count_file, transcript_count_file, output)
    } else if (mode == "qc_view") {
      qc_stats_function(session, gene_count_file, transcript_count_file, input, output, session_id)
    }
  })

  observeEvent(input$view_stats_button, {
    message("Validating the view stats part of the integration module input...")
    error_msg <- is_view_stats_input_valid(input)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "seurat-status-msg-container"))
    message("Validated!")

    # Create a directory to store the sample files
    session_id <- session$token
    sample_files_dir <- file.path(session_id, "sample_files")
    if (!dir.exists(sample_files_dir)) {
      dir.create(sample_files_dir, recursive = TRUE)
    }

    # Rename the temporary gene count files to their original filenames
    gene_count_files <- input$gene_count_files
    renamed_files <- file.path(dirname(gene_count_files$datapath), gene_count_files$name)

    if (all(mapply(file.exists, gene_count_files$datapath))) {
      mapply(file.rename, from = gene_count_files$datapath, to = renamed_files)
    }

    # Copy the gene count files over to the created directory
    if (all(mapply(file.exists, renamed_files))) {
      file.copy(renamed_files, sample_files_dir)
    }

    # Rename the temporary transcript count files to their original filenames
    transcript_count_files <- input$transcript_count_files
    renamed_files <- file.path(dirname(transcript_count_files$datapath), transcript_count_files$name)

    if (all(mapply(file.exists, transcript_count_files$datapath))) {
      mapply(file.rename, from = transcript_count_files$datapath, to = renamed_files)
    }

    # Copy the transcript count files over to the created directory
    if (all(mapply(file.exists, renamed_files))) {
      file.copy(renamed_files, sample_files_dir)
    }

    seurat_data_storage$mode <- "view"
    samples_available <- gsub("_gene_count\\.csv$", "", input$gene_count_files$name)

    update_sample_selector(session, samples_available)

    # TODO: Fix the 'gene selector update' event firing twice
    sample_id <- samples_available[1]
    gene_count_file <- file.path(session_id, "sample_files", paste0(sample_id, "_gene_count.csv"))
    transcript_count_file <- file.path(session_id, "sample_files", paste0(sample_id, "_transcript_count.csv"))
    view_stats_function(session, gene_count_file, transcript_count_file, output)
  })

  observeEvent(input$view_qc_stats_button, {
    message("Validating the QC view stats part of the integration module input...")
    error_msg <- is_qc_view_stats_input_valid(input)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "seurat-status-msg-container"))
    message("Validated!")

    # Create a directory to store the sample files
    session_id <- session$token
    sample_files_dir <- file.path(session_id, "sample_files")
    if (!dir.exists(sample_files_dir)) {
      dir.create(sample_files_dir, recursive = TRUE)
    }

    # Rename the temporary gene count files to their original filenames
    gene_count_files <- input$gene_count_files
    renamed_files <- file.path(dirname(gene_count_files$datapath), gene_count_files$name)

    if (all(mapply(file.exists, gene_count_files$datapath))) {
      mapply(file.rename, from = gene_count_files$datapath, to = renamed_files)
    }

    # Copy the gene count files over to the created directory
    if (all(mapply(file.exists, renamed_files))) {
      file.copy(renamed_files, sample_files_dir)
    }

    # Rename the temporary transcript count files to their original filenames
    transcript_count_files <- input$transcript_count_files
    renamed_files <- file.path(dirname(transcript_count_files$datapath), transcript_count_files$name)

    if (all(mapply(file.exists, transcript_count_files$datapath))) {
      mapply(file.rename, from = transcript_count_files$datapath, to = renamed_files)
    }

    # Copy the transcript count files over to the created directory
    if (all(mapply(file.exists, renamed_files))) {
      file.copy(renamed_files, sample_files_dir)
    }

    seurat_data_storage$mode <- "qc_view"
    samples_available <- gsub("_gene_count\\.csv$", "", input$gene_count_files$name)

    update_sample_selector(session, samples_available)

    # TODO: Fix the 'gene selector update' event firing twice
    sample_id <- samples_available[1]
    gene_count_file <- file.path(session_id, "sample_files", paste0(sample_id, "_gene_count.csv"))
    transcript_count_file <- file.path(session_id, "sample_files", paste0(sample_id, "_transcript_count.csv"))
    qc_stats_function(session, gene_count_file, transcript_count_file, input, output, session_id)
  })

  observeEvent(input$perform_integration_button, {
    message("Validating the 'perform integration' part of the integration module input...")
    error_msg <- is_perform_integ_input_valid(input)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "seurat-status-msg-container"))
    message("Validated!")

    # Create a directory to store the sample files
    session_id <- session$token
    sample_files_dir <- file.path(session_id, "sample_files")
    if (!dir.exists(sample_files_dir)) {
      dir.create(sample_files_dir, recursive = TRUE)
    }

    # Rename the temporary gene count files to their original filenames
    gene_count_files <- input$gene_count_files
    renamed_gene_files <- file.path(dirname(gene_count_files$datapath), gene_count_files$name)

    if (all(mapply(file.exists, gene_count_files$datapath))) {
      mapply(file.rename, from = gene_count_files$datapath, to = renamed_gene_files)
    }

    # Copy the gene count files over to the created directory
    if (all(mapply(file.exists, renamed_gene_files))) {
      file.copy(renamed_gene_files, sample_files_dir)
    }

    # Rename the temporary transcript count files to their original filenames
    transcript_count_files <- input$transcript_count_files
    renamed_transcript_files <- file.path(dirname(transcript_count_files$datapath), transcript_count_files$name)

    if (all(mapply(file.exists, transcript_count_files$datapath))) {
      mapply(file.rename, from = transcript_count_files$datapath, to = renamed_transcript_files)
    }

    # Copy the transcript count files over to the created directory
    if (all(mapply(file.exists, renamed_transcript_files))) {
      file.copy(renamed_transcript_files, sample_files_dir)
    }

    seurat_data_storage$mode <- "perform_integ"

    integ_samples_function(session, renamed_gene_files, renamed_transcript_files, input, output, session_id)
  })

  observeEvent(input$rename_clusters_button, {
    if (((seurat_data_storage$mode != "qc_view") && (seurat_data_storage$mode != "perform_integ")) || (total_cluster_num <= 0) || (global_npc_input <= 1) || isFALSE(global_seurat_obj)) {
      msg <- "Error: This error message should not appear! Please report this issue to the GenomeProtSC team."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    new_cluster_names <- c()
    for (i in 1:total_cluster_num) {
      new_cluster_name <- trimws(input[[paste0("cluster_name_", i)]])
      if (new_cluster_name == "") {
        msg <- "Error: At least one of the new cluster names is empty or consists entirely of whitespace. Please ensure new cluster names contain alphanumeric characters, no commas, and no double quotes."
        session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
        return()
      }
      if (grepl(",", new_cluster_name, fixed = TRUE)) {
        msg <- "Error: At least one of the new cluster names contains commas. Please ensure new cluster names do not contain commas."
        session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
        return()
      }
      if (grepl('"', new_cluster_name, fixed = TRUE)) {
        msg <- "Error: At least one of the new cluster names contains double quotes. Please ensure new cluster names do not contain double quotes."
        session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
        return()
      }
      new_cluster_names <- c(new_cluster_names, new_cluster_name)
    }

    if (length(unique(new_cluster_names)) != total_cluster_num) {
      msg <- "Error: There are duplicate cluster names. Please ensure new cluster names are unique."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    file_available_seurat(FALSE)
    session$sendCustomMessage("clearStatusMessage", list(container = "seurat-status-msg-container"))
    session$sendCustomMessage("disableButtonOnly", list(id = "view_stats_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "view_qc_stats_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "perform_integration_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "find_marker_genes_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "find_all_marker_genes_button"))
    session$sendCustomMessage("disableButton", list(id = "rename_clusters_button", spinnerId = "seurat-loading-container"))
    remove_cluster_name_inputs(session)

    output$umap_plot <- NULL
    output$marker_gene_plot <- NULL

    hide_elements(seurat_module_rename_clusters_parts)
    hide_elements(c("umap_plot_heading", "umap_plot"))
    hide_marker_gene_outputs_seurat_module()

    if (seurat_data_storage$mode == "qc_view") {
      names(new_cluster_names) <- levels(global_seurat_obj)
      global_seurat_obj <<- RenameIdents(global_seurat_obj, new_cluster_names)
    }
    else if (seurat_data_storage$mode == "perform_integ") {
      levels(global_seurat_obj@meta.data$harm_cluster) <<- new_cluster_names
    }

    ### Update what's being shown

    if (seurat_data_storage$mode == "qc_view") {
      umap_plot <- get_umap_plot(global_seurat_obj)
    }
    else if (seurat_data_storage$mode == "perform_integ") {
      umap_plot <- get_umap_harm_plot(global_seurat_obj)
    }

    output$umap_plot <- renderPlot({
      suppressWarnings(print(umap_plot))
    })

    top_level_dir <- getwd()
    output_files_dir <- outdir_seurat
    setwd(output_files_dir)

    if (seurat_data_storage$mode == "qc_view") {
      umap_plot_file <- write_umap_plot(global_seurat_obj)
    }
    else if (seurat_data_storage$mode == "perform_integ") {
      umap_plot_file <- write_umap_harm_plot(global_seurat_obj)
    }

    seurat_file <- "seurat_object.rds"
    saveRDS(global_seurat_obj, file = seurat_file)

    for (i in seq_along(new_cluster_names)) {
      cluster_name <- as.character(new_cluster_names[i])
      insertUI(
        selector = '#cluster_list',
        ui = textInput(inputId = paste0("cluster_name_", i), label = paste0("Rename cluster name '", cluster_name, "' to..."), placeholder = "Insert new cluster name here")
      )
    }

    files_to_zip <- c(umap_plot_file, seurat_file)
    zipfile_path <- file.path("..", seurat_file_to_download())

    zip(zipfile = zipfile_path, files = files_to_zip)
    setwd(top_level_dir)

    show_rename_clusters_part_and_umap_outputs_seurat_module()

    session$sendCustomMessage("enableButtonOnly", list(id = "view_stats_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "view_qc_stats_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "perform_integration_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "find_marker_genes_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "find_all_marker_genes_button"))
    session$sendCustomMessage("enableButton", list(id = "rename_clusters_button", spinnerId = "seurat-loading-container"))

    if (file.exists(file.path(session_id, seurat_file_to_download()))) {
      file_available_seurat(TRUE)
    }
  })

  observeEvent(input$find_marker_genes_button, {
    if (((seurat_data_storage$mode != "qc_view") && (seurat_data_storage$mode != "perform_integ")) || (total_cluster_num <= 0) || (global_npc_input <= 1) || isFALSE(global_seurat_obj)) {
      msg <- "Error: This error message should not appear! Please report this issue to the GenomeProtSC team."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    ident.1 <- trimws(input$ident.1)
    ident.2 <- trimws(input$ident.2)
    p_adj.thresh <- input$p_adj.thresh
    logfc.threshold <- input$logfc.threshold
    min.pct <- input$min.pct
    min.diff.pct <- input$min.diff.pct
    only.pos <- input$only.pos

    if (!(is.finite(p_adj.thresh) && (p_adj.thresh > 0) && (p_adj.thresh <= 1))) {
      msg <- "Error: p_adj.thresh must be a positive number no larger than 1."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    if (!(is.finite(logfc.threshold) && (logfc.threshold > 0))) {
      msg <- "Error: logfc.threshold must be a positive number."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    if (!(is.finite(min.pct) && (min.pct > 0) && (min.pct <= 1))) {
      msg <- "Error: min.pct must be a positive number no larger than 1."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    if (!(is.finite(min.diff.pct) && (((min.diff.pct > 0) && (min.diff.pct <= 1)) || (min.diff.pct < 0)))) {
      msg <- "Error: min.diff.pct must either be a positive number no larger than 1 or a negative number."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    if (min.diff.pct < 0) {
      min.diff.pct = -Inf
    }

    cluster_names <- as.character(Idents(global_seurat_obj))

    # TODO: Properly support Harmony clusters

    if (grepl(",", ident.1, fixed = TRUE)) {
      ident.1 <- unique(trimws(base::strsplit(ident.1, ",")[[1]]))
    }

    if (any(!(ident.1 %in% cluster_names))) {
      msg <- "Error: ident.1 must contain existing cluster names."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    if (ident.2 == "") {
      ident.2 = NULL
    } else if (grepl(",", ident.2, fixed = TRUE)) {
      ident.2 <- unique(trimws(base::strsplit(ident.2, ",")[[1]]))
    }

    if (!is.null(ident.2) && any(!(ident.2 %in% cluster_names))) {
      msg <- "Error: ident.2 must be either contain existing cluster names or be empty."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    if (!is.null(ident.2) && !identical(base::intersect(ident.1, ident.2), character(0))) {
      msg <- "Error: ident.1 and ident.2 must each contain unique cluster names."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    file_available_seurat(FALSE)
    session$sendCustomMessage("clearStatusMessage", list(container = "seurat-status-msg-container"))
    session$sendCustomMessage("disableButtonOnly", list(id = "view_stats_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "view_qc_stats_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "perform_integration_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "rename_clusters_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "find_all_marker_genes_button"))
    session$sendCustomMessage("disableButton", list(id = "find_marker_genes_button", spinnerId = "seurat-loading-container"))

    output$marker_gene_plot <- NULL
    hide_marker_gene_outputs_seurat_module()

    markers <- FindMarkers(global_seurat_obj, ident.1 = ident.1, ident.2 = ident.2, random.seed = 1, logfc.threshold = logfc.threshold,
                           min.pct = min.pct, min.diff.pct = min.diff.pct, only.pos = only.pos)
    markers$gene <- rownames(markers)

    # Apply the adjusted p-value threshold
    markers <- markers %>% dplyr::filter(p_val_adj < p_adj.thresh)

    top20 <- markers$gene[1:20]

    marker_gene_plot <- get_marker_gene_plot(global_seurat_obj, top20)
    output$marker_gene_plot <- renderPlot({
      suppressWarnings(print(marker_gene_plot))
    })
    show_marker_gene_outputs_seurat_module()

    ### Update the zip file

    top_level_dir <- getwd()
    output_files_dir <- outdir_seurat
    setwd(output_files_dir)

    # Indicate which clusters the marker genes are found from
    markers$cluster <- paste(ident.1, collapse = ",")

    # If not comparing against every other cluster, specify which clusters are compared against
    if (!is.null(ident.2)) {
      markers$against <- paste(ident.2, collapse = ",")
    }

    marker_genes_file <- "marker_genes.csv"
    write.csv(markers, marker_genes_file, row.names = FALSE)

    marker_gene_plot_file <- write_marker_gene_plot(global_seurat_obj, top20)

    files_to_zip <- c(marker_genes_file, marker_gene_plot_file)
    zipfile_path <- file.path("..", seurat_file_to_download())

    zip(zipfile = zipfile_path, files = files_to_zip)
    setwd(top_level_dir)

    session$sendCustomMessage("enableButtonOnly", list(id = "view_stats_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "view_qc_stats_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "perform_integration_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "rename_clusters_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "find_all_marker_genes_button"))
    session$sendCustomMessage("enableButton", list(id = "find_marker_genes_button", spinnerId = "seurat-loading-container"))

    if (file.exists(file.path(session_id, seurat_file_to_download()))) {
      file_available_seurat(TRUE)
    }
  })

  observeEvent(input$find_all_marker_genes_button, {
    if (((seurat_data_storage$mode != "qc_view") && (seurat_data_storage$mode != "perform_integ")) || (total_cluster_num <= 0) || (global_npc_input <= 1) || isFALSE(global_seurat_obj)) {
      msg <- "Error: This error message should not appear! Please report this issue to the GenomeProtSC team."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    # TODO: Properly support Harmony clusters

    p_adj.thresh <- input$p_adj.thresh_all
    logfc.threshold <- input$logfc.threshold_all
    min.pct <- input$min.pct_all
    min.diff.pct <- input$min.diff.pct_all
    only.pos <- input$only.pos_all

    if (!(is.finite(p_adj.thresh) && (p_adj.thresh > 0) && (p_adj.thresh <= 1))) {
      msg <- "Error: p_adj.thresh must be a positive number no larger than 1."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    if (!(is.finite(logfc.threshold) && (logfc.threshold > 0))) {
      msg <- "Error: logfc.threshold must be a positive number."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    if (!(is.finite(min.pct) && (min.pct > 0) && (min.pct <= 1))) {
      msg <- "Error: min.pct must be a positive number no larger than 1."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    if (!(is.finite(min.diff.pct) && (((min.diff.pct > 0) && (min.diff.pct <= 1)) || (min.diff.pct < 0)))) {
      msg <- "Error: min.diff.pct must either be a positive number no larger than 1 or a negative number."
      session$sendCustomMessage("showStatusMessage", list(message = msg, container = "seurat-status-msg-container", color = "red"))
      return()
    }

    if (min.diff.pct < 0) {
      min.diff.pct = -Inf
    }

    file_available_seurat(FALSE)
    session$sendCustomMessage("clearStatusMessage", list(container = "seurat-status-msg-container"))
    session$sendCustomMessage("disableButtonOnly", list(id = "view_stats_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "view_qc_stats_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "perform_integration_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "rename_clusters_button"))
    session$sendCustomMessage("disableButtonOnly", list(id = "find_marker_genes_button"))
    session$sendCustomMessage("disableButton", list(id = "find_all_marker_genes_button", spinnerId = "seurat-loading-container"))

    output$marker_gene_plot <- NULL
    hide_marker_gene_outputs_seurat_module()

    markers <- FindAllMarkers(global_seurat_obj, random.seed = 1, logfc.threshold = logfc.threshold,
                              min.pct = min.pct, min.diff.pct = min.diff.pct, only.pos = only.pos)

    # Apply the adjusted p-value threshold
    markers <- markers %>% dplyr::filter(p_val_adj < p_adj.thresh)

    markers %>%
      group_by(cluster) %>%
      slice_head(n = 20) %>%
      ungroup() -> top20
    top20 <- top20$gene

    marker_gene_plot <- get_marker_gene_plot(global_seurat_obj, top20)
    output$marker_gene_plot <- renderPlot({
      suppressWarnings(print(marker_gene_plot))
    })
    show_marker_gene_outputs_seurat_module()

    ### Update the zip file

    top_level_dir <- getwd()
    output_files_dir <- outdir_seurat
    setwd(output_files_dir)

    marker_genes_file <- "marker_genes.csv"
    write.csv(markers, marker_genes_file, row.names = FALSE)

    marker_gene_plot_file <- write_marker_gene_plot(global_seurat_obj, top20)

    files_to_zip <- c(marker_genes_file, marker_gene_plot_file)
    zipfile_path <- file.path("..", seurat_file_to_download())

    zip(zipfile = zipfile_path, files = files_to_zip)
    setwd(top_level_dir)

    session$sendCustomMessage("enableButtonOnly", list(id = "view_stats_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "view_qc_stats_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "perform_integration_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "rename_clusters_button"))
    session$sendCustomMessage("enableButtonOnly", list(id = "find_marker_genes_button"))
    session$sendCustomMessage("enableButton", list(id = "find_all_marker_genes_button", spinnerId = "seurat-loading-container"))

    if (file.exists(file.path(session_id, seurat_file_to_download()))) {
      file_available_seurat(TRUE)
    }
  })

  # END SEURAT FILTERING AND INTEGRATION MODULE

  # FRAGPIPE MODULE

  # create reactive value for the VSN peptide counts file
  file_available_vsn <- reactiveVal(FALSE)

  # perform peptide counts normalization when submit is pressed
  observeEvent(input$peptide_counts_submit_button, {
    file_available_vsn(FALSE)

    # ensure download button remains greyed out (if submit is re-pressed)
    # shinyjs::disable("vsn_download_button")
    # shinyjs::runjs("document.getElementById('vsn_download_button').style.backgroundColor = '#d3d3d3';")
    # disable submit button after it is pressed
    # session$sendCustomMessage("disableButton", list(id = "vsn_download_button", spinnerId = "vsn-loading-container"))

    error_msg <- normalize_peptide_counts(input, session_id)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "vsn-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "vsn-status-msg-container"))

    # check if the zip file is created
    if (file.exists(file.path(session_id, "vsn_peptide_counts.csv"))) {
      file_available_vsn(TRUE)
    }
  })

  # END FRAGPIPE MODULE

  # INTEGRATION MODULE

  file_available_integ <- reactiveVal(FALSE)
  file_available_integ_marker <- reactiveVal(FALSE)

  observeEvent(input$integ_submit_button, { 
    message("Validating integration module input...")
    error_msg <- is_integ_input_valid(input)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "integ-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "integ-status-msg-container"))
    message("Validated!")

    session$sendCustomMessage("disableButton", list(id = "integ_submit_button", spinnerId = "integ-loading-container")) # disable submit button

    # run integration server
    integration_server(input, output, session)

    # check if the zip file is created
    if (file.exists(file.path(session_id, "integration_results.zip"))) {
      file_available_integ(TRUE)
    }
  })

  observeEvent(input$integ_marker_submit_button, {
    message("Validating integration module input (for marker genes)...")
    error_msg <- is_integ_marker_input_valid(input)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "integ-marker-status-msg-container", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "integ-marker-status-msg-container"))
    message("Validated!")

    session$sendCustomMessage("disableButton", list(id = "integ_marker_submit_button", spinnerId = "integ-marker-loading-container")) # disable submit button

    # run integration server for marker genes
    integration_marker_server(input, session)

    # check if the zip file is created
    if (file.exists(file.path(session_id, "updated_combined_annotations.zip"))) {
      file_available_integ_marker(TRUE)
      shinyjs::enable("integ_marker_download_button")
      shinyjs::runjs("document.getElementById('integ_marker_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "integ_marker_submit_button", spinnerId = "integ-marker-loading-container")) # re-enable submit button
    }
  })

  observe({
    if (file_available_integ()) {
      shinyjs::enable("integ_download_button")
      shinyjs::runjs("document.getElementById('integ_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "integ_submit_button", spinnerId = "integ-loading-container")) # re-enable submit button
    }
  })

  output$integ_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_integration_results.zip")
    },
    content = function(file) {
      file.copy(file.path(session_id, "integration_results.zip"), file)
    }
  )

  observe({
    if (file_available_integ_marker()) {
      shinyjs::enable("integ_marker_download_button")
      shinyjs::runjs("document.getElementById('integ_marker_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "integ_marker_submit_button", spinnerId = "integ-marker-loading-container")) # re-enable submit button
    }
  })

  output$integ_marker_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_updated_combined_annotations.zip")
    },
    content = function(file) {
      file.copy(file.path(session_id, "updated_combined_annotations.zip"), file)
    }
  )

  # END INTEGRATION MODULE

  # VISUALIZATION MODULE

  vis_seurat_obj <- FALSE
  transcript_expression_info <- c()

  update_transcript_selector_isovis <- function(session, transcripts_available) {
    transcripts_available <- sort(transcripts_available)
    updateSelectInput(session, "transcript_selector_isovis", choices = transcripts_available)
  }

  # enable download once files are available
  observe({
    if (file_available_vsn()) {
      shinyjs::enable("vsn_download_button")
      shinyjs::runjs("document.getElementById('vsn_download_button').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "vsn_download_button", spinnerId = "vsn-loading-container")) # re-enable submit button
    }
  })

  # download handler for the VSN peptide counts file
  output$vsn_download_button <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_vsn_peptide_counts.csv")
    },
    content = function(file) {
      file.copy(file.path(session_id, "vsn_peptide_counts.csv"), file)
    }
  )

  # Hide the transcript selector and plot outputs by default
  hide("transcript_selector_isovis")
  hide("transcript_expression_info_heading_isovis")
  hide("transcript_expression_info_container_isovis")
  hide("marker_gene_info_heading_isovis")
  hide("marker_gene_info_container_isovis")
  hide("umap_vis_plot_heading_isovis")
  hide("umap_vis_plot_isovis")
  hide("gene_expression_plot_heading_isovis")
  hide("gene_expression_plot_isovis")
  hide("gene_expression_dot_plot_heading_isovis")
  hide("gene_expression_dot_plot_isovis")
  hide("transcript_expression_plot_heading_isovis")
  hide("transcript_expression_plot_isovis")
  hide("transcript_expression_dot_plot_heading_isovis")
  hide("transcript_expression_dot_plot_isovis")

  # set output dir
  outdir_vis_isovis <- file.path(session_id, "visualization_isovis")

  # create output dir
  dir.create(outdir_vis_isovis, recursive = TRUE)

  # create reactive value for the new visualization module's output file
  file_available_vis_isovis <- reactiveVal(FALSE)

  # enable download once files are available
  observe({
    if (file_available_vis_isovis()) {
      shinyjs::enable("vis_download_button_isovis")
      shinyjs::runjs("document.getElementById('vis_download_button_isovis').style.backgroundColor = '#4CAF50';")
      session$sendCustomMessage("enableButton", list(id = "vis_download_button_isovis", spinnerId = "vis-loading-container_isovis")) # re-enable submit button
    }
  })

  # download handler for the visualization_results.zip file
  output$vis_download_button_isovis <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_", format(Sys.time(), "%H%M"), "_visualization_results.zip")
    },
    content = function(file) {
      file.copy(file.path(outdir_vis_isovis, "visualization_results_isovis.zip"), file)
    }
  )

  observeEvent(input$transcript_selector_isovis, {
    req(input$selected_isovis_gene)

    top_level_dir <- getwd()
    setwd(outdir_vis_isovis)

    files_to_zip_vis <- c()

    if (file.exists("UMAP.pdf")) {
      file.remove("UMAP.pdf")
    }

    if (file.exists("UMAP_Harmony.pdf")) {
      file.remove("UMAP_Harmony.pdf")
    }

    gene_to_plot <- input$selected_isovis_gene

    if (!isFALSE(vis_seurat_obj)) {
      hide("transcript_selector_isovis")
      output$transcript_expression_info_container_isovis <- renderTable(data.frame(c()), rownames = FALSE, colnames = FALSE)

      transcript_features <- FALSE
      if ("assays" %in% slotNames(vis_seurat_obj)) {
        assays <- vis_seurat_obj@assays
        if (!is.null(assays[["iso"]])) {
          iso_assay <- assays$iso
          if ("features" %in% slotNames(iso_assay)) {
            transcript_features <- rownames(iso_assay@features)

            # Only consider transcripts where the number of dashes is equal to the number of dashes in the gene name plus one
            dash_counts <- lengths(regmatches(gene_to_plot, gregexpr('-', gene_to_plot))) + 1
            same_dash_counts_indices <- which(lengths(regmatches(transcript_features, gregexpr('-', transcript_features))) == dash_counts)
            transcript_features <- transcript_features[same_dash_counts_indices]

            transcript_features <- grep(paste0("(^|-|\\b)", gene_to_plot, "($|\\b)"), transcript_features, value = TRUE)
          }
        }
      }

      if ("meta.data" %in% slotNames(vis_seurat_obj)) {
        vis_seurat_obj_metadata <- vis_seurat_obj@meta.data
        if ("harm_cluster" %in% colnames(vis_seurat_obj_metadata)) {
          umap_plot <- get_umap_harm_plot(vis_seurat_obj)
          output$umap_vis_plot_isovis <- renderPlot({
            suppressWarnings(print(umap_plot))
          })
          showElement("umap_vis_plot_heading_isovis")
          showElement("umap_vis_plot_isovis")

          umap_file <- write_umap_harm_plot(vis_seurat_obj)
          files_to_zip_vis <- c(files_to_zip_vis, umap_file)
        } else if ("seurat_clusters" %in% colnames(vis_seurat_obj_metadata)) {
          umap_plot <- get_umap_plot(vis_seurat_obj)
          output$umap_vis_plot_isovis <- renderPlot({
            suppressWarnings(print(umap_plot))
          })
          showElement("umap_vis_plot_heading_isovis")
          showElement("umap_vis_plot_isovis")

          umap_file <- write_umap_plot(vis_seurat_obj)
          files_to_zip_vis <- c(files_to_zip_vis, umap_file)
        }
      }

      if (gene_to_plot %in% Features(vis_seurat_obj)) {
        gene_expression_plot <- get_feature_plot(vis_seurat_obj, gene_to_plot)
        output$gene_expression_plot_isovis <- renderPlot({
          suppressWarnings(print(gene_expression_plot))
        })

        gene_expression_dot_plot <- get_dot_plot(vis_seurat_obj, gene_to_plot)
        output$gene_expression_dot_plot_isovis <- renderPlot({
          suppressWarnings(print(gene_expression_dot_plot))
        })

        showElement("gene_expression_plot_heading_isovis")
        showElement("gene_expression_plot_isovis")
        showElement("gene_expression_dot_plot_heading_isovis")
        showElement("gene_expression_dot_plot_isovis")

        gene_expression_plot_file <- write_gene_expression_plot(vis_seurat_obj, input$selected_isovis_gene)
        gene_expression_dot_plot_file <- write_gene_expression_dot_plot(vis_seurat_obj, input$selected_isovis_gene)
        files_to_zip_vis <- c(files_to_zip_vis, gene_expression_plot_file, gene_expression_dot_plot_file)
      }

      transcript_feature <- sub("-.*", "", input$transcript_selector_isovis)
      if (!is.null(transcript_expression_info)) {
        if (!(transcript_feature %in% transcript_expression_info$transcript_id)) {
          transcript_expression_df <- data.frame(c(transcript_feature, "Lack of evidence", "None"))
          rownames(transcript_expression_df) <- c("Selected transcript", "Level of evidence suggesting translation", "Cluster(s) the transcript is expressed in")

          output$transcript_expression_info_container_isovis <- renderTable(transcript_expression_df, rownames = TRUE, colnames = FALSE)
          showElement("transcript_expression_info_heading_isovis")
          showElement("transcript_expression_info_container_isovis")
        } else {
          transcript_feature_index <- which(transcript_expression_info$transcript_id == transcript_feature)

          level_of_evidence <- transcript_expression_info$translation_evidence_level[transcript_feature_index]
          if (level_of_evidence == 2) {
            level_of_evidence <- "Definite evidence"
          } else if (level_of_evidence == 1) {
            level_of_evidence <- "Potential evidence"
          } else {
            level_of_evidence <- "Lack of evidence"
          }

          cluster_string <- transcript_expression_info$expressed_in_clusters[transcript_feature_index]
          if (cluster_string == "") {
            cluster_string <- "None"
          }

          transcript_expression_df <- data.frame(c(transcript_feature, level_of_evidence, cluster_string))
          rownames(transcript_expression_df) <- c("Selected transcript", "Level of evidence suggesting translation", "Cluster(s) the transcript is expressed in")

          output$transcript_expression_info_container_isovis <- renderTable(transcript_expression_df, rownames = TRUE, colnames = FALSE)
          showElement("transcript_expression_info_heading_isovis")
          showElement("transcript_expression_info_container_isovis")
        }
      } else {
        hide("transcript_expression_info_heading_isovis")
        hide("transcript_expression_info_container_isovis")
        output$transcript_expression_info_container_isovis <- renderTable(data.frame(c()), rownames = FALSE, colnames = FALSE)
      }

      transcript_expression_plot <- get_feature_plot(vis_seurat_obj, input$transcript_selector_isovis)
      output$transcript_expression_plot_isovis <- renderPlot({
        suppressWarnings(print(transcript_expression_plot))
      })

      transcript_expression_dot_plot <- get_dot_plot(vis_seurat_obj, input$transcript_selector_isovis, assay = "iso")
      output$transcript_expression_dot_plot_isovis <- renderPlot({
        suppressWarnings(print(transcript_expression_dot_plot))
      })

      showElement("transcript_selector_isovis")
      showElement("transcript_expression_plot_heading_isovis")
      showElement("transcript_expression_plot_isovis")
      showElement("transcript_expression_dot_plot_heading_isovis")
      showElement("transcript_expression_dot_plot_isovis")

      transcript_expression_plot_file <- write_transcript_expression_plot(vis_seurat_obj, input$transcript_selector_isovis)
      transcript_expression_dot_plot_file <- write_transcript_expression_dot_plot(vis_seurat_obj, input$transcript_selector_isovis)
      files_to_zip_vis <- c(files_to_zip_vis, transcript_expression_plot_file, transcript_expression_dot_plot_file)
    }

    zipfile_path_vis <- "visualization_results_isovis.zip"
    if (file.exists(zipfile_path_vis)) {
      file.remove(zipfile_path_vis)
    }

    if (length(files_to_zip_vis) > 0) {
      zip(zipfile = zipfile_path_vis, files = files_to_zip_vis)
    }
    setwd(top_level_dir)

    if (file.exists(file.path(outdir_vis_isovis, zipfile_path_vis))) {
      file_available_vis_isovis(TRUE)
    }
  })

  update_marker_gene_pane <- function() {
    req(input$marker_gene_info)

    marker_gene_info <- input$marker_gene_info

    cluster <- marker_gene_info$cluster
    against <- marker_gene_info$against
    p_val_adj <- marker_gene_info$p_val_adj
    avg_log2fc <- marker_gene_info$avg_log2fc
    pct1 <- marker_gene_info$pct1
    pct2 <- marker_gene_info$pct2

    if (is.null(p_val_adj) || is.null(avg_log2fc) || is.null(pct1) || is.null(pct2) || is.null(cluster)) {
      hide("marker_gene_info_heading_isovis")
      hide("marker_gene_info_container_isovis")
      output$marker_gene_info_container_isovis <- renderTable(data.frame(c()), rownames = FALSE, colnames = FALSE)
      return()
    }

    p_val_adj <- unlist(p_val_adj)
    avg_log2fc <- unlist(avg_log2fc)
    pct1 <- unlist(pct1)
    pct2 <- unlist(pct2)
    cluster <- unlist(cluster)
    if (!is.null(against)) {
      against <- unlist(against)
    }

    num_clusters <- length(p_val_adj)
    if ((length(avg_log2fc) != num_clusters) || (length(pct1) != num_clusters) || (length(pct2) != num_clusters)) {
      hide("marker_gene_info_heading_isovis")
      hide("marker_gene_info_container_isovis")
      output$marker_gene_info_container_isovis <- renderTable(data.frame(c()), rownames = FALSE, colnames = FALSE)
      return()
    }
    if ((num_clusters == 0) || ((num_clusters > 1) && ((length(cluster) != num_clusters) || !is.null(against)))) {
      hide("marker_gene_info_heading_isovis")
      hide("marker_gene_info_container_isovis")
      output$marker_gene_info_container_isovis <- renderTable(data.frame(c()), rownames = FALSE, colnames = FALSE)
      return()
    }

    included_rownames <- c()
    content <- c()

    if (num_clusters == 1) {
      cluster <- paste(cluster, collapse = ", ")
      included_rownames <- c(included_rownames, "Cluster(s) containing marker gene (1st group)")
      content <- c(content, cluster)

      included_rownames <- c(included_rownames, "Cluster(s) used for comparison (2nd group)")
      if (!is.null(against)) {
        against <- paste(against, collapse = ", ")
        content <- c(content, against)
      } else {
        content <- c(content, "All other clusters")
      }

      included_rownames <- c(included_rownames, "Adjusted p-value")
      content <- c(content, format(p_val_adj[1], digits = 5, scientific = (p_val_adj[1] <= 1e-5)))

      included_rownames <- c(included_rownames, "log2(fold change of average expression between the two groups)")
      content <- c(content, format(avg_log2fc[1], digits = 5, scientific = (avg_log2fc[1] <= 1e-5)))

      included_rownames <- c(included_rownames, "Percentage of cells where the gene is detected in the 1st group")
      content <- c(content, pct1[1] * 100)

      included_rownames <- c(included_rownames, "Percentage of cells where the gene is detected in the 2nd group")
      content <- c(content, pct2[1] * 100)
    } else {
      for (i in 1:num_clusters) {
        cluster_name <- cluster[i]
        p_val_adj_val <- p_val_adj[i]
        avg_log2fc_val <- avg_log2fc[i]
        pct1_val <- pct1[i]
        pct2_val <- pct2[i]

        included_rownames <- c(included_rownames, paste0("Adjusted p-value (", cluster_name, ')'))
        content <- c(content, format(p_val_adj_val, digits = 5, scientific = (p_val_adj_val <= 1e-5)))

        included_rownames <- c(included_rownames, paste0("log2(fold change of average expression between ", cluster_name, " and all other clusters)"))
        content <- c(content, format(avg_log2fc_val, digits = 5, scientific = (avg_log2fc_val <= 1e-5)))

        included_rownames <- c(included_rownames, paste0("Percentage of cells where the gene is detected in ", cluster_name))
        content <- c(content, pct1_val * 100)

        included_rownames <- c(included_rownames, paste0("Percentage of cells where the gene is detected in all clusters except ", cluster_name))
        content <- c(content, pct2_val * 100)
      }
    }

    if (!is.null(input$selected_isovis_gene)) {
      included_rownames <- c("Gene", included_rownames)
      content <- c(input$selected_isovis_gene, content)
    }

    marker_gene_df <- data.frame(content)
    rownames(marker_gene_df) <- included_rownames
    output$marker_gene_info_container_isovis <- renderTable(marker_gene_df, rownames = TRUE, colnames = FALSE)
    showElement("marker_gene_info_heading_isovis")
    showElement("marker_gene_info_container_isovis")
  }

  # Handling user-uploaded files
  observeEvent(input$vis_submit_button_isovis, {
    error_msg <- is_vis_isovis_input_valid(input)
    if (error_msg != "") {
      message(error_msg)
      session$sendCustomMessage("showStatusMessage", list(message = error_msg, container = "vis-status-msg-container_isovis", color = "red"))
      return()
    }
    session$sendCustomMessage("clearStatusMessage", list(container = "vis-status-msg-container_isovis"))

    top_level_dir <- getwd()
    setwd(outdir_vis_isovis)

    files_to_zip_vis <- c()

    if (file.exists("UMAP.pdf")) {
      file.remove("UMAP.pdf")
    }

    if (file.exists("UMAP_Harmony.pdf")) {
      file.remove("UMAP_Harmony.pdf")
    }

    seurat_obj <- readRDS(input$seurat_rds_file_isovis$datapath)
    if ("meta.data" %in% slotNames(seurat_obj)) {
      seurat_obj_metadata <- seurat_obj@meta.data
      vis_seurat_obj <<- seurat_obj
      if ("harm_cluster" %in% colnames(seurat_obj_metadata)) {
        umap_plot <- get_umap_harm_plot(seurat_obj)
        output$umap_vis_plot_isovis <- renderPlot({
          suppressWarnings(print(umap_plot))
        })
        showElement("umap_vis_plot_heading_isovis")
        showElement("umap_vis_plot_isovis")

        umap_file <- write_umap_harm_plot(vis_seurat_obj)
        files_to_zip_vis <- c(files_to_zip_vis, umap_file)
      } else if ("seurat_clusters" %in% colnames(seurat_obj_metadata)) {
        umap_plot <- get_umap_plot(seurat_obj)
        output$umap_vis_plot_isovis <- renderPlot({
          suppressWarnings(print(umap_plot))
        })
        showElement("umap_vis_plot_heading_isovis")
        showElement("umap_vis_plot_isovis")

        umap_file <- write_umap_plot(vis_seurat_obj)
        files_to_zip_vis <- c(files_to_zip_vis, umap_file)
      }
    }

    if (!is.null(input$transcript_info_file_isovis)) {
      transcript_expression_info <<- read.csv(input$transcript_info_file_isovis$datapath)
      output$transcript_expression_info_container_isovis <- renderTable(data.frame(c()), rownames = FALSE, colnames = FALSE)
    }

    zipfile_path_vis <- "visualization_results_isovis.zip"
    if (file.exists(zipfile_path_vis)) {
      file.remove(zipfile_path_vis)
    }

    if (length(files_to_zip_vis) > 0) {
      zip(zipfile = zipfile_path_vis, files = files_to_zip_vis)
    }
    setwd(top_level_dir)

    if (file.exists(file.path(outdir_vis_isovis, zipfile_path_vis))) {
      file_available_vis_isovis(TRUE)
    }
  })

  observeEvent(input$selected_isovis_gene, {
    req(vis_seurat_obj)

    gene_to_plot <- input$selected_isovis_gene
    if (nchar(gene_to_plot) == 0) {
      hide("transcript_selector_isovis")
      hide("transcript_expression_info_heading_isovis")
      hide("transcript_expression_info_container_isovis")
      hide("marker_gene_info_heading_isovis")
      hide("marker_gene_info_container_isovis")
      hide("umap_vis_plot_heading_isovis")
      hide("umap_vis_plot_isovis")
      hide("gene_expression_plot_heading_isovis")
      hide("gene_expression_plot_isovis")
      hide("gene_expression_dot_plot_heading_isovis")
      hide("gene_expression_dot_plot_isovis")
      hide("transcript_expression_plot_heading_isovis")
      hide("transcript_expression_plot_isovis")
      hide("transcript_expression_dot_plot_heading_isovis")
      hide("transcript_expression_dot_plot_isovis")

      output$transcript_expression_info_container_isovis <- renderTable(data.frame(c()), rownames = FALSE, colnames = FALSE)
      output$marker_gene_info_container_isovis <- renderTable(data.frame(c()), rownames = FALSE, colnames = FALSE)
      output$umap_vis_plot_isovis <- renderPlot(NULL)
      output$gene_expression_plot_isovis <- renderPlot(NULL)
      output$gene_expression_dot_plot_isovis <- renderPlot(NULL)
      output$transcript_expression_plot_isovis <- renderPlot(NULL)
      output$transcript_expression_dot_plot_isovis <- renderPlot(NULL)
      return()
    }

    top_level_dir <- getwd()
    setwd(outdir_vis_isovis)

    files_to_zip_vis <- c()

    if (file.exists("UMAP.pdf")) {
      file.remove("UMAP.pdf")
    }

    if (file.exists("UMAP_Harmony.pdf")) {
      file.remove("UMAP_Harmony.pdf")
    }

    transcript_features <- FALSE
    if ("assays" %in% slotNames(vis_seurat_obj)) {
      assays <- vis_seurat_obj@assays
      if (!is.null(assays[["iso"]])) {
        iso_assay <- assays$iso
        if ("features" %in% slotNames(iso_assay)) {
          hide("transcript_selector_isovis")
          transcript_features <- rownames(iso_assay@features)

          # Only consider transcripts where the number of dashes is equal to the number of dashes in the gene name plus one
          dash_counts <- lengths(regmatches(gene_to_plot, gregexpr('-', gene_to_plot))) + 1
          same_dash_counts_indices <- which(lengths(regmatches(transcript_features, gregexpr('-', transcript_features))) == dash_counts)
          transcript_features <- transcript_features[same_dash_counts_indices]

          transcript_features <- grep(paste0("(^|-|\\b)", gene_to_plot, "($|\\b)"), transcript_features, value = TRUE)
        }
      }
    }

    if ("meta.data" %in% slotNames(vis_seurat_obj)) {
      vis_seurat_obj_metadata <- vis_seurat_obj@meta.data
      if ("harm_cluster" %in% colnames(vis_seurat_obj_metadata)) {
        umap_plot <- get_umap_harm_plot(vis_seurat_obj)
        output$umap_vis_plot_isovis <- renderPlot({
          suppressWarnings(print(umap_plot))
        })
        showElement("umap_vis_plot_heading_isovis")
        showElement("umap_vis_plot_isovis")

        umap_file <- write_umap_harm_plot(vis_seurat_obj)
        files_to_zip_vis <- c(files_to_zip_vis, umap_file)
      } else if ("seurat_clusters" %in% colnames(vis_seurat_obj_metadata)) {
        umap_plot <- get_umap_plot(vis_seurat_obj)
        output$umap_vis_plot_isovis <- renderPlot({
          suppressWarnings(print(umap_plot))
        })
        showElement("umap_vis_plot_heading_isovis")
        showElement("umap_vis_plot_isovis")

        umap_file <- write_umap_plot(vis_seurat_obj)
        files_to_zip_vis <- c(files_to_zip_vis, umap_file)
      }
    }

    if (input$selected_isovis_gene %in% Features(vis_seurat_obj)) {
      gene_expression_plot <- get_feature_plot(vis_seurat_obj, input$selected_isovis_gene)
      output$gene_expression_plot_isovis <- renderPlot({
        suppressWarnings(print(gene_expression_plot))
      })

      gene_expression_dot_plot <- get_dot_plot(vis_seurat_obj, input$selected_isovis_gene)
      output$gene_expression_dot_plot_isovis <- renderPlot({
        suppressWarnings(print(gene_expression_dot_plot))
      })

      showElement("gene_expression_plot_heading_isovis")
      showElement("gene_expression_plot_isovis")
      showElement("gene_expression_dot_plot_heading_isovis")
      showElement("gene_expression_dot_plot_isovis")

      gene_expression_plot_file <- write_gene_expression_plot(vis_seurat_obj, input$selected_isovis_gene)
      gene_expression_dot_plot_file <- write_gene_expression_dot_plot(vis_seurat_obj, input$selected_isovis_gene)
      files_to_zip_vis <- c(files_to_zip_vis, gene_expression_plot_file, gene_expression_dot_plot_file)
    }

    if (!isFALSE(transcript_features) && (length(transcript_features) != 0)) {
      update_transcript_selector_isovis(session, transcript_features)
      showElement("transcript_selector_isovis")
    }

    update_marker_gene_pane()

    zipfile_path_vis <- "visualization_results_isovis.zip"
    if (file.exists(zipfile_path_vis)) {
      file.remove(zipfile_path_vis)
    }

    zip(zipfile = zipfile_path_vis, files = files_to_zip_vis)
    setwd(top_level_dir)

    if (file.exists(file.path(outdir_vis_isovis, zipfile_path_vis))) {
      file_available_vis_isovis(TRUE)
    }
  })

  # END VISUALIZATION MODULE

  # remove session id tmp directory created each time app is run
  session$onSessionEnded(function() {
    if (dir.exists(session_id)) {
      unlink(session_id, recursive = TRUE)
    }
  })
}