library(FLAMES)
library(SingleCellExperiment)
library(optparse)

option_list = list(
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="Genome file", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL,
              help="GTF file", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Directory containing FASTQs", metavar="character"),
  make_option(c("-c", "--cell_count"), type="numeric", default=NULL,
              help="Expected number of cells per sample", metavar="numeric"),
  make_option(c("-n", "--ndr"), type="numeric", default=NULL,
              help="Bambu NDR (0.0 <= NDR <= 1.0)", metavar="numeric"),
  make_option(c("-t", "--threads"), type="numeric", default=8,
              help="Number of CPU threads (default: 8)", metavar="numeric"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

genome_file <- opt$genome
annotation_file <- opt$annotation
cell_num <- opt$cell_count
ndr <- opt$ndr
cpu_threads <- opt$threads
output <- opt$output
fastqs_dir <- opt$input

fastqs_in <- list.files(fastqs_dir, full.names = TRUE)

if (!is.numeric(ndr) || !is.finite(ndr) || (ndr < 0) || (ndr > 1)) {
  ndr <- NULL
}
message(paste0("Bambu NDR used: ", ndr))

if (!is.numeric(cpu_threads) || !is.finite(cpu_threads) || (cpu_threads <= 0) || (cpu_threads >= 24)) {
  cpu_threads <- 8
}

flames_config <- NULL

if (!is.null(ndr)) {
  flames_config <- FLAMES::create_config(
    outdir = output,
    threads = cpu_threads,
    pipeline_parameters.demultiplexer = "BLAZE",
    bambu_isoform_identification = TRUE,
    multithread_isoform_identification = TRUE,
    isoform_parameters.bambu_ndr = ndr
  )
} else {
  flames_config <- FLAMES::create_config(
    outdir = output,
    threads = cpu_threads,
    pipeline_parameters.demultiplexer = "BLAZE",
    bambu_isoform_identification = TRUE,
    multithread_isoform_identification = TRUE
  )

  # Set bambu_ndr to NULL
  temp_json <- readChar(flames_config, file.info(flames_config)$size)
  temp_json <- gsub('\\"bambu_ndr\\": \\[0.5\\]', '"bambu_ndr": null', temp_json)
  writeChar(temp_json, flames_config, eos = NULL)
}

if (length(fastqs_in) > 1) {
  message("Multiple fastq files detected")
  pipeline <- FLAMES::MultiSampleSCPipeline(
    config_file = flames_config,
    outdir = output,
    fastq = fastqs_in,
    annotation = annotation_file,
    genome_fa = genome_file,
    barcodes_file = NULL,
    expect_cell_number = cell_num
  )
} else {
  message("Single fastq file detected")
  pipeline <- FLAMES::SingleCellPipeline(
    config_file = flames_config,
    outdir = output,
    fastq = fastqs_in,
    annotation = annotation_file,
    genome_fa = genome_file,
    barcodes_file = NULL,
    expect_cell_number = cell_num
  )
}

pipeline <- FLAMES::run_FLAMES(pipeline)