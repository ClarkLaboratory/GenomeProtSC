# inputs: MQ/FragPipe peptides.tsv output, GTF, metadata, FASTA, Seurat object RDS file
# outputs: BED12/GTF files of peptides, ORFs and transcripts, database of peptides with info on locations etc, summary file of peptides

library(stringi)
source("global.R")
source("R/integration_functions.R")

option_list = list(
  make_option(c("-p", "--proteomics"), type="character", default=NULL,
              help="Proteomics data file", metavar="character"),
  make_option(c("-f", "--fasta"), type="character", default=NULL,
              help="Custom FASTA used for proteomics", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Custom metadata used for proteomics", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="GTF used to generate custom FASTA", metavar="character"),
  make_option(c("-r", "--rds"), type="character", default=NULL,
              help="Seurat object file", metavar="character"),
  make_option(c("-s", "--savepath"), type="character", default=NULL,
              help="Output directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

options(scipen=999)

proteomics_import_file <- opt$proteomics
fasta_import_file <- opt$fasta
metadata_import_file <- opt$metadata
gtf_import_file <- opt$gtf
rds_import_file <- opt$rds
output_directory <- opt$savepath

if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

# ------------- import files ------------- #

pd <- suppressWarnings(import_proteomics_data(proteomics_import_file))

gtf <- makeTxDbFromGFF(gtf_import_file, format="gtf") # make txdb of gtf

orf_df <- import_orf_metadata(metadata_import_file)

seurat_obj <- readRDS(rds_import_file)
all_transcript_id_gene_pairs <- Features(seurat_obj, assay = "iso")
all_transcript_ids <- sub("-.*", "", all_transcript_id_gene_pairs)

gtf_for_exporting <- import(gtf_import_file, format="gtf")

md <- import_fasta(fasta_import_file, pd, gtf, gtf_for_exporting)

# ------------- run analysis ------------- #

# ------- map ORF transcript coords to spliced genomic coords ------- #

# create unique ID of ORF in transcript
md$orf_tx_id <- paste0(md$protein_name, "_", md$transcript_id)

# filter for unique orf and transcript for mapping the coordinates
orf_transcript_coords_df <- md %>% dplyr::select(orf_tx_id, txstart, txend, transcript_id, gene_id, strand)
orf_transcript_coords_df <- orf_transcript_coords_df[!(base::duplicated(orf_transcript_coords_df)),] # remove duplicates

# make GRanges from df of ORF transcript coordinates
orf_transcript_coords <- makeGRangesFromDataFrame(orf_transcript_coords_df,
                                                  keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
                                                  seqnames.field="transcript_id", start.field="txstart", end.field="txend", strand.field="strand",
                                                  starts.in.df.are.0based=FALSE, na.rm=TRUE)

names(orf_transcript_coords) <- c(orf_transcript_coords$orf_tx_id) # set names as unique IDs
mcols(orf_transcript_coords)$gene_id <- c(orf_transcript_coords_df$gene_id) # set gene_id

# get exons for mapping coordinates to genome
exons <- exonsBy(gtf, "tx", use.names=TRUE) # get exon data per transcript
exons_filt <- exons[names(exons) %in% orf_transcript_coords_df$transcript_id] # filter for only transcripts with mapped peptides

orf_tx_names <- as.character(seqnames(orf_transcript_coords)) # get tx names

# match names of transcripts, return index of match
names(orf_transcript_coords) <- match(orf_tx_names, names(exons_filt)) 

# create vector of unique names and gene ID for later
orf_ids <- orf_transcript_coords$orf_tx_id
orf_gene_ids <- orf_transcript_coords$gene_id

# use ORFik to map transcript coordinates to spliced genomic coordinates
orf_in_genomic <- ORFik::pmapFromTranscriptF(orf_transcript_coords, exons_filt, removeEmpty = T)

# map back to GRangesList, with group information, add gene_id back
orf_in_genomic@unlistData$PID <- orf_ids[groupings(orf_in_genomic)]
orf_in_genomic@unlistData$gene_id <- orf_gene_ids[groupings(orf_in_genomic)]

# unlist to add exon_number for GTF export
orf_in_genomic_gr <- unlist(orf_in_genomic, use.names=F) # convert to GRanges

# create vector of exon number per peptide and transcript
exon_number_vec <- ave(seq_along(orf_in_genomic_gr), mcols(orf_in_genomic_gr)$PID, FUN = seq_along)

# add to GRanges
mcols(orf_in_genomic_gr)$exon_number <- exon_number_vec
# re-list
orf_in_genomic <- split(orf_in_genomic_gr, ~ mcols(orf_in_genomic_gr)$PID)

# ------- map peptide transcript coords to spliced genomic coords ------- #

# use ORF transcript coords to determine peptide transcript coords
peptide_transcript_coords <- extract_peptide_coords(md)

# get vector of transcript names
peptide_tx_names <- as.character(seqnames(peptide_transcript_coords)) # get tx names

# match names of transcripts, return index of match
names(peptide_transcript_coords) <- match(peptide_tx_names, names(exons_filt)) 

# create vectors of unique names and gene ID for later
pep_ids <- peptide_transcript_coords$peptide
pep_PID_ids <- peptide_transcript_coords$PID
pep_gene_ids <- peptide_transcript_coords$gene_id

# use ORFik to map transcript coordinates to spliced genomic coordinates
pep_in_genomic <- ORFik::pmapFromTranscriptF(peptide_transcript_coords, exons_filt, removeEmpty = F)

# map back to GRangesList, with group information, add other info back
pep_in_genomic@unlistData$peptide <- pep_ids[groupings(pep_in_genomic)]
pep_in_genomic@unlistData$PID <- pep_PID_ids[groupings(pep_in_genomic)]
pep_in_genomic@unlistData$gene <- pep_gene_ids[groupings(pep_in_genomic)]

# unlist to add exon_number for GTF export
pep_in_genomic_gr <- unlist(pep_in_genomic, use.names=F) # convert to GRanges

# rename with transcript and peptide
tx_pep_names <- c(paste0(names(pep_in_genomic_gr), "_", pep_in_genomic_gr$peptide))
names(pep_in_genomic_gr) <- tx_pep_names # set names
pep_in_genomic_gr$tx_pid_grouping <- paste0(pep_in_genomic_gr$PID, "_", names(pep_in_genomic_gr))

# remove 0 ranges 
pep_in_genomic_gr <- subset(pep_in_genomic_gr, (start(pep_in_genomic_gr) != 0 & end(pep_in_genomic_gr) != 0)  )

# create vector of exon number per peptide and transcript
exon_number_vec <- ave(seq_along(pep_in_genomic_gr), pep_in_genomic_gr$tx_pid_grouping, FUN = seq_along)

# add to GRanges
mcols(pep_in_genomic_gr)$exon_number <- exon_number_vec
mcols(pep_in_genomic_gr)$tx_pid_grouping <- NULL
# re-list
pep_in_genomic <- split(pep_in_genomic_gr, ~ names(pep_in_genomic_gr))

# ------------- export BED12 files ------------- #

# export bed12 of peptides
ORFik::export.bed12(pep_in_genomic, file.path(output_directory, "peptides.bed12"), rgb = 0)

# export bed12 of ORFs
# currently export with PID_transcript as the name, means ORFs are often present multiple times
ORFik::export.bed12(orf_in_genomic, file.path(output_directory, "ORFs.bed12"), rgb = 0)

# format GTF of all transcripts that had mapped peptides
# gtf_for_exporting <- import(gtf_import_file, format="gtf")
gtf_filtered <- gtf_for_exporting[mcols(gtf_for_exporting)$transcript_id %in% md$transcript]
gtf_filtered$group_id <- "transcripts"

# reformat exons for bed12
gtf_as_bed12 <- gtf_filtered[mcols(gtf_filtered)$type == "exon"]

names(gtf_as_bed12) <- paste0(gtf_as_bed12$transcript_id, "_", gtf_as_bed12$gene_id)

# convert to grl
tx_in_genomic <- split(gtf_as_bed12, ~ names(gtf_as_bed12))

# export bed12 of transcripts
ORFik::export.bed12(tx_in_genomic, file.path(output_directory, "transcripts.bed12"), rgb = 0)

# ------- transcripts and IsoVis GTF -------- #

# for combined GTF

# create object unformatted for IsoVis
orf_in_genomic_gr_isovis <- orf_in_genomic_gr

# format GTF to combine with transcripts and peptides last
orf_in_genomic_gr$source <- c("custom")
orf_in_genomic_gr$type <- c("CDS")
orf_in_genomic_gr$phase <- 0
orf_in_genomic_gr$ORF_id <- names(orf_in_genomic_gr)
orf_in_genomic_gr$transcript_id <- names(orf_in_genomic_gr)

names(orf_in_genomic_gr) <- NULL
orf_in_genomic_gr$group_id <- "ORFs"

# for IsoVis

# re-filter initial GTF of transcripts and exons
gtf_isovis <- gtf_for_exporting[mcols(gtf_for_exporting)$transcript_id %in% md$transcript]

# set standard GTF columns
orf_in_genomic_gr_isovis$source <- c("custom")
orf_in_genomic_gr_isovis$type <- c("CDS")
orf_in_genomic_gr_isovis$phase <- 0
orf_in_genomic_gr_isovis$transcript_id <- names(orf_in_genomic_gr)

# get PID_transcript column as df
pids <- data.frame(ids = orf_in_genomic_gr_isovis$PID)
# separate into just protein ID based on last occurrence of an '_'
pids <- pids %>% separate(ids, into = "protein_id", sep = "\\_(?!.*_)", remove=F)
# paste as protein_id
orf_in_genomic_gr_isovis$protein_id <- paste0(pids$protein_id)
# remove other columns
names(orf_in_genomic_gr_isovis) <- NULL
orf_in_genomic_gr_isovis$PID <- NULL

# combine transcripts, exons and CDS
isovis_export <- c(gtf_isovis, orf_in_genomic_gr_isovis)
# sort
isovis_export <- sortSeqlevels(isovis_export)
isovis_export <- sort(isovis_export)
# export GTF compatible with IsoVis
rtracklayer::export(isovis_export, file.path(output_directory, "transcripts_and_ORFs_for_isovis.gtf"), format="gtf")

# ------- summary file of peptide mappings -------- #

# convert to df
mcols(pep_in_genomic_gr)$txname <- names(pep_in_genomic_gr)
results_pept_df <- pep_in_genomic_gr %>% as_tibble()
results_pept_df <- separate(results_pept_df, txname, into = c("transcript_id", "peptide"), sep = "_", remove = TRUE)
results_pept_df$gene_id <- results_pept_df$gene
results_pept_df$gene <- NULL

# group by peptide and transcript to summarise based on how many exons peptide spans
results_pept_df_unique <- results_pept_df %>% 
  dplyr::group_by(peptide, transcript_id) %>% 
  dplyr::slice_max(exon_number) %>% 
  dplyr::mutate(number_exons = exon_number) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-start, -end, -width, -exon_number)

# merge results with metadata
metadata_to_merge <- md %>% 
  dplyr::select(PID, peptide, transcript_id, gene_id, gene_name, protein_length, tx_len)

peptide_result <- left_join(results_pept_df_unique, metadata_to_merge, by=c("PID", "peptide", "gene_id", "transcript_id"))

peptide_result <- peptide_result[!(base::duplicated(peptide_result)),]

peptide_result <- left_join(peptide_result, orf_df, by=c("PID", "transcript_id"))

peptide_result <- peptide_result[!(base::duplicated(peptide_result)),]

peptide_result$seqnames <- NULL

peptide_result <- peptide_result %>% 
  dplyr::mutate(longest_orf_in_transcript = case_when(
    longest_orf_in_transcript == "Y" ~ TRUE,
    longest_orf_in_transcript == "N" ~ FALSE
  ))

setDT(peptide_result)

peptide_result[, c("peptide_ids_gene", "peptide_ids_orf", "peptide_ids_transcript", "shared_novel_protein_peptide") := 
                 .(length(unique(gene_id)) == 1,
                   length(unique(gene_id)) == 1 & length(unique(PID)) == 1,
                   length(unique(gene_id)) == 1 & length(unique(PID)) == 1 & length(unique(transcript_id)) == 1,
                   length(unique(PID)) > 1 & all(startsWith(PID, "ORF"))),
               by = peptide]

peptide_result[, orf_identified := any(peptide_ids_orf == TRUE), by = PID]
peptide_result[, gene_identified := any(peptide_ids_gene == TRUE), by = gene_id]
peptide_result[, transcript_identified := any(peptide_ids_transcript == TRUE), by = transcript_id]

# get missing peptides back in the output
missing_peptides <- pd %>% dplyr::filter(!(peptide %in% peptide_result$peptide))

# ensure same col names
columns_to_add <- setdiff(names(peptide_result), names(missing_peptides))

# add missing columns to missing df with NA values
missing_peptides[columns_to_add] <- NA

# ensure the column order is same before rbind
missing_peptides <- missing_peptides[, names(peptide_result)]

combined_peptide_result <- rbind(peptide_result, missing_peptides)

combined_peptide_result$PID <- gsub(",", ".", combined_peptide_result$PID)
combined_peptide_result <- combined_peptide_result[!(base::duplicated(combined_peptide_result)),]
combined_peptide_result$transcript_length <- combined_peptide_result$tx_len

combined_peptide_result <- combined_peptide_result %>% 
  mutate(simplified_biotype = case_when(
    transcript_biotype %in% c("protein_coding", "protein_coding_LoF", "protein_coding_CDS_not_defined") ~ "protein_coding",
    transcript_biotype %in% c("polymorphic_pseudogene", "pseudogene", "processed_pseudogene", 
                              "transcribed_unprocessed_pseudogene", "transcribed_processed_pseudogene",
                              "unprocessed_pseudogene", "transcribed_unitary_pseudogene", "translated_processed_pseudogene",
                              "translated_unprocessed_pseudogene", "unitary_pseudogene") ~ "pseudogene",
    transcript_biotype %in% c("nonsense_mediated_decay", "non_stop_decay") ~ "NMD",
    transcript_biotype %in% c("retained_intron") ~ "retained_intron",
    transcript_biotype %in% c("lncRNA") ~ "lncRNA",
    transcript_biotype %in% c("novel") ~ "novel",
    TRUE ~ "other"
  ))

# rearrange columns for output
combined_peptide_result <- combined_peptide_result %>% dplyr::select(peptide,accession,PID,transcript_id,gene_id,gene_name,strand,
                                                                     number_exons,transcript_length,transcript_biotype,simplified_biotype,
                                                                     protein_length,orf_genomic_coordinates,orf_type,localisation,uniprot_status,
                                                                     openprot_id,`molecular_weight(kDA)`,isoelectric_point,hydrophobicity,
                                                                     aliphatic_index,longest_orf_in_transcript,peptide_ids_gene,peptide_ids_orf,
                                                                     peptide_ids_transcript,shared_novel_protein_peptide,orf_identified,
                                                                     gene_identified,transcript_identified)

# export summary data
write.csv(combined_peptide_result, file.path(output_directory, "peptide_info.csv"), row.names=F, quote=F)

# ------- include information of all other expressed transcripts that do not have mapped peptides ------- #

gtf_remaining <- gtf_for_exporting[!(mcols(gtf_for_exporting)$transcript_id %in% md$transcript)]
gtf_remaining <- gtf_remaining[mcols(gtf_remaining)$transcript_id %in% all_transcript_ids]
gtf_remaining$group_id <- "transcripts"

# ------- combined GTF -------- #

# include orf_status and peptide_status in GTF mcols
results_to_merge_with_granges <- left_join(results_pept_df, peptide_result, by=c("transcript_id", "peptide", "strand", "PID", "gene_id"))
results_to_merge_with_granges <- results_to_merge_with_granges[!(duplicated(results_to_merge_with_granges)),]
results_to_merge_with_granges <- results_to_merge_with_granges %>% 
  dplyr::select(-openprot_id, -`molecular_weight(kDA)`, -isoelectric_point, -hydrophobicity, -aliphatic_index)
results_to_merge_with_granges$naming <- paste0(results_to_merge_with_granges$transcript_id, "_", results_to_merge_with_granges$peptide)

# make GRanges from df of ORF transcript coordinates
pep_in_genomic_gr_export <- makeGRangesFromDataFrame(results_to_merge_with_granges,
                                                     keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
                                                     seqnames.field="seqnames", start.field="start", end.field="end", strand.field="strand",
                                                     starts.in.df.are.0based=FALSE, na.rm=TRUE)

names(pep_in_genomic_gr_export) <- c(pep_in_genomic_gr_export$naming) # set names

# add mcols
pep_in_genomic_gr_export$source <- "custom"
pep_in_genomic_gr_export$type <- "exon"
pep_in_genomic_gr_export$group_id <- "peptides"

# export annotations for vis
gtf_for_exporting$group_id <- "transcripts"
combined <- c(pep_in_genomic_gr_export, orf_in_genomic_gr, gtf_for_exporting)
combined_annotations_temp_gtf_path <- file.path(output_directory, "combined_annotations_temp.gtf")
rtracklayer::export(combined, combined_annotations_temp_gtf_path, format = "gtf")

# Prepend a comment to the combined annotations GTF file to signify that it comes from GenomeProt(SC)
combined_annotations_gtf_path <- file.path(output_directory, "combined_annotations.gtf")
system(paste0("echo '##GenomeProt' > ", combined_annotations_gtf_path))
system(paste0("cat ", combined_annotations_temp_gtf_path, " >> ", combined_annotations_gtf_path))
file.remove(combined_annotations_temp_gtf_path)

# ------- level of translational evidence TXT ------- #

# If a transcript is uniquely identified by a peptide, there is definite evidence supporting it being translated
transcript_ids_with_definite_evidence <- unique((peptide_result %>% dplyr::filter(isTRUE(peptide_ids_orf)))[["transcript_id"]]) # isTRUE()? == "TRUE"? == TRUE?

# Any transcript that has a peptide mapping has potential evidence supporting its translation
transcript_ids_with_potential_evidence <- unique(peptide_result[["transcript_id"]])

# Remove transcript IDs that are definitely translated
transcript_ids_with_potential_evidence <- setdiff(transcript_ids_with_potential_evidence, transcript_ids_with_definite_evidence)

if ("meta.data" %in% slotNames(seurat_obj)) {
  seurat_obj_metadata <- seurat_obj@meta.data
  if ("harm_cluster" %in% colnames(seurat_obj_metadata)) {
    cluster_numbers <- seurat_obj@meta.data$harm_cluster
    cluster_names <- levels(cluster_numbers)
  } else if ("seurat_clusters" %in% colnames(seurat_obj_metadata)) {
    cluster_numbers <- seurat_obj$seurat_clusters
    cluster_names <- levels(seurat_obj)
  }
}

# Cache the cell barcode indices for each cluster
indices_for_cluster_list <- list()
for (j in seq_along(cluster_names)) {
  indices_for_cluster <- which(cluster_numbers == j - 1)
  indices_for_cluster_list[[j]] <- indices_for_cluster
}

output_txt_table <- data.frame()

# For each feature (i.e. expressed transcript) found in the scRNA-seq experiment, determine the cell clusters that have expressed it
expressed_vals_info_all <- FetchData(seurat_obj[["iso"]], vars = all_transcript_id_gene_pairs, layer = "counts")
for (i in seq_along(all_transcript_id_gene_pairs)) {
  transcript_id_gene_pair <- all_transcript_id_gene_pairs[i]
  transcript_id <- all_transcript_ids[i]
  expressed_vals_info <- expressed_vals_info_all[[transcript_id_gene_pair]]
  expressed_clusters <- c()
  for (j in seq_along(cluster_names)) {
    indices_for_cluster <- indices_for_cluster_list[[j]]
    expression_vals <- expressed_vals_info[indices_for_cluster]
    if (any(expression_vals > 0)) {
      cluster_name <- cluster_names[j]
      expressed_clusters <- c(expressed_clusters, cluster_name)
    }
  }
  expressed_clusters_string <- paste(expressed_clusters, collapse = ", ")
  new_row <- c(transcript_id, expressed_clusters_string)
  output_txt_table <- rbind(output_txt_table, new_row)
}

# Add a column that represents the level of evidence supporting a transcript being translated
# 2 = Definitely translated, 1 = Potentially translated, 0 or missing = Lack of evidence
output_txt_table <- cbind(output_txt_table, rep(0, nrow(output_txt_table)))

colnames(output_txt_table) <- c("transcript_id", "expressed_in_clusters", "translation_evidence_level")

# Add information for transcripts that are found in the GTF

# 2: Do this for transcript IDs with definite translation evidence
for (i in seq_along(transcript_ids_with_definite_evidence)) {
  transcript_id <- transcript_ids_with_definite_evidence[i]
  if (transcript_id %in% output_txt_table$transcript_id) {
    j <- which(output_txt_table$transcript_id == transcript_id)[1]
    output_txt_table$translation_evidence_level[j] <- 2
  } else {
    new_row <- c(transcript_id, "", 2)
    output_txt_table <- rbind(output_txt_table, new_row)
  }
}

# 1: Do this for transcript IDs with potential translation evidence
for (i in seq_along(transcript_ids_with_potential_evidence)) {
  transcript_id <- transcript_ids_with_potential_evidence[i]
  if (transcript_id %in% output_txt_table$transcript_id) {
    j <- which(output_txt_table$transcript_id == transcript_id)[1]
    output_txt_table$translation_evidence_level[j] <- 1
  } else {
    new_row <- c(transcript_id, "", 1)
    output_txt_table <- rbind(output_txt_table, new_row)
  }
}

# Write the file
write.csv(output_txt_table, file = file.path(output_directory, "transcript_expression_info.csv"), row.names = FALSE)