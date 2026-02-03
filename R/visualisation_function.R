library(Seurat)
library(Matrix)
library(dplyr)
library(data.table)
library(tibble)
library(ggplot2)

get_stats_msg <- function(seurat_obj) {
  # Get the number of cells
  num_cells <- ncol(seurat_obj)

  # Get statistics about the number of molecules detected per cell
  nCount_RNA <- seurat_obj$nCount_RNA
  min_nCount_RNA <- min(nCount_RNA)
  max_nCount_RNA <- max(nCount_RNA)
  median_nCount_RNA <- median(nCount_RNA)
  stdev_nCount_RNA <- sd(nCount_RNA)

  # Get statistics about the number of genes detected per cell
  nFeature_RNA <- seurat_obj$nFeature_RNA
  min_nFeature_RNA <- min(nFeature_RNA)
  max_nFeature_RNA <- max(nFeature_RNA)
  median_nFeature_RNA <- median(nFeature_RNA)
  stdev_nFeature_RNA <- sd(nFeature_RNA)

  # Get statistics about mitochondrial genes expressed per cell
  percent_mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  min_percent_mt <- min(percent_mt)
  max_percent_mt <- max(percent_mt)
  median_percent_mt <- median(percent_mt)
  stdev_percent_mt <- sd(percent_mt)

  msg <- paste0("Number of cells: ", num_cells, "\n",
                "\n",
                "Number of molecules detected per cell (nCount_RNA):\n",
                "Min: ", min_nCount_RNA, "\n",
                "Max: ", max_nCount_RNA, "\n",
                "Median: ", median_nCount_RNA, "\n",
                "Stdev: ", stdev_nCount_RNA, "\n",
                "\n",
                "Number of genes detected per cell (nFeature_RNA):\n",
                "Min: ", min_nFeature_RNA, "\n",
                "Max: ", max_nFeature_RNA, "\n",
                "Median: ", median_nFeature_RNA, "\n",
                "Stdev: ", stdev_nFeature_RNA, "\n",
                "\n",
                "Percentage of molecules from mitochondrial genes per cell (percent.mt):\n",
                "Min: ", min_percent_mt, "\n",
                "Max: ", max_percent_mt, "\n",
                "Median: ", median_percent_mt, "\n",
                "Stdev: ", stdev_percent_mt)

  return(msg)
}

get_stats_df <- function(seurat_obj) {
  # Get the number of cells
  num_cells <- ncol(seurat_obj)

  # Get statistics about the number of molecules detected per cell
  nCount_RNA <- seurat_obj$nCount_RNA
  min_nCount_RNA <- min(nCount_RNA)
  max_nCount_RNA <- max(nCount_RNA)
  median_nCount_RNA <- median(nCount_RNA)
  stdev_nCount_RNA <- sd(nCount_RNA)

  # Get statistics about the number of genes detected per cell
  nFeature_RNA <- seurat_obj$nFeature_RNA
  min_nFeature_RNA <- min(nFeature_RNA)
  max_nFeature_RNA <- max(nFeature_RNA)
  median_nFeature_RNA <- median(nFeature_RNA)
  stdev_nFeature_RNA <- sd(nFeature_RNA)

  # Get statistics about mitochondrial genes expressed per cell
  percent_mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  min_percent_mt <- min(percent_mt)
  max_percent_mt <- max(percent_mt)
  median_percent_mt <- median(percent_mt)
  stdev_percent_mt <- sd(percent_mt)

  minimums <- c(num_cells, min_nCount_RNA, min_nFeature_RNA, min_percent_mt)
  maximums <- c(NA, max_nCount_RNA, max_nFeature_RNA, max_percent_mt)
  medians <- c(NA, median_nCount_RNA, median_nFeature_RNA, median_percent_mt)
  stdevs <- c(NA, stdev_nCount_RNA, stdev_nFeature_RNA, stdev_percent_mt)

  for (i in 2:4) {
    minimums[i] <- format(as.numeric(minimums[i]), digits = 6, scientific = ((minimums[i] <= 1e-6) || (minimums[i] >= 1e6)))
    maximums[i] <- format(as.numeric(maximums[i]), digits = 6, scientific = ((maximums[i] <= 1e-6) || (maximums[i] >= 1e6)))
    medians[i] <- format(as.numeric(medians[i]), digits = 6, scientific = ((medians[i] <= 1e-6) || (medians[i] >= 1e6)))
    stdevs[i] <- format(as.numeric(stdevs[i]), digits = 6, scientific = ((stdevs[i] <= 1e-6) || (stdevs[i] >= 1e6)))
  }

  df <- data.frame(minimums, maximums, medians, stdevs)
  rownames(df) <- c("Number of cells", "Number of molecules detected per cell", "Number of genes detected per cell", "Percentage of molecules from mitochondrial genes per cell")
  colnames(df) <- c("min", "max", "median", "stdev")

  return(df)
}

get_seurat_obj_from_file <- function(gene_counts_file) {
  gene_counts <- fread(gene_counts_file)

  # Turn any missing gene count into 0
  for (i in names(gene_counts)) {
    gene_counts[is.na(get(i)), (i):=0]
  }

  gene_counts <- gene_counts %>%
    remove_rownames %>%      # Row names are originally just numbers for each row
    column_to_rownames("V1") # Replace those numbers by the gene symbol or gene ID

  project_name <- gsub("_gene_count\\.csv$", "", basename(gene_counts_file))

  seurat_obj <- CreateSeuratObject(gene_counts, project = project_name)
  return(seurat_obj)
}

get_stats_msg_from_file <- function(gene_counts_file) {
  seurat_obj <- get_seurat_obj_from_file(gene_counts_file)
  stats_msg <- get_stats_msg(seurat_obj)
  return(stats_msg)
}

get_scatter_plot <- function(seurat_obj) {
  scatter_plot_obj <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident") +
    geom_smooth(method = "lm") + NoLegend() +
    labs(title = "Number of genes (nFeature_RNA) vs number of molecules (nCount_RNA)")
  return(scatter_plot_obj)
}

get_scatter_plot_from_file <- function(gene_counts_file) {
  seurat_obj <- get_seurat_obj_from_file(gene_counts_file)
  scatter_plot_obj <- get_scatter_plot(seurat_obj)
  return(scatter_plot_obj)
}

get_vln_plot <- function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  # use group.by to show overall info rather than cluster-level info (after running FindClusters())
  vln_plot_obj <- VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
  return(vln_plot_obj)
}

get_vln_plot_from_file <- function(gene_counts_file) {
  seurat_obj <- get_seurat_obj_from_file(gene_counts_file)
  vln_plot_obj <- get_vln_plot(seurat_obj)
  return(vln_plot_obj)
}

get_elbow_plot <- function(qc_seurat_obj, ndims) {
  num_cells <- ncol(qc_seurat_obj)
  elbow_plot_obj <- ElbowPlot(qc_seurat_obj, ndims = min(num_cells - 1, ndims), reduction = "pca")
  return(elbow_plot_obj)
}

get_umap_plot <- function(qc_seurat_obj) {
  umap_plot_obj <- DimPlot(qc_seurat_obj, reduction = "umap", label = TRUE, seed = 1)
  return(umap_plot_obj)
}

get_umap_harm_plot <- function(seurat_obj) {
  umap_harm_plot_obj <- DimPlot(seurat_obj, reduction = "umap.harm", group.by = c("orig.ident", "harm_cluster"), label = TRUE, seed = 1)
  return(umap_harm_plot_obj)
}

get_marker_gene_plot <- function(seurat_obj, features) {
  marker_gene_plot_obj <- DoHeatmap(seurat_obj, features = features)
  return(marker_gene_plot_obj)
}

get_feature_plot <- function(seurat_obj, feature) {
  feature_plot_obj <- FeaturePlot(seurat_obj, feature)
  return(feature_plot_obj)
}

get_dot_plot <- function(seurat_obj, feature, assay = "RNA") {
  dot_plot_obj <- DotPlot(seurat_obj, feature, assay = assay)
  return(dot_plot_obj)
}

write_stats_msg <- function(seurat_obj) {
  stats_msg <- get_stats_msg(seurat_obj)
  txt_filename <- "Statistics.txt"
  writeLines(stats_msg, txt_filename)
  return(txt_filename)
}

write_scatter_plot <- function(seurat_obj) {
  scatter_plot_obj <- get_scatter_plot(seurat_obj)
  pdf_filename <- "Scatter plot.pdf"
  ggsave(filename = pdf_filename, plot = scatter_plot_obj, width = 10, height = 6, device = "pdf")
  return(pdf_filename)
}

write_vln_plot <- function(seurat_obj) {
  vln_plot_obj <- get_vln_plot(seurat_obj)
  pdf_filename <- "Violin plot.pdf"
  ggsave(filename = pdf_filename, plot = vln_plot_obj, width = 10, height = 6, device = "pdf")
  return(pdf_filename)
}

write_elbow_plot <- function(qc_seurat_obj, ndims) {
  elbow_plot_obj <- get_elbow_plot(qc_seurat_obj, ndims)
  pdf_filename <- "Elbow plot.pdf"
  ggsave(filename = pdf_filename, plot = elbow_plot_obj, width = 10, height = 6, device = "pdf")
  return(pdf_filename)
}

write_umap_plot <- function(qc_seurat_obj) {
  umap_plot_obj <- get_umap_plot(qc_seurat_obj)
  pdf_filename <- "UMAP.pdf"
  ggsave(filename = pdf_filename, plot = umap_plot_obj, width = 10, height = 6, device = "pdf")
  return(pdf_filename)
}

write_umap_harm_plot <- function(seurat_obj) {
  umap_harm_plot_obj <- get_umap_harm_plot(seurat_obj)
  pdf_filename <- "UMAP_Harmony.pdf"
  ggsave(filename = pdf_filename, plot = umap_harm_plot_obj, width = 10, height = 6, device = "pdf")
  return(pdf_filename)
}

write_marker_gene_plot <- function(seurat_obj, features) {
  marker_gene_plot_obj <- get_marker_gene_plot(seurat_obj, features)
  pdf_filename <- "Top 20 marker genes.pdf"
  ggsave(filename = pdf_filename, plot = marker_gene_plot_obj, width = 10, height = 6, device = "pdf")
  return(pdf_filename)
}

write_gene_expression_plot <- function(seurat_obj, feature) {
  feature_plot_obj <- get_feature_plot(seurat_obj, feature)
  pdf_filename <- paste0("Gene expression UMAP plot (", feature, ").pdf")
  ggsave(filename = pdf_filename, plot = feature_plot_obj, width = 10, height = 6, device = "pdf")
  return(pdf_filename)
}

write_transcript_expression_plot <- function(seurat_obj, feature) {
  feature_plot_obj <- get_feature_plot(seurat_obj, feature)
  pdf_filename <- paste0("Transcript expression UMAP plot (", feature, ").pdf")
  ggsave(filename = pdf_filename, plot = feature_plot_obj, width = 10, height = 6, device = "pdf")
  return(pdf_filename)
}

write_gene_expression_dot_plot <- function(seurat_obj, feature) {
  dot_plot_obj <- get_dot_plot(seurat_obj, feature)
  pdf_filename <- paste0("Gene expression dot plot (", feature, ").pdf")
  ggsave(filename = pdf_filename, plot = dot_plot_obj, width = 10, height = 6, device = "pdf")
  return(pdf_filename)
}

write_transcript_expression_dot_plot <- function(seurat_obj, feature) {
  dot_plot_obj <- get_dot_plot(seurat_obj, feature, assay = "iso")
  pdf_filename <- paste0("Transcript expression dot plot (", feature, ").pdf")
  ggsave(filename = pdf_filename, plot = dot_plot_obj, width = 10, height = 6, device = "pdf")
  return(pdf_filename)
}

write_plot <- function(plot_obj, pdf_filename) {
  ggsave(filename = pdf_filename, plot = plot_obj, width = 10, height = 6, device = "pdf")
  return(pdf_filename)
}

perform_qc <- function(seurat_obj, min_nCount_RNA, max_nCount_RNA, min_nFeature_RNA, max_nFeature_RNA, max_percent_mt, npc_input, k_input, cluster_res_input) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  qc_seurat_obj <- subset(seurat_obj, subset = nCount_RNA >= min_nCount_RNA & nCount_RNA <= max_nCount_RNA &
                           nFeature_RNA >= min_nFeature_RNA & nFeature_RNA <= max_nFeature_RNA &
                           percent.mt <= max_percent_mt, return.null = TRUE)

  # If there are no cells left, stop performing further QC
  if (is.null(qc_seurat_obj)) {
    return(FALSE)
  }

  # Edge case: If the filtering is too stringent and only one cell is left, the remaining steps will fail!
  # Stop performing further QC if less than two cells are left
  if (ncol(qc_seurat_obj) < 2) {
    return(FALSE)
  }

  # Normalize data
  qc_seurat_obj <- NormalizeData(qc_seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

  # Identify highly variable features
  qc_seurat_obj <- FindVariableFeatures(qc_seurat_obj, selection.method = "vst", nfeatures = 2000)

  # Apply linear transformation
  all_genes <- rownames(qc_seurat_obj)
  qc_seurat_obj <- ScaleData(qc_seurat_obj, features = all_genes)

  # Perform PCA
  message("Doing PCA...")
  num_cells <- ncol(qc_seurat_obj)
  #npc <- min(num_cells - 1, 50)
  npc <- min(num_cells - 1, npc_input)
  qc_seurat_obj <- RunPCA(qc_seurat_obj, features = VariableFeatures(object = qc_seurat_obj), npcs = npc)
  message("After PCA...")

  cluster_res <- cluster_res_input
  k <- k_input

  # Cluster cells
  qc_seurat_obj <- FindNeighbors(qc_seurat_obj, dims = 1:npc, k.param = k)
  qc_seurat_obj <- FindClusters(qc_seurat_obj, resolution = cluster_res)

  qc_seurat_obj <- RunUMAP(qc_seurat_obj, n.neighbors = min(npc, 30), dims = 1:npc, seed.use = 42)

  return(qc_seurat_obj)
}

perform_qc_from_file <- function(gene_counts_file, min_nCount_RNA, max_nCount_RNA, min_nFeature_RNA, max_nFeature_RNA, max_percent_mt, npc_input, k_input, cluster_res_input) {
  seurat_obj <- get_seurat_obj_from_file(gene_counts_file)
  qc_seurat_obj <- perform_qc(seurat_obj, min_nCount_RNA, max_nCount_RNA, min_nFeature_RNA, max_nFeature_RNA, max_percent_mt, npc_input, k_input, cluster_res_input)
  return(qc_seurat_obj)
}