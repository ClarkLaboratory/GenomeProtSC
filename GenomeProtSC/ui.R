library(shiny)
library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
  dashboardHeader(title = "GenomeProtSC",
                  dropdownMenu(type = "messages",
                               tags$li(HTML('<li><a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/physiology/Parker-laboratory-Metabolic-Proteomics" target="_blank"><i class="fa fa-user"></i><h4>About us</h4><p>Parker Laboratory</p></a></li>')),
                               tags$li(HTML('<li><a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab" target="_blank"><i class="fa fa-user"></i><h4>About us</h4><p>Clark Laboratory</p></a></li>')),
                               tags$li(HTML('<li><a href="mailto:genomeprot@outlook.com" target="_blank"><i class="fa fa-question"></i><h4>Support</h4><p>genomeprot@outlook.com</p></a></li>'))
                  ), titleWidth = 300),
  # tabs
  dashboardSidebar(
    sidebarMenu(menuItem("Welcome", tabName = "welcome"),
                menuItem("README", tabName = "readme"),
                menuItem("1a. Run FLAMES", tabName = "flames", icon = icon("fire")),
                menuItem("1b. Filter & integrate scRNA-seq samples", tabName = "seurat_integration", icon = icon("gear")),
                menuItem("1c. Generate proteome database", tabName = "prot_db", icon = icon("database")),
                menuItem("2. Run proteomics analysis", tabName = "analyse_proteomics", icon = icon("search")),
                menuItem("3. Integrate proteomics & transcriptomics", tabName = "integration", icon = icon("arrows-turn-to-dots")),
                menuItem("4. Visualize results with IsoVis", tabName = "visualization_isovis", icon = icon("chart-simple"))
    ),
    width = 300
  ),
  # body
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$link(rel = "shortcut icon", href = "favicon.ico"),
      tags$link(rel = "apple-touch-icon", sizes = "180x180", href = "favicon.ico"),
      tags$link(rel = "icon", type = "image/png", sizes = "32x32", href = "/favicon-32x32.png"),
      tags$link(rel = "icon", type = "image/png", sizes = "16x16", href = "/favicon-16x16.png"),
      tags$style(HTML("
        .spinner {
          margin: 0 auto;
          width: 30px;
          height: 30px;
          border: 6px solid #ccc;
          border-top: 6px solid #333;
          border-radius: 50%;
          animation: spin 1s linear infinite;
        }

        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }

        .loading-container {
          display: none;
          text-align: center;
          margin-top: 20px;
        }

        #downloadResults {
          background-color: #4CAF50; /* Green */
          border: none;
          color: white;
          padding: 15px 32px;
          text-align: center;
          text-decoration: none;
          display: inline-block;
          font-size: 12px;
        }

        #downloadResults:disabled {
          background-color: #d3d3d3; /* Gray */
          color: #a9a9a9; /* Dark gray */
        }

        .spacing {
          margin-top: 20px;
        }

        /* Remove the 20px bottom margin for progress bars to save more space */
        .progress {
          margin-bottom: 0px;
        }

        /* Shrink the bottom margin for form groups from 20px to 10px to save more space */
        .form-group {
          margin-bottom: 10px;
        }
      ")),
      tags$script(HTML("
        Shiny.addCustomMessageHandler('clearClusterList', function(params) {
          var cluster_list = document.getElementById('cluster_list');
          if (!cluster_list)
              return;
          Shiny.unbindAll();
          cluster_list.innerHTML = '';
          Shiny.bindAll();
        });

        Shiny.addCustomMessageHandler('disableButtonOnly', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = true;
          button.style.backgroundColor = 'grey';
          button.style.borderColor = 'grey';
        });

        Shiny.addCustomMessageHandler('disableButton', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = true;
          button.style.backgroundColor = 'grey';
          button.style.borderColor = 'grey';
          document.getElementById(params.spinnerId).style.display = 'block';
        });

        Shiny.addCustomMessageHandler('enableButtonOnly', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = false;
          button.style.backgroundColor = '';
          button.style.borderColor = '';
        });

        Shiny.addCustomMessageHandler('enableButton', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = false;
          button.style.backgroundColor = '';
          button.style.borderColor = '';
          document.getElementById(params.spinnerId).style.display = 'none';
        });

        Shiny.addCustomMessageHandler('showStatusMessage', function(params) {
          var message = params.message;
          var color = params.color;
          var container = params.container;
          document.getElementById(container).innerText = message;
          document.getElementById(container).style.color = color;
        });

        Shiny.addCustomMessageHandler('clearStatusMessage', function(params) {
          var container = params.container;
          document.getElementById(container).innerText = '';
          document.getElementById(container).style.color = '';
        });

        var is_first_resize_ignored = false;
        window.addEventListener('message', handleMessage, false);
        function handleMessage(event)
        {
          let event_data = event.data;

          if (typeof(event_data) === 'string')
          {
            if (event_data.startsWith('To copy: '))
            {
              let text_to_copy = event_data.substring('To copy: '.length);
              let text_elem = document.getElementById('isovis_clipboard');
              if (text_elem)
              {
                text_elem.innerText = text_to_copy;
                text_elem.value = text_to_copy;
                text_elem.textContent = text_to_copy;
                text_elem.select();
                document.execCommand('copy');
                let isovis_window = document.getElementById('isovis_window');
                if (isovis_window)
                  isovis_window.focus();
              }
            }
            else if (event_data.startsWith('To resize: '))
            {
              if (!is_first_resize_ignored)
              {
                is_first_resize_ignored = true;
                return;
              }
              let new_window_height = Number.parseInt(event_data.substring('To resize: '.length));
              if (!Number.isNaN(new_window_height) && (new_window_height > 0))
              {
                if (new_window_height > window.innerHeight - 4)
                  new_window_height = window.innerHeight;
                else
                  new_window_height += 4;

                console.log('New iframe height:', new_window_height);

                let isovis_window = document.getElementById('isovis_window');
                if (isovis_window)
                    isovis_window.style.height = `${new_window_height}px`;
              }
            }
            else
              Shiny.onInputChange('selected_isovis_gene', event_data);
          }
          else if (typeof(event_data) === 'object')
            Shiny.onInputChange('marker_gene_info', event_data);
        }
      "))
    ),
    tabItems(
      tabItem(tabName = "welcome",
              fluidRow(
                column(12,
                       div(class = "box box-primary", style = "padding-right: 5%; padding-left: 5%; font-size:110%", 
                           div(class = "box-body", shiny::includeMarkdown("welcome-page-text.md")),
                           img(src = "images/workflow.png", width = "100%"),
                       )
                )
              )
      ),
      tabItem(tabName = "readme",
              fluidRow(
                column(12,
                       div(class = "box box-primary", style = "padding-right: 5%; padding-left: 5%; font-size:110%", 
                           div(class = "box-body", shiny::includeMarkdown("README.md")),
                       )
                )
              )
      ),
      tabItem(tabName = "flames",
              h2("Run FLAMES on your long-read single-cell FASTQ data"),
              h5("Outputs: A gene count matrix (CSV) and a transcript count matrix (CSV) for each sample, plus an isoform annotation GTF file containing all isoforms found to be expressed by the samples."),
              p("Notes: This step takes a long time to run and requires a minimum of 8 CPUs for larger datasets. For FASTQs, avoid complex file names as these are used as column headers later. A minimum of 50 cells in total are required for Seurat cell clustering."),
              fluidRow(
                column(6,
                       # Sequencing options
                       h3("Upload files:"),
                       fileInput("user_fastq_files", label = "Single-cell FASTQs (one per sample/replicate; only .fastq and .fastq.gz files are accepted)", multiple = TRUE, accept = c(".fastq", ".fastq.gz"), buttonLabel = "Browse..."),
                       numericInput("expected_cell_count", label = "Expected cells per sample (approx.):", value = 1000, step = 1),
                       numericInput("bambu_ndr", label = "Bambu NDR (0.0 <= NDR <= 1.0; enter any other value to let Bambu decide):", value = 0.5),
                       selectInput("organism", label = "Organism:", 
                                   choices = list("Human (H. sapiens)" = "human",
                                                  "Mouse (M. musculus)" = "mouse"), 
                                   selected = "human"),
                       actionButton("db_submit_button", "Submit", class = "btn btn-info")
                ),
                column(6,
                       h3("Download your results:"),
                       downloadButton("db_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "db-loading-container", class = "loading-container", div(class = "spinner")),
                       div(id = "db-status-msg-container")
                )
              )
      ),
      tabItem(tabName = "seurat_integration",
              h2("Use Seurat to filter and integrate scRNA-seq samples"),
              h5("Outputs: A text file showing the statistics of the sample being looked at, a PDF for each plot drawn, and an RDS file of the Seurat object being shown. For quality control and sample integration specifically, a tab-separated text file showing the sample, cell barcode and cluster name of each cell in the Seurat object is also included in the output."),
              fluidRow(
                column(4,
                       h3("Upload files:"),
                       fileInput("gene_count_files", label = "Gene count CSVs from the FLAMES module (must be named as '<sample ID>_gene_count.csv'):", multiple = TRUE, accept = c(".csv"), buttonLabel = "Browse..."),
                       fileInput("transcript_count_files", label = "Transcript count CSVs from the FLAMES module (must be named as '<sample ID>_transcript_count.csv'):", multiple = TRUE, accept = c(".csv"), buttonLabel = "Browse..."),
                       actionButton("view_stats_button", "View statistics of gene count CSVs", class = "btn btn-info"),
                       h1(),

                       actionButton("show_qc_part", "Show 'quality control parameters' menu ▲", class = "btn btn-secondary"),
                       actionButton("hide_qc_part", "Hide 'quality control parameters' menu ▼", class = "btn btn-secondary"),

                       h3(id = "qc_heading", "Quality control parameters:"),
                       actionButton("show_adv_qc_params_button", "Show advanced quality control parameters", class = "btn btn-warning"),
                       actionButton("hide_adv_qc_params_button", "Hide advanced quality control parameters", class = "btn btn-warning"),
                       numericInput("min_nCount_RNA",
                                    label = "Minimum number of molecules detected per cell",
                                    value = 500),
                       numericInput("max_nCount_RNA",
                                    label = "Maximum number of molecules detected per cell",
                                    value = 10000),
                       numericInput("min_nFeature_RNA",
                                    label = "Minimum number of genes detected per cell",
                                    value = 500),
                       numericInput("max_nFeature_RNA",
                                    label = "Maximum number of genes detected per cell",
                                    value = 10000),
                       numericInput("max_percent_mt",
                                    label = "Maximum percentage of molecules from mitochondrial genes per cell",
                                    value = 10),
                       numericInput("npc_input",
                                    label = "Number of principal components (PCs) when performing principal component analysis (note: this value is best determined by looking at the elbow plot; the actual value used will be at most one less than the total number of cells in the filtered object)",
                                    value = 20),
                       numericInput("k_input",
                                    label = "Number of neighbours when constructing the k-nearest neighbours graph",
                                    value = 20),
                       numericInput("cluster_res_input",
                                    label = "Cluster resolution when finding clusters",
                                    value = 0.9),
                       actionButton("view_qc_stats_button", "View statistics of gene count CSVs after quality control", class = "btn btn-info"),
                       h1(),

                       actionButton("show_integ_part", "Show 'Harmony sample integration parameters' menu ▲", class = "btn btn-secondary"),
                       actionButton("hide_integ_part", "Hide 'Harmony sample integration parameters' menu ▼", class = "btn btn-secondary"),

                       h3(id = "harmony_integ_heading", "Harmony sample integration parameters:"),
                       actionButton("show_adv_harmony_integ_params_button", "Show advanced Harmony integration parameters", class = "btn btn-warning"),
                       actionButton("hide_adv_harmony_integ_params_button", "Hide advanced Harmony integration parameters", class = "btn btn-warning"),
                       numericInput("tau",
                                    label = "tau: Expected number of cells per cluster; used to protect against over-clustering small datasets with large ones",
                                    value = 0),
                       numericInput("theta",
                                    label = "theta: Control cluster diversity (0 = no diversity, larger = more diverse clusters)",
                                    value = 2),
                       numericInput("npcs",
                                    label = "npcs: Number of principal components (PCs) to compute on the gene and transcript counts",
                                    value = 50),
                       numericInput("lambda",
                                    label = "lambda: Ridge regression penalty (larger = more protection against over-correction)",
                                    value = 1),
                       numericInput("sigma",
                                    label = "sigma: Width of soft k-means clusters (smaller = soft clustering becomes more similar to hard clustering, larger = cells assigned to more clusters)",
                                    value = 0.1),
                       numericInput("nclust",
                                    label = "nclust: Number of clusters (set to -1 for Harmony's default value)",
                                    value = -1),
                       numericInput("block.size",
                                    label = "block.size: Proportion of cells to update during clustering (between 0 and 1; larger values may be faster but less accurate)",
                                    value = 0.05),
                       numericInput("max.iter.harmony",
                                    label = "max.iter.harmony: Maximum number of rounds to run Harmony",
                                    value = 10),
                       numericInput("max.iter.cluster",
                                    label = "max.iter.cluster: Maximum number of rounds to run clustering at each round of Harmony",
                                    value = 20),
                       numericInput("epsilon.cluster",
                                    label = "epsilon.cluster: Convergence tolerance for the clustering round of Harmony (non-numeric or non-positive input = never stop early)",
                                    value = 0.00001),
                       numericInput("epsilon.harmony",
                                    label = "epsilon.harmony: Convergence tolerance for Harmony (non-numeric or non-positive input = never stop early)",
                                    value = 0.0001),
                       actionButton("perform_integration_button", "Perform sample integration", class = "btn btn-info"),
                       p(id = "harmony_integ_footer_1", "The default values of the Harmony integration parameters are taken from Seurat's 'Harmony Integration' documentation webpage: https://satijalab.org/seurat/reference/harmonyintegration"),
                       p(id = "harmony_integ_footer_2", "Harmony paper: Korsunsky, I., Millard, N., Fan, J. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16, 1289–1296 (2019). https://doi.org/10.1038/s41592-019-0619-0"),
                       h1(),

                       actionButton("show_rename_clusters_part", "Show 'rename clusters' menu ▲", class = "btn btn-secondary"),
                       actionButton("hide_rename_clusters_part", "Hide 'rename clusters' menu ▼", class = "btn btn-secondary"),

                       h3(id = "rename_clusters_heading", "Rename clusters:"),
                       p(id = "rename_clusters_info", "Ensure new cluster names do not contain any commas (',') or double quotes ('\"')."),
                       div(id = "cluster_list"),
                       actionButton("rename_clusters_button", "Update Seurat object and plots with renamed clusters", class = "btn btn-info"),
                       h1(),

                       actionButton("show_marker_genes_part", "Show 'find marker genes of one cell cluster' menu ▲", class = "btn btn-secondary"),
                       actionButton("hide_marker_genes_part", "Hide 'find marker genes of one cell cluster' menu ▼", class = "btn btn-secondary"),

                       h3(id = "find_marker_genes_heading", "Find marker genes of specific cell clusters:"),
                       p(id = "find_marker_genes_info", "To specify more than one cluster, separate cluster names by commas (',')."),
                       textInput("ident.1",
                                 label = "ident.1: Cell cluster(s) to find marker genes from",
                                 value = ""),
                       textInput("ident.2",
                                 label = "ident.2: Cell cluster(s) to be compared against (if empty, ident.2 represents all cells other than ident.1)",
                                 value = ""),
                       numericInput("p_adj.thresh",
                                    label = "p_adj.thresh: Only return markers with an adjusted p-value below a threshold",
                                    value = 0.05),
                       numericInput("logfc.threshold",
                                    label = "logfc.threshold: Only test genes with at least a specific log fold change difference between two groups of cells",
                                    value = 0.1),
                       numericInput("min.pct",
                                    label = "min.pct: Only test genes with at least a specific fraction of cells in ident.1 and ident.2",
                                    value = 0.01),
                       numericInput("min.diff.pct",
                                    label = "min.diff.pct: Only test genes with at least a specific fractional difference between ident.1 and ident.2 (default: negative value = don't use differences to filter genes for testing)",
                                    value = -1),
                       checkboxInput("only.pos",
                                     label = "only.pos: Only return positive markers",
                                     value = FALSE),
                       actionButton("find_marker_genes_button", "Update output zip with list of marker genes", class = "btn btn-info"),
                       h1(),

                       actionButton("show_all_marker_genes_part", "Show 'find marker genes of all cell clusters' menu ▲", class = "btn btn-secondary"),
                       actionButton("hide_all_marker_genes_part", "Hide 'find marker genes of all cell clusters' menu ▼", class = "btn btn-secondary"),

                       h3(id = "find_all_marker_genes_heading", "Find marker genes of all cell clusters:"),
                       numericInput("p_adj.thresh_all",
                                    label = "p_adj.thresh: Only return markers with an adjusted p-value below a threshold",
                                    value = 0.05),
                       numericInput("logfc.threshold_all",
                                    label = "logfc.threshold: Only test genes with at least a specific log fold change difference between two groups of cells",
                                    value = 0.1),
                       numericInput("min.pct_all",
                                    label = "min.pct: Only test genes with at least a specific fraction of cells in ident.1 and ident.2",
                                    value = 0.01),
                       numericInput("min.diff.pct_all",
                                    label = "min.diff.pct: Only test genes with at least a specific fractional difference between ident.1 and ident.2 (default: negative value = don't use differences to filter genes for testing)",
                                    value = -1),
                       checkboxInput("only.pos_all",
                                     label = "only.pos: Only return positive markers",
                                     value = FALSE),
                       actionButton("find_all_marker_genes_button", "Update output zip with list of all marker genes", class = "btn btn-info")
                ),
                column(8,
                       selectInput("sample_selector", "Select a sample:", choices = NULL),
                       h3("Download your results:"),
                       downloadButton("seurat_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "seurat-loading-container", class = "loading-container", div(class = "spinner")),
                       div(id = "seurat-status-msg-container"),

                       actionButton("show_statistics_output_part", "Show statistics of samples after quality control ▲", class = "btn btn-secondary"),
                       actionButton("hide_statistics_output_part", "Hide statistics of samples after quality control ▼", class = "btn btn-secondary"),
                       h4(id = "statistics_heading", "Statistics of the sample / integrated samples:"),
                       tableOutput("seurat_statistics_msg_container"),
                       h4(id = "scatter_plot_heading", "Scatter plot of unique genes vs reads per cell in the sample:"),
                       plotOutput("scatter_plot"),
                       h4(id = "vln_plot_heading", "Violin plot of single-cell statistics of the sample:"),
                       plotOutput("vln_plot"),
                       h4(id = "elbow_plot_heading", "Elbow plot of sample from principal component analysis (PCA):"),
                       plotOutput("elbow_plot"),

                       h4(id = "umap_plot_heading", "UMAP plot of clusters:"),
                       plotOutput("umap_plot"),
                       h4(id = "marker_gene_plot_heading", "Heatmap of top 20 marker genes:"),
                       plotOutput("marker_gene_plot")
                )
              )
      ),
      tabItem(tabName = "prot_db",
              h2("Generate a proteome database from Bambu transcript annotations"),
              h5("Outputs: 'proteome_database.fasta', which stores the amino acid sequences of all open reading frames (ORFs) found in the isoform annotation GTF file; 'proteome_database_metadata.txt', which contains information on each ORF found; and 'proteome_database_transcripts.gtf', which is the annotation of all transcripts used to generate 'proteome_database.fasta'."),
              fluidRow(
                column(6,
                       # Sequencing options
                       h3("Upload files:"),
                       fileInput("prot_db_gtf_file", label = "Isoform annotation GTF file from the FLAMES module", NULL, buttonLabel = "Browse...", multiple = FALSE, accept = c(".gtf")),
                       selectInput("prot_db_organism", label = "Organism:",
                                   choices = list("Human (H. sapiens)" = "human",
                                                  "Mouse (M. musculus)" = "mouse"),
                                   selected = "human"),
                       h3("Protein database options:"),
                       numericInput("min_orf_length",
                                    label = "ORF length (amino acids):",
                                    value = 30),
                       h5(tags$b("Find short (10 to 'ORF length' amino acids) ORFs in UTRs of reference transcripts:")),
                       checkboxInput("user_find_utr_5_orfs", label = "Upstream 5' ORFs",
                                     value = FALSE, width = NULL),
                       checkboxInput("user_find_utr_3_orfs", label = "Downstream 3' ORFs",
                                     value = FALSE, width = NULL),
                       actionButton("prot_db_submit_button", "Submit", class = "btn btn-info")
                ),
                column(6,
                       h3("Download your results:"),
                       downloadButton("prot_db_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "prot-db-loading-container", class = "loading-container", div(class = "spinner")),
                       div(id = "prot-db-status-msg-container")
                )
              )
      ),
      tabItem(tabName = "analyse_proteomics",
              h2("Perform proteomics searches using your custom database"),
              h5("This module is under development. Currently, users need to run FragPipe externally and return here to normalize peptide counts once complete."),
              fluidRow(
                column(6,
                       h3("Apply variance stabilizing normalization to peptide counts:"),
                       fileInput("peptide_counts_file", label = "Upload peptide counts to normalize (e.g. 'report.pr_matrix.tsv'):", buttonLabel = "Browse...", multiple = FALSE, accept = c(".csv", ".tsv", ".txt")),
                       actionButton("peptide_counts_submit_button", "Normalize peptide counts", class = "btn btn-info")
                ),
                column(6,
                       h3("Download your results:"),
                       downloadButton("vsn_download_button", "Download normalized peptide counts (CSV)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "vsn-loading-container", class = "loading-container", div(class = "spinner")),
                       div(id = "vsn-status-msg-container")
                )
              )
      ),
      tabItem(tabName = "integration", 
              h2("Integrate proteomics results with transcriptomics"),
              h5("Outputs: A combined GTF and BED12s of peptides, ORFs and transcripts for visualization. Also includes a CSV file listing, for each transcript, the level of evidence supporting its translation and the cell clusters that have expressed it, plus summary data and a report."),
              fluidRow(
                column(6,
                       h3("Upload transcriptomics and proteomics data files to integrate:"),
                       fileInput("user_proteomics_file", "Upload proteomics results:", NULL, buttonLabel = "Browse...", multiple = FALSE, accept = c(".csv", ".tsv", ".txt")),
                       fileInput("user_fasta_file", "Upload 'proteome_database.fasta':", NULL, buttonLabel = "Browse...", multiple = FALSE, accept = c(".fasta")),
                       fileInput("user_metadata_file", "Upload 'proteome_database_metadata.txt':", NULL, buttonLabel = "Browse...", multiple = FALSE, accept = c(".txt")),
                       fileInput("user_post_gtf_file", "Upload 'proteome_database_transcripts.gtf':", NULL, buttonLabel = "Browse...", multiple = FALSE, accept = c(".gtf")),
                       fileInput("user_rds_file", "Upload Seurat object RDS file:", NULL, buttonLabel = "Browse...", multiple = FALSE, accept = c(".rds")),
                       actionButton("integ_submit_button", "Submit", class = "btn btn-info"),

                       h1(),

                       h3("Upload a combined annotations GTF file and a marker gene info file to integrate:"),
                       fileInput("user_gtf_file", "Upload 'combined_annotations.gtf':", NULL, buttonLabel = "Browse...", multiple = FALSE, accept = c(".gtf")),
                       fileInput("user_marker_file", "Upload 'marker_genes.csv' from module 1b:", NULL, buttonLabel = "Browse...", multiple = FALSE, accept = c(".csv")),
                       actionButton("integ_marker_submit_button", "Submit", class = "btn btn-info"),
                ),
                column(6,
                       h3("Download transcriptomics and proteomics integration results:"),
                       downloadButton("integ_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "integ-loading-container", class = "loading-container", div(class = "spinner")),
                       div(id = "integ-status-msg-container"),

                       h1(),

                       h3("Download marker gene integration results:"),
                       downloadButton("integ_marker_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "integ-marker-loading-container", class = "loading-container", div(class = "spinner")),
                       div(id = "integ-marker-status-msg-container")
                )
              )
      ),
      tabItem(tabName = "visualization_isovis",
              fluidRow(
                column(12, 
                       h2("Visualize results obtained from the pipeline with IsoVis", align = "left"),
                       tags$textarea(id = "isovis_clipboard", rows = 1, cols = 1, "", style = "position: absolute; z-index: -1"), # Use a textarea HTML element hidden behind the IsoVis window to store tooltip text 'copied' from IsoVis and to do the actual copying
                       tags$iframe(id = "isovis_window", src = "https://isomix.org/isovis/", style = "width: 100%; height: 80vh;")
                )
              ),
              fluidRow(
                column(6,
                       h3("Upload files to visualize outside of IsoVis:"),
                       h4("Certain files cannot be visualized inside IsoVis, including Seurat cell cluster UMAP plots. They have to be uploaded into GenomeProtSC instead."),
                       h4("Once uploaded, information relevant to the gene selected in IsoVis will be automatically shown from these files."),
                       fileInput("transcript_info_file_isovis", label = "Upload 'transcript_expression_info.csv' file:", buttonLabel = "Browse...", multiple = FALSE, accept = c(".csv")),
                       fileInput("seurat_rds_file_isovis", label = "Seurat object to view (required if 'transcript_expression_info.csv' is uploaded):", buttonLabel = "Browse...", multiple = FALSE, accept = c(".rds")),
                       actionButton("vis_submit_button_isovis", "Submit", class = "btn btn-info"),
                ),
                column(6,
                       h3("Download your results:"),
                       downloadButton("vis_download_button_isovis", "Download plots shown below (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "vis-loading-container_isovis", class = "loading-container", div(class = "spinner")),
                       div(id = "vis-status-msg-container_isovis")
                )
              ),
              fluidRow(
                column(12,
                       selectInput("transcript_selector_isovis", "Select a transcript:", choices = NULL),

                       h4(id = "transcript_expression_info_heading_isovis", "Isoform expression information"),
                       tableOutput("transcript_expression_info_container_isovis"),
                       h4(id = "marker_gene_info_heading_isovis", "Marker gene information"),
                       tableOutput("marker_gene_info_container_isovis"),
                       h4(id = "umap_vis_plot_heading_isovis", "UMAP plot of cell clusters"),
                       plotOutput("umap_vis_plot_isovis"),
                       h4(id = "gene_expression_plot_heading_isovis", "UMAP plot of gene expression in cells"),
                       plotOutput("gene_expression_plot_isovis"),
                       h4(id = "gene_expression_dot_plot_heading_isovis", "Dot plot of gene expression in cell clusters"),
                       plotOutput("gene_expression_dot_plot_isovis"),
                       h4(id = "transcript_expression_plot_heading_isovis", "UMAP plot of transcript expression in cells"),
                       plotOutput("transcript_expression_plot_isovis"),
                       h4(id = "transcript_expression_dot_plot_heading_isovis", "Dot plot of transcript expression in cell clusters"),
                       plotOutput("transcript_expression_dot_plot_isovis")
                )
              )
      )
    )
  ),  
  skin = "blue"
)
