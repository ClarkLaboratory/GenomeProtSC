# GenomeProtSC: an integrated proteogenomics analysis platform for long-read single-cell data

## ⚠️ WARNING ⚠️: This repository is currently undergoing heavy updates and may contain code that does not fully work!

## Contents

- [Installation](#installation)
- [General usage](#general-usage)
  - [1. Proteome database generation](#1-proteome-database-generation)
  - [2. Proteomics](#2-proteomics)
  - [3. Integration](#3-integration)
  - [4. Visualization](#4-visualization)
- [Detailed input and output descriptions](#detailed-input-and-output-descriptions)
  - [1a. FLAMES](#1a-flames)
  - [1b. Seurat](#1b-seurat)
  - [1c. Proteome database generation](#1c-proteome-database-generation)
  - [2. Proteomics](#2-proteomics-1)
  - [3. Integration](#3-integration-1)
  - [4. Visualization](#4-visualization-1)

## Installation
Note: Option 1 and 2 are under development and not yet available.

### Option 1 (recommended): Run the R Shiny application with Docker
Make sure you have [Docker](https://docs.docker.com/engine/install/) installed and the application running in the background before you begin.

Open your terminal application and run:
```
docker run --rm -p 3838:3838 josieg/genomeprotsc:v1
```
This will take approximately 10-20 minutes to download the Docker image the first time the app is run.
The --rm removes the container after it's stopped and the -p 3838:3838 maps your local port 3838 to the same port inside the container.

To **access the local R Shiny application**, navigate to this link on your web browser http://0.0.0.0:3838.

You can now upload all files and run the steps in your web browser. Although the app is running through a web browser, no files are being uploaded to the internet and everything will be run locally.

To stop the container, close the web browser tab and head back to the terminal where Docker is running and press ctrl+c.

### Option 2 (recommmended for downstream analysis only): Access GenomeProtSC online
https://genomeprotSC.researchsoftware.unimelb.edu.au/

### Option 3: Install and run the R Shiny application locally

#### Ensure your operating system is up to date:

```
sudo apt update
sudo apt install software-properties-common
```

#### Install R:

Click the following link and download the version for your operating system (Windows, macOS or Linux): https://cran.r-project.org/index.html

#### Install RStudio:

Click the following link and download RStudio for Desktop: https://posit.co/download/rstudio-desktop/

#### Install Python and cd-hit:

```
# if you don't already have python3
sudo apt install python
sudo apt install python3-pip

pip install biopython
pip install py-cdhit
pip install peptides

# install cd-hit
sudo apt install cd-hit
```

If you do not have sudo permissions, you can instead use a package manager such as Conda or compile the dependencies from source.

#### Install the required R and Bioconductor packages:

```
# either run the following commands with 'sudo' or open RStudio and install
# install Bioconductor and FLAMES (replace '3.22' with the latest Bioconductor version)

sudo R -e 'install.packages("BiocManager")'
sudo R -e 'BiocManager::install(version = "3.22")'
sudo R -e 'BiocManager::install("FLAMES")'

# install base R packages
sudo R -e 'install.packages(c("shiny", "shinyjs", "shinythemes", "shinydashboard", "data.table", "dplyr", "tidyr", "readr", "tibble", "purrr", "forcats", "phylotools", "markdown", "rmarkdown", "ggplot2", "ggrepel", "devtools", "optparse", "reshape2", "stringr", "stringi", "RColorBrewer", "scales", "gplots"))'

# install other Bioconductor packages
sudo R -e 'BiocManager::install(c("retracklayer", "ORFik", "GenomicAlignments", "GenomicFeatures", "GenomicRanges", "Biostrings", "mygene", "ORFik", "patchwork", "Rsamtools", "SummarizedExperiment", "tximport", "patchwork", "vsn"), ask = F)'

# install required genomes from Bioconductor
sudo R -e 'BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm39"))'
```

Note that FLAMES takes a significant amount of time to install.

#### Export genome files:

```
# in R or R studio, export BSGenome objects as FASTA files
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm39)

# files need to be saved in the GenomeProtSC/refs directory for FLAMES

# human
genomedb <- BSgenome.Hsapiens.UCSC.hg38
export(genomedb, "path/to/GenomeProtSC/refs/human.fasta", verbose = T, compress = F, format = "fasta")

# mouse
genomedb <- BSgenome.Mmusculus.UCSC.mm39
export(genomedb, "path/to/GenomeProtSC/refs/mouse.fasta", verbose = T, compress = F, format = "fasta")
```

#### Clone this repository:

```
git clone https://github.com/ClarkLaboratory/GenomeProtSC.git
```

#### Decompress the UniProt + OpenProt reference files:

```
cd GenomeProtSC/data
gunzip openprot_uniprotDb_hs.txt.gz
gunzip openprot_uniprotDb_mm.txt.gz
```

#### Run the GenomeProtSC app:

If using RStudio:
- Open either the `server.R` file or the `ui.R` file inside the GenomeProtSC folder in RStudio
- Click the 'Run App' button in the top right corner

If using the command line:

```
# provide the path to the GenomeProtSC app
Rscript -e "shiny::runApp('/path/to/app/GenomeProtSC/', host = '0.0.0.0', port = 3838)"
```

## General usage

GenomeProtSC is an integrated proteogenomics platform. The application consists of four components: 1) proteome database generation, 2) proteomics (performed partially externally to GenomeProtSC), 3) integration, and 4) visualization. The proteome database generation component is split into three modules: 1a) FLAMES long-read scRNA-seq module, 1b) Seurat module, and 1c) Proteome database generation module.

## 1. Proteome database generation

The first module processes long-read single-cell data with FLAMES and generates a custom proteome database to perform proteomics searches. The module accepts RNA sequencing FASTQ files from 10X single-cell long-read sequencing. The main outputs from this module are the results from FLAMES (single-cell transcriptomics), a FASTA file with candidate protein sequences and a metadata file with details of each candidate protein.

GenomeProtSC currently supports open reading frame (ORF) identification and database generation for human and mouse data. Users can specify an option to include short upstream ORF (uORF) and downstream ORF (dORF) protein sequences > 10 amino acids (AA). Protein sequences are generated based on a user defined minimum length set to > 30 AA by default.

### 1a. Run FLAMES on long-read single-cell FASTQ data

#### Inputs:

- FASTQ file(s) containing sequencing reads (one per sample; can be gzipped)
- Expected number of cells per sample (approximate)
- Bambu NDR (novel discovery rate; affects the threshold where novel transcripts are called)
- Reference annotation GTF file (automatically supplied by GenomeProtSC; currently limited to GENCODE/Ensembl human and mouse GTFs)

#### Outputs (zipped):

- Isoform annotation GTF file containing all isoforms found to be expressed in the samples
- Quantification files (one gene count CSV & one transcript count CSV per sample)

### 1b. Use Seurat to filter, integrate and cluster scRNA-seq samples

#### Inputs:

- Gene count CSV(s) from the FLAMES module; one per sample
- Transcript count CSV(s) from the FLAMES module; one per sample
- Various parameters

#### Outputs (zipped):

- Text file containing statistics of the sample being investigated
- Figures (in PDF format) for each plot generated by GenomeProtSC, including violin plots of the number of genes and molecules detected per cell; a scatter plot showing the number of genes per cell against the number of molecules per cell; and a UMAP plot showing cell clusters (if quality control or sample integration were performed).
- Seurat R object
- If quality control filtering and/or sample integration were performed, a tab-separated text file showing the sample, cell barcode and cluster name of each cell in the Seurat object is included in the output.
- If marker genes were calculated, a CSV file containing marker gene information is included in the output.

### 1c. Generate the proteome database

#### Inputs:

- Isoform annotation GTF file from the FLAMES module
- The species the long-read scRNA-seq data came from (currently limited to human and mouse)
- ORF length (the minimum number of amino acids ORFs should have)
- Whether to look for ORFs present in the 5' and/or 3' untranslated regions of reference transcripts

#### Outputs (zipped):

- The generated proteome database (a FASTA file containing the amino acid sequences of all ORFs found in the isoform annotation GTF file)
- A metadata file (contains information on each ORF found)
- A GTF annotation of all transcripts used to generate the proteome database

## 2. Proteomics

The matching of peptides to ORFs is currently performed externally to GenomeProtSC using the database produced by module 1c. It is recommended to use FragPipe to process the proteomics data with the generated proteome database.

### Apply variance stabilizing normalization to peptide counts

#### Input:

A peptide counts file obtained from peptide quantification with FragPipe

#### Output:

A normalized peptide counts CSV file

## 3. Integration

This module integrates proteomics and transcriptomics data and consists of two steps. The first step involves associating peptides back to transcript isoforms and mapping them to spliced genomic coordinates for downstream visualization, with the main output being a GTF annotation file containing the integration of identified transcripts, peptides and ORFs. The second (optional) step is to integrate marker gene information into the GTF annotation file outputted in the first step.

#### Inputs:

- Proteomics results from FragPipe
- All 3 output files from the proteome database generation module (module 1c)
- Seurat R object from the Seurat module (module 1b)
- Integrated combined GTF annotation (generated by the first step and required for the optional second step)
- Marker gene CSV file (module 1b, optional)

#### Outputs (zipped):

- An integrated combined GTF annotation of identified transcripts, peptides and ORFs
- A CSV file containing peptide-to-transcript mappings with spliced genomic coordinates
- BED12 files for visualizing transcripts, peptides and ORFs in the UCSC Genome Browser
- A GTF file for visualizing transcripts and ORFs in IsoVis
- An HTML report summarizing the identified transcripts, peptides and ORFs
- A transcript expression information file (CSV) showing, for each transcript, the cell clusters it is expressed in and the level of evidence supporting its translation

## 4. Visualization

The visualization module uses IsoVis to generate peptide mapping plots along transcript isoforms with quantitative peptide intensities and uses Seurat to create UMAP plots of the expression of selected genes and isoforms in cell clusters. This allows users to visualize transcript and peptide abundance across different experimental conditions. This module requires the combined GTF annotation file generated in the integration module, and optionally the Seurat R object from module 1b, the peptide intensities from external proteomics analysis, and the transcript expression information CSV file from the integration module. Inside IsoVis, we have included gene filtering options to quickly search for features of interest, such as genes with unknown functions.

#### Inputs (in IsoVis):

- The combined GTF annotation file from the integration module (required)
- Peptide intensities (optional)

#### Inputs (in GenomeProtSC):

- The transcript expression information CSV file from the integration module (optional)
- Seurat R object from the Seurat module (required if transcript expression information is uploaded)

#### Outputs (zipped):

- Figures (in PDF format) for each plot generated in GenomeProtSC, including a UMAP plot showing cell cluster annotations, and plots of gene and transcript expression levels at the cell / cluster level

## Detailed input and output descriptions

### 1a. FLAMES

| Input | Type | Required? | Description |
|-|-|-|-|
| Sequencing data | (gzipped) FASTQ(s) | Yes | Long-read single-cell sequencing reads. Upload one per sample. |

| Output | Type | Filename | Description |
|-|-|-|-|
| Gene counts | CSV | <sample_name>_gene_count.csv (e.g. if the user uploaded early.fastq, they will get early_gene_count.csv) | Gene count file with cell barcodes as columns and genes as rows. One for each sample. |
| Transcript counts | CSV | <sample_name>_transcript_count.csv | Transcript count file with cell barcodes as columns and transcripts as rows. One for each sample. |
| Isoform annotation | GTF | isoform_annotated.gtf | Annotation of all isoforms found to be expressed by the samples. |

### 1b. Seurat

| Input | Type | Required? | Description |
|-|-|-|-|
| Gene counts | CSV(s) | Yes | Generated in module 1a. Gene count file with cell barcodes as columns and genes as rows. Upload one for each sample. |
| Transcript counts | CSV(s) | Yes | Generated in module 1a. Transcript count file with cell barcodes as columns and transcripts as rows. Upload one for each sample. |

| Output | Type | Filename | Description |
|-|-|-|-|
| Sample statistics | TXT | Statistics.txt | A textual description of the statistics of the number of genes and molecules detected per cell, and the percentage of features from mitochondrial genes per cell. |
| Exported figures | PDFs | <plot_name>.pdf | One PDF for each figure drawn. |
| Seurat object | RDS | seurat_object.rds | Seurat object of the selected sample. Can be loaded into R with the `readRDS()` function. |
| Samples, barcodes and clusters | Tab-separated text (TXT) | sample_barcode_cluster.txt | The sample, barcode and cluster of each cell. Applicable if quality control or sample integration were done. |
| Marker gene information | CSV | marker_genes.csv | Metrics for each marker gene found (e.g. adjusted p-value, proportion of cells expressing the gene in the cluster it is a marker of, average log2(fold change)) and the clusters they are marker genes of. Applicable if marker genes were calculated. |

### 1c. Proteome database generation

| Input | Type | Required? | Description |
|-|-|-|-|
| Isoform annotation | GTF | Yes | Generated in module 1a. Annotation of all isoforms found to be expressed in the scRNA-seq samples. |

| Output | Type | Filename | Description |
|-|-|-|-|
| Proteome database | FASTA | proteome_database.fasta | Amino acid sequences of all ORFs in the input. |
| Proteome database metadata | TXT | proteome_database_metadata.txt | Information on each ORF in the database. |
| Proteome database transcripts | GTF | proteome_database_transcripts.gtf | Annotation of all transcripts used to generate the database. |

Proteome FASTA examples (unannotated ORFs are denoted by "ORF_"):
```
>protein_accession|CO=genomic_coordinates GA=gene_accession GN=gene_name TA=transcript_accession
MCGNNMSAPMPAVVPAARKATAAVIFLHGLGDTGHGWAEAFAGIKSPHIKYICPHAPVMPVTLNMNMAMPSWFDIVGLSPDSQEDESGIKQAAETVKALIDQEVKNGIPSNRIILGGFSQGPINSANRDISVLQCHGDCDPLVPLMFGSLTVERLKALINPANVTFKIYEGMMHSSCQQEMMDVKHFIDKLLPPID
>P10711|CO=chr1:4928137-4966584 GA=ENSMUSG00000033813.16 GN=Tcea1 TA=ENSMUST00000081551.14
MEDEVVRIAKKMDKMVQKKNAAGALDLLKELKNIPMTLELLQSTRIGMSVNALRKQSTDEEVTSLAKSLIKSWKKLLDGPSTDKDPEEKKKEPAISSQNSPEAREESSSSSNVSSRKDETNARDTYVSSFPRAPSTSDSVRLKCREMLAAALRTGDDYVAIGADEEELGSQIEEAIYQEIRNTDMKYKNRVRSRISNLKDAKNPNLRKNVLCGNIPPDLFARMTAEEMASDELKEMRKNLTKEAIREHQMAKTGGTQTDLFTCGKCKKKNCTYTQVQTRSADEPMTTFVVCNECGNRWKFC
>ORF_3|CO=chr2:53029193-53081430 GA=ENSMUSG00000061136.17 GN=Prpf40a TA=ENSMUST00000209364.3
MQATPSEAGGESPQSCLSVSRSDWTVGKPVSLLAPLIPPRSSGQPLPFGPGGRQPLRSLLVGMCSGSGRRRSSLSPTMRPGTGAERGGLMMGHPGMHYAPMGMHPMGQRANMPPVPHGMMPQMMPPMGG
```

**Note:** UniProt or RefSeq accessions are retained for annotated proteins.

#### Open reading frame (ORF) category definitions:

| Type  | Definition    |
|-|-|
| CDS | Annotated in UniProt or RefSeq |
| 5UTR | Coordinates are within the 5' UTR region of an mRNA transcript |
| 3UTR | Coordinates are within the 3' UTR region of an mRNA transcript |
| 5UTR:CDS | Start site is within the 5' UTR region and stop site is within the CDS region of an mRNA transcript |
| gene_overlap | Encoded by a transcript that overlaps a region with annotated protein-coding genes |
| intergenic | Encoded by a transcript that does not overlap a region with annotated protein-coding genes |

### 2. Proteomics

We recommend installing and running [FragPipe](https://github.com/Nesvilab/FragPipe) for analysing mass spectrometry-based proteomics data.

#### Peptide detection and quantification with FragPipe:

| Input  | Type | Required? | Description |
|-|-|-|-|
| Mass spectrometry data | mzML, RAW | Yes | Mass spectrometry files. |
| Proteome database (proteome_database.fasta) | FASTA | Yes | Generated in module 1c. Amino acid sequences of all ORFs in the isoforms found to be expressed by the scRNA-seq samples. |

The FragPipe output file generated should be either `peptides.txt` or `peptide.tsv`. Additionally, if FragPipe were configured to perform peptide quantification, `report.pr_matrix.tsv` would also be output.

#### Applying variance stabilizing normalization to peptide counts:

| Input  | Type | Required? | Description |
|-|-|-|-|
| Peptide counts | CSV / TXT / TSV | Yes | Generated from peptide quantification with FragPipe. |

| Output | Type | Description |
|-|-|-|
| Normalized peptide counts | CSV | Peptide counts normalized with the variance stabilizing normalization algorithm. |

### 3. Integration

#### Integration of transcriptomics and proteomics data:

| Input  | Type | Required? | Description |
|-|-|-|-|
| Proteomics peptide data | CSV / TSV / TXT | Yes | Generated in module 2. Peptide results. Typically 'peptides.txt', 'peptide.tsv' or 'report.pr_matrix.tsv'. |
| Proteome database (proteome_database.fasta) | FASTA | Yes | Generated in module 1c. Amino acid sequences of all ORFs in the isoforms found to be expressed by the scRNA-seq samples. |
| Proteome database metadata (proteome_database_metadata.txt) | TXT | Yes | Generated in module 1c. Information on each ORF in the database. |
| Proteome database transcripts (proteome_database_transcripts.gtf) | GTF | Yes | Generated in module 1c. Annotation of all transcripts used to generate the database. |
| Seurat R object | RDS | Yes | Generated in module 1b. Seurat object of an scRNA-seq sample. Can be loaded into R with the `readRDS()` function. |

| Output | Type | Filename | Description |
|-|-|-|-|
| Peptide information | CSV | peptide_info.csv | Main results file with peptide mapping data. |
| Report | HTML | summary_report.html | Summary report. |
| Combined annotation data | GTF | combined_annotations.gtf | Annotations of transcripts, peptides, and ORFs for visualization. |
| Transcripts with CDS annotation data | GTF | transcripts_and_ORFs_for_isovis.gtf | Annotations of transcripts and ORFs for visualization in IsoVis. |
| Transcript coordinates | BED12 | transcripts.bed12 | Transcript spliced genomic coordinates. |
| Peptide coordinates | BED12 | peptides.bed12 | Peptide spliced genomic coordinates. |
| ORF coordinates | BED12 | ORFs.bed12 | ORF spliced genomic coordinates. |
| Transcript expression information | CSV | transcript_expression_info.csv | The level of evidence suggesting translation of each transcript and the cell clusters found to be expressing each transcript. |

#### Integration of combined annotations and marker genes:

| Input | Type | Required? | Description |
|-|-|-|-|
| Combined annotation data | GTF | Yes | Generated from the integration of transcriptomics and proteomics data in module 3. |
| Marker gene information | CSV | Yes | Generated in module 1b. |

| Output | Type | Filename | Description |
|-|-|-|-|
| Combined annotation data with marker genes | GTF (zipped) | combined_annotations_with_marker_genes.gtf | Combined annotation data containing marker gene information. |

#### Description of `peptides_info.csv` output:

| Column name | Class | Description |
|-|-|-|
| peptide | character | Peptide sequence |
| accession | character | Protein accession |
| PID | character | protein_accession\|CO=genomic_coordinates (included for compatibility with FASTA header) |
| transcript_id | character | ENSEMBL or novel transcript ID |
| gene_id | character | ENSEMBL gene ID |
| gene_name | character | Gene name/symbol |
| strand | character | Strand (+ or -) |
| number_exons | integer | Number of exons spanned by the peptide |
| transcript_length | integer | Transcript length (nt) |
| transcript_biotype | character | Transcript biotype from GTF |
| simplified_biotype | character | Simplified transcript biotype |
| protein_length | integer | Protein length (AA) |
| orf_genomic_coordinates | numeric | Genomic coordinates of ORF |
| orf_type | character | Annotated, unannotated, or variant protein |
| localisation | character | ORF location in the genome (see ORF category definitions) |
| uniprot_status | character | Review status in UniProt (reviewed / unreviewed) |
| openprot_id | character | OpenProt ID (if present) |
| molecular_weight(kDA) | numeric | Molecular weight of protein (kDa) |
| isoelectric_point | numeric | Isoelectric point of ORF calculated using pKa scale EMBOSS (Rice et al., 2000) |
| hydrophobicity | numeric | Hydrophobicity profile of ORF calculated using Kyte-Doolittle scale (Kyte et al., 1982) |
| aliphatic_index | numeric | Aliphatic index of ORF (Ikai 1980) |
| longest_orf_in_transcript | true/false | Longest ORF in the transcript (longest within CDS regions for known proteins) |
| peptide_ids_gene | true/false | Is peptide uniquely mapped to gene? |
| peptide_ids_orf | true/false | Is peptide uniquely mapped to ORF? |
| peptide_ids_transcript | true/false | Is peptide uniquely mapped to transcript? |
| shared_novel_protein_peptide | true/false | Is peptide shared with other novel proteins? |
| orf_identified | true/false | Is ORF identified with unique peptide evidence? |
| gene_identified | true/false | Is gene identified with unique peptide evidence? |
| transcript_identified | true/false | Is transcript identified with unique peptide evidence? |

#### Description of `transcript_expression_info.csv` output:

| Column name | Class | Description |
|-|-|-|
| transcript_id | character | ENSEMBL or novel transcript ID |
| expressed_in_clusters | character | List of cell clusters expressing the transcript |
| translation_evidence_level | integer | A number representing the level of evidence suggesting the transcript was translated. 0, 1 and 2 refer to a lack of evidence, potential evidence, and definite evidence respectively. |

### 4. Visualization

| Input (in GenomeProtSC) | Type | Required? | Description |
|-|-|-|-|
| Transcript expression information | CSV | No | Generated in module 3. The level of evidence suggesting translation of each transcript and the cell clusters found to be expressing each transcript. |
| Seurat R object | RDS | Only if the above file is uploaded | Generated in module 1b. Seurat object of an scRNA-seq sample. Can be loaded into R with the `readRDS()` function. |

| Input (in IsoVis) | Type | Required? | Description |
|-|-|-|-|
| Combined annotations | GTF | Yes | Generated in module 3. Annotations of transcripts, peptides and ORFs for visualization. Can contain marker gene information if the user integrated such information from module 1b. Upload as a **stack data** file in IsoVis. |
| Peptide intensities | CSV / TXT / TSV | No | Generated in module 2 if peptide quantification were done. Typically 'report.pr_matrix.tsv'. Upload as a **peptide counts data** file in IsoVis. |

| Output | Type | Filename | Description |
|-|-|-|-|
| Exported figures (zipped) | PDFs | <plot_name>.gtf | One PDF for each figure drawn in GenomeProtSC. |