# Description

This repo contains data and code used to analyse RNA-seq and single-cell RNA-seq datasets 
in collaboration with Luisa Barleben (Marina Kolesnichenko's Lab). 

# Publications

The latest public version of the manuscript can be found [here](https://www.researchsquare.com/article/rs-5431057/v1). 

Short summary of the manuscript: 
"Epithelial NF-κB signaling plays a paradoxically protective role in colitis. Using ulcerative colitis (UC) patient biopsies and mouse models, we show that NF-κB-driven pro-inflammatory signals from intestinal epithelial cells recruit regulatory T cells and aid recovery. In UC, NF-κB selectively activates genes that shape immune responses. These findings uncover a key epithelial mechanism essential for resolving inflammation." 

# Dependencies

RNA-seq and scRNA-seq analyses requires the following R packages. You can install them using the commands below.

```
install.packages(c("data.table", "ggplot2", "ggpubr", "pbapply", "ggrepel", "pheatmap", "openxlsx", "effectsize", "BiocManager", "devtools", "Seurat", "gprofiler2", "stringr", "pheatmap", "ggrepel", "ggridges"))
BiocManager::install(c("RUVSeq", "EDASeq", "DESeq2", "progeny", "singscore", 'celldex', 'SingleR', 'AUCell', 'BiocParallel', 'ComplexHeatmap'))
devtools::install_github('immunogenomics/presto') #optional, makes marker finding in Seurat extra fast 
devtools::install_github("jinworks/CellChat") 
```

You will also need `guix` to install other tools (PiGx-RNAseq pipeline), BBMAP aligner, Samtools. See below. 

# Data

## RNA-seq 

Bulk RNAseq of mouse colon crypts under six different conditions was carried out using Illumina TruSeq strand-specific protocol in paired-end mode.
The sample conditions consisted of a combination of two distinct NF-κB states (suppressed or WT) and three different DSS treatments (untreated, colitis, recovery). 
Each condition consisted of 4 replicates, thus a total of 24 samples were obtained. The raw sequences were processed and analysed using the [PiGx RNA-seq pipeline](https://github.com/BIMSBbioinfo/pigx_rnaseq).

  - **Raw Data**: Raw fastq reads can be downloaded from [here](https://bimsbstatic.mdc-berlin.de/akalin/buyar/marina/manuscript_data_colitis/rnaseq_reads.tgz). 
  - **Annotations**: DNA, cDNA, ang GTF files used in the analysis can be downloaded from [here](https://bimsbstatic.mdc-berlin.de/akalin/buyar/marina/manuscript_data_colitis/annotations.tgz). 
  - **Sample Sheet**: Contains the experimental setup. See `./data/rnaseq/sample_sheet.csv`
  - **PiGx RNA-Seq Settings File**: The settings file required to run pigx-rnaseq pipeline. See `./data/rnaseq/settings.yaml` 
  - **Processed count tables**: Raw/normalized feature count tables can be found under `./data/rnaseq/feature_counts`
  
## single-cell RNA-seq

  - **Raw Reads and CellRanger Output**: Raw fastq reads along with CellRanger outputs can be downloaded from [here](https://bimsbstatic.mdc-berlin.de/akalin/buyar/marina/manuscript_data_colitis/scrnaseq_reads.tgz)
  - **Processed Seurat Object**:
      - **RDS format**: https://bimsbstatic.mdc-berlin.de/akalin/buyar/marina/manuscript_data_colitis/seu.RDS
      - **cloupe format**: https://bimsbstatic.mdc-berlin.de/akalin/buyar/marina/manuscript_data_colitis/seu.cloupe
  - **Processed CellChat Object**: https://bimsbstatic.mdc-berlin.de/akalin/buyar/marina/manuscript_data_colitis/cellChat_by_condition.RDS

# Analysis 

## RNA-seq 

1. Raw RNA-seq data was processed using PiGx-RNAseq pipeline. 

```
guix shell pigx-rnaseq
pigx-rnaseq -s ./data/rnaseq/settings.yaml ./data/rnaseq/sample_sheet.csv
```

The processed count tables can be found under `./data/rnaseq/feature_counts`. 

2. The raw RNA-seq counts were further corrected for batch effects using RUV-seq and differential expression analysis was carried out using DESeq2. 

Compile a report exploration of the count tables, batch correction, gene-set analysis, and DE analysis:

```
srcdir=$(pwd)/src
datadir=$(pwd)/data 
outdir=$(pwd)

Rscript -e "library(rmarkdown); rmarkdown::render('./src/analysis.Rmd', \
  params = list(datadir = '${datadir}', srcdir = '${srcdir}', outdir = '${outdir}'), \
  output_file = 'analysis.html', \
  output_dir = '${outdir}')"
```

3. Results: 

- HTML report:  `./results/rnaseq/analysis.html`. 
- The DEseq2 results: `./results/rnaseq/DE_results.xlsx`. 

## scRNA-seq

This experiment has an EGFP-knockin expectedly in 1-3% of the cells. The initial CellRanger output for count tables had discarded a significant portion of potential EGFP containing cells. 
So, we re-aligned the reads against EGFP construct sequence. Then, we extract the cellular barcodes of the reads that contain EGFP hits. 
Then we rescue those cells from the CellRanger outputs and obtained a Seurat object that includes the initial CellRanger processed cells along with candidate EGFP-containing cells.

1. Download and export the contents of the single-cell read files

```
cd ./data/scrnaseq
wget https://bimsbstatic.mdc-berlin.de/akalin/buyar/marina/manuscript_data_colitis/scrnaseq_reads.tgz
tar -xzvf scrnaseq_reads.tgz
```

2. Find EGFP reads: re-align the reads against EGFP construct. 

```
# go to folder and activate environment
cd single_cell_analysis 
guix package --manifest=guix.scm --profile=colitis
source colitis/etc/profile
# assuming you downloaded and extracted the reads archive (see Step 1)
mkdir find_egpf_reads; cd find_egfp_reads
bash ../../src/realign_egfp.sh ../../data/scrnaseq/egfp.fa ../../data/scrnaseq/reads 
```

3. Find EGFP Cells

Based on the reads that contain EGFP fragments (from step 2); find the cell barcodes corresponding to these reads.
EGFP fragments are in _R2 read files; we need to map them to _R1 read files and extract their read ids, seqs, and barcodes.
Then, we can match those barcodes to quantified barcodes from cellranger.

```
source colitis/etc/profile 
cd single_cell_analysis   # unless you are not there already 
mkdir find_egfp_cells; cd find_egfp_cells
# usage: Rscript /path/to/src/subset_fastq_egfp_ids.R  /path/to/scrnaseq-reads /path/to/find_egpf_reads 
Rscript ../../src/subset_fastq_egfp_ids.R ../../data/scrnaseq/reads ../find_egfp_reads 
```

4. Map the EGFP cell barcodes to all barcodes quantified by cellranger
   
We combine the initial CellRanger count tables with the EGFP cell counts, and create a Seurat object. 
Then we processes it, annotate cell labels, assign gene-set scores using AUCell and prints the object seu.RDS. 

```
cd single_cell_analysis  # unless you are not there already 
# Usage: Rscript /path/to/src/process_seurat_with_egfp_cells.R /path/to/data/scrnaseq/reads /path/to/output/from/step3 /path/to/data/genesets/cancersea.mouse.gmt
Rscript ../src/process_seurat_with_egfp_cells.R ../data/scrnaseq/reads find_egfp_cells ../data/genesets/cancersea.mouse.gmt
```

5. Run CellChat to detect cell-cell communications

```
cd single_cell_analysis   # unless you are not there already 
mkdir find_cell_communication; cd find_cell_communication
Rscript ../src/run_cellchat.R seu.RDS #output from Step 4
```

6. Prepare report for EGFP analysis 

Compile the rmarkdown script at `src/manuscript_figures_egfp_analysis.Rmd`

```
cd single_cell_analysis
srcdir=`readlink -f ../src` 
seurat_file=`readlink -f seu.RDS`
datadir=`readlink -f ../data`
outdir=$(pwd)

Rscript -e "library(rmarkdown); rmarkdown::render('../src/manuscript_figures_egfp_analysis.Rmd', \
  params = list(datadir = '${datadir}', srcdir = '${srcdir}', seurat_file = '${seurat_file}'), \
  output_file = 'manuscript_figures_egfp_analysis.html', \
  output_dir = '${outdir}')"
```


7. Prepare report for CellChat analysis 

```
cd single_cell_analysis
cellchat_file=`readlink -f cellChat_by_condition.RDS`
outdir=$(pwd)

Rscript -e "library(rmarkdown); rmarkdown::render('../src/manuscript_figures.cellChat.Rmd', \
  params = list(cellchat_file = '${cellchat_file}'), \
  output_file = 'manuscript_figures.cellChat.html', \
  output_dir = '${outdir}')"
```


8. Results

- EGFP analysis report: `./results/manuscript_figures_egfp_analysis.html`. 
- CellChat analysis report: `./results/manuscript_figures.cellChat.html`. 










  
