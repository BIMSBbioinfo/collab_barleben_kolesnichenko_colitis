# Description

This repo contains data and code used to analyse RNA-seq and single-cell RNA-seq datasets 
in collaboration with Luisa Barleben (Marina Kolesnichenko's Lab). 

# Publications

The latest public version of the manuscript can be found [here](https://www.researchsquare.com/article/rs-5431057/v1). 

Short summary of the manuscript: 
"Epithelial NF-κB signaling plays a paradoxically protective role in colitis. Using ulcerative colitis (UC) patient biopsies and mouse models, we show that NF-κB-driven pro-inflammatory signals from intestinal epithelial cells recruit regulatory T cells and aid recovery. In UC, NF-κB selectively activates genes that shape immune responses. These findings uncover a key epithelial mechanism essential for resolving inflammation." 

# Data

## RNA-seq 

Bulk RNAseq of mouse colon crypts under six different conditions was carried out using Illumina TruSeq strand-specific protocol in paired-end mode.
The sample conditions consisted of a combination of two distinct NF-κB states (suppressed or WT) and three different DSS treatments (untreated, colitis, recovery). 
Each condition consisted of 4 replicates, thus a total of 24 samples were obtained. The raw sequences were processed and analysed using the [PiGx RNA-seq pipeline](https://github.com/BIMSBbioinfo/pigx_rnaseq).

  - **Raw Data**: Raw fastq reads along with the DNA, cDNA, ang GTF files used in the analysis can be downloaded from [here](insert link)
  - **Sample Sheet**: Contains the experimental setup. See ./data/rnaseq/sample_sheet.csv
  - **PiGx RNA-Seq Settings File**: The settings file required to run pigx-rnaseq pipeline. See ./data/rnaseq/settings.yaml 
  - **Processed count tables**: Raw/normalized feature count tables can be found under ./data/rnaseq/feature_counts
  
## single-cell RNA-seq

  - **Processed Seurat Object**:
      - **RDS format**: 
      - **cloupe format**: 
      
  - **Processed CellChat Object**:


# Analysis 

## RNA-seq 

1. Raw RNA-seq data was processed using PiGx-RNAseq pipeline. 

```
guix shell pigx-rnaseq
pigx-rnaseq -s ./data/rnaseq/settings.yaml ./data/rnaseq/sample_sheet.csv
```

The processed count tables can be found under `./data/rnaseq/feature_counts`. 

2. The raw RNA-seq counts were further corrected for batch effects using RUV-seq and differential expression analysis was carried out using DESeq2. 

This analysis requires the following R packages. You can install them using the commands below.
```
install.packages(c("data.table", "ggplot2", "ggpubr", "pbapply", "ggrepel", "pheatmap", "openxlsx", "effectsize", "BiocManager"))
BiocManager::install(c("RUVSeq", "EDASeq", "DESeq2", "progeny", "singscore"))
```

Compile a report exploration of the count tables, batch correction, gene-set analysis, and DE analysis:

```
srcdir=$(pwd)/src
datadir=$(pwd)/data 
outdir=$(pwd)

/opt/R/4.5/bin/Rscript -e "library(rmarkdown); rmarkdown::render('./src/analysis.Rmd', \
  params = list(datadir = '${datadir}', srcdir = '${srcdir}', outdir = '${outdir}'), \
  output_file = 'analysis.html', \
  output_dir = '${outdir}')"
```

- HTML report:  `./results/rnaseq/analysis.html`. 
- The DEseq2 results: `./results/rnaseq/DE_results.xlsx' 










  
