# Given EGFP reads and their corresponding cell barcodes; map the barcodes to quantified cells from cellranger
# then explore the single-cell data for different values of mismatches allowed when mapping egfp barcodes 
library(Seurat)
library(AUCell)
library(SingleR)
library(celldex)
library(data.table)
library(BiocParallel)
library(ggpubr)
library(ggplot2)
data.table::setDTthreads(4)
ggplot2::theme_set(ggpubr::theme_pubclean())


args <- commandArgs(trailingOnly=T)

readsDir <- args[1] #data/scrnaseq/reads 
egfp_cells_dir <- args[2] #find_egfp_cells
geneset_file <- args[3] #data/genesets/cancersea.mouse.gmt                      
outdir <- getwd()
cellranger_outdir <- file.path(readsDir, 'cellranger_output')


# given a seurat object, using singleR to annotate cells 
annotate_cells <- function(seu) {
  mouse <- celldex::ImmGenData() #celldex::MouseRNAseqData()
  require(SingleR)
  main_labels <- SingleR(test = Seurat::GetAssayData(seu), 
                  ref = mouse, labels = mouse$label.main, 
                  BPPARAM=BiocParallel::MulticoreParam(24))
  fine_labels <- SingleR(test = Seurat::GetAssayData(seu), 
                         ref = mouse, labels = mouse$label.fine, 
                         BPPARAM=BiocParallel::MulticoreParam(24))
  return(list('celltype' = main_labels$pruned.labels, 'cell_subtype' = fine_labels$pruned.labels))
}

process_seurat <- function(exp, ...) {
  # see https://satijalab.org/seurat/articles/pbmc3k_tutorial
  # process seurat object
  seu <- Seurat::CreateSeuratObject(counts = exp, ...)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seu)
  seu <- ScaleData(seu, features = all.genes)
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  seu <- FindNeighbors(seu, dims = 1:10)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:10)
  return(seu)
}
# score cells for cancersea state signatures 
import_gmt <- function(f) {
  sets <- list()
  lines <- readLines(f)
  for (l in lines) {
    w <- strsplit(l, "\t")[[1]]
    sets[[w[1]]] <- w[-c(1,2)]
  }
  return(sets)
}


score_gene_sets <- function(gene_set_path, seu) {
  genesets <- import_gmt(gene_set_path)
  cells_AUC <- AUCell_run(Seurat::GetAssayData(seu), genesets, 
                          BPPARAM = BiocParallel::MulticoreParam(24))
  scores <- data.table(cells_AUC@assays@data$AUC, keep.rownames = T)
  scores <- melt.data.table(scores, id.vars = 'rn')
  dt <- data.table('variable' = colnames(seu), 'sample'  = seu$sample, 
                   'egfp_status' = seu$egfp_status, 'celltype' = seu$celltype)
  scores <- merge.data.table(scores, dt, by = 'variable')
  return(scores)
}

get_basis_matrix <- function(M, factors) {
  B <- pbapply::pbapply(M, 2, function(x) {
    tapply(X = x, factors, FUN = mean)
  })
  return(B)
}

# import cellranger output 

# get expression data for the samples 
expL <- sapply(simplify = F, c('PoolA', 'PoolB', 'PoolC'), function(x) {
  message(date(), " => importing raw counts for ",x)
  # import unfiltered expression matrix; find egfp cells keep them. 
  exp <- Seurat::Read10X(file.path(cellranger_outdir, x, 'outs', 'raw_feature_bc_matrix'))
  # find cells with egfp; save them on top of filtered cells
  egfp_cells <- data.table::fread(file.path(egfp_cells_dir, 
                                            paste0('egfp_reads2cells_', x, '.tsv')))
  barcodes <- unique(paste0(egfp_cells$barcode, "-1"))
  exp_egfp <- exp[,intersect(barcodes, colnames(exp))]
  
  message(date(), " => Found ",ncol(exp_egfp), " egfp cells in unfiltered matrix for ",x)
  
  # import filtered expression data from cellranger output 
  message(date(), " => Importing filtered count matrix")
  exp_filtered <- Seurat::Read10X(file.path(cellranger_outdir, x, 'outs', 'filtered_feature_bc_matrix'))
  exp_filtered <- exp_filtered[,setdiff(colnames(exp_filtered), colnames(exp_egfp))]
  
  # merge egfp cells with the filtered cells 
  message(date(), " => Combining egfp cells with filtered count matrix")
  exp <- cbind(exp_egfp, exp_filtered)
  return(list('exp' = exp, 'egfp_cells' = colnames(exp_egfp)))
})

# combine exp matrices
features <- Reduce(intersect, lapply(expL, function(x) rownames(x$exp)))
counts <- do.call(cbind, lapply(names(expL), function(x) { 
  M <- expL[[x]][['exp']][features,] 
  colnames(M) <- paste0(x, "_", colnames(M))
  return(M)
  }))
samples <- do.call(c, lapply(names(expL), function(x) rep(x, ncol(expL[[x]][['exp']]))))
egfp_cells <- paste(do.call(c, lapply(names(expL), function(x) paste0(x, "_", expL[[x]]$egfp_cells))))

# add "EGFP" as a feature also to the count matrix.
# count: 1 if egfp is available
egfp_counts <- ifelse(colnames(counts) %in% egfp_cells, 1, 0)
counts <- rbind(counts, matrix(egfp_counts, nrow = 1, dimnames = list('rownames' = 'Egfp')))  

seu <- process_seurat(counts, min.cells = 3, min.features = 50)
seu$sample <- gsub("_.+$", "", colnames(seu))
seu$egfp_status <- ifelse(colnames(seu) %in% egfp_cells, 'EGFP_POS', 'EGFP_NEG')
message(date(), " => breakdown of egfp status in final object" )
print(table(seu$sample, seu$egfp_status))

# annotate cell types 
labels <- annotate_cells(seu)
seu$celltype <- labels$celltype
seu$cell_subtype <- labels$cell_subtype

# cancersea state signature scores 
scores <- score_gene_sets(geneset_file, seu)

# add scores to the seurat object
m <- dcast.data.table(scores, variable ~ rn, value.var = 'value')
m <- as.matrix(data.frame(m[,-1], row.names = m[[1]], check.names = F))
# reorder rows to match seu samples
m <- m[colnames(seu),]
seu@meta.data <- cbind(seu@meta.data, m)

# save the processed seurat object
saveRDS(seu, file = file.path(outdir, 'seu.RDS'))

message(date(), " => Finished processing Seurat object")



