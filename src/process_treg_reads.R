# combine cellranger output into a single object, process the seurat object, annotate cell types 

# then explore the single-cell data for different values of mismatches allowed when mapping egfp barcodes 
library(Seurat)
library(celldex)
library(data.table)
library(BiocParallel)
library(ggpubr)
library(ggplot2)
library(SingleR)
library(openxlsx)
data.table::setDTthreads(4)
ggplot2::theme_set(ggpubr::theme_pubclean())

args <- commandArgs(trailingOnly=T)

cellranger_outdir <- args[1]
samples <- c('A5391_SP255_003_Wt_colon', 'A5391_SP255_006_ko_colon', 'A5391_SP255_009_Wt_spleen',  'A5391_SP255_012_ko_spleen')

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

get_basis_matrix <- function(M, factors) {
  B <- pbapply::pbapply(M, 2, function(x) {
    tapply(X = x, factors, FUN = mean)
  })
  return(B)
}

# get expression data for the samples 
expL <- sapply(simplify = F, samples, function(x) {
  message(date(), " => importing raw counts for ",x)
  # import filtered expression matrix 
  exp <- Seurat::Read10X(file.path(cellranger_outdir, x, 'filtered_feature_bc_matrix'))
  return(exp)
})

features <- Reduce(intersect, lapply(expL,rownames))
counts <- do.call(cbind, lapply(names(expL), function(x) { 
  M <- expL[[x]][features,] 
  colnames(M) <- paste0(x, "_", colnames(M))
  return(M)
}))
# create seurat object and process 
seu <- process_seurat(counts, min.cells = 3, min.features = 200)
seu$sample <- toupper(sub("A5391_SP255_[0-9]{3}_(.+)_.+$", "\\1", colnames(seu)))
seu$sample <- factor(seu$sample,
                     levels = c("KO_COLON", "WT_COLON", "KO_SPLEEN", "WT_SPLEEN"))
seu$tissue <- gsub("(KO|WT)_", "", seu$sample)
seu$condition <- gsub("_.+$", "", seu$sample)

labels <- annotate_cells(seu)
seu$celltype <- labels$celltype
seu$cell_subtype <- labels$cell_subtype
# remove rare annotations 
keep <- which(seu$celltype %in% names(which(table(seu$celltype) > 50)))
seu <- seu[,keep]

M <- Seurat::GetAssayData(seu)
# define T-regulatory cells as those that express Foxp3 
tregs <- names(which(M['Foxp3',] > 0)) 
# confine the tregs to those that are in clusters 1/3/5/7 
# change this criteria when re-doing the analysis 
tregs <- intersect(tregs, colnames(M)[which(seu$seurat_clusters %in% c(1,3,5,7))])
seu@meta.data[tregs,]$celltype <- 'Tregs'

saveRDS(object = seu, file = 'seu.RDS')

require(loupeR)
seu_cloupe <- loupeR::create_loupe_from_seurat(seu)
file.rename('converted.cloupe', 'seu.cloupe')

# save cell type/subtype percentages to table
m <- as.matrix(table(seu$celltype, seu$sample))
m <- m[order(rowSums(m), decreasing = T),]
m <- apply(m, 2, function(x) x)
m_pc <- apply(m, 2, function(x) round(x/sum(x) * 100, 2))

OUT <- createWorkbook() # for complete tables
addWorksheet(OUT, "Cell_Type_Counts")
writeData(OUT, sheet = "Cell_Type_Counts", x = data.table(m, keep.rownames = T))
addWorksheet(OUT, "Cell_Type_Percentages")
writeData(OUT, sheet = "Cell_Type_Percentages", x = data.table(m_pc, keep.rownames = T))

# save cell type/subtype percentages to table
m <- as.matrix(table(seu$cell_subtype, seu$sample))
m <- m[order(rowSums(m), decreasing = T),]
m <- apply(m, 2, function(x) x)
m_pc <- apply(m, 2, function(x) round(x/sum(x) * 100, 2))

addWorksheet(OUT, "Cell_Subtype_Counts")
writeData(OUT, sheet = "Cell_Subtype_Counts", x = data.table(m, keep.rownames = T))
addWorksheet(OUT, "Cell_Subtype_Percentages")
writeData(OUT, sheet = "Cell_Subtype_Percentages", x = data.table(m_pc, keep.rownames = T))
saveWorkbook(OUT, file.path("Cell_Types_Counts.xlsx"), overwrite = T)





