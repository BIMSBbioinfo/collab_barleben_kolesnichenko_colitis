# find cell-cell communications using cellchat
library(CellChat)
library(Seurat)

# https://github.com/jinworks/CellChat
# tutorial https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html

args <- commandArgs(trailingOnly = T)

message(date() ," => importing seurat object")
s <- readRDS(args[1]) #read seurat object
print(s)

# simplify cell type labels (group rare ones into "other")
s$celltype <- ifelse(s$celltype %in% names(which(table(s$celltype) < 100)), 
                       'other', s$celltype)
# update group labels 
conditions <- list('PoolA' = 'Recovery', 
                   'PoolB' = 'DSS', 
                   'PoolC' = 'Untreated')
s$condition <- paste(conditions[s$sample])

message(date(), " => Finetuning cell type labels for epithelial and immune cells")
M <- Seurat::GetAssayData(s)

s$celltype2 <- s$celltype
# define T-regulatory cells as those that express either Cd25 (Il2ra) or Foxp3 
tregs <- union(names(which(M['Il2ra',] > 0)), names(which(M['Foxp3',] > 0)))
s@meta.data[tregs,]$celltype2 <- 'Tregs'

# keep immune cells and Cdh1+ expression epithelial cells 
immune_cells <- names(s$celltype2[s$celltype2 %in% c('B cells', 'DC', 'Macrophages', 'Mast cells', 
                                  'NKT', 'T cells', 'Tgd', 'Tregs')])

# remove epithelial cells that don't express Cdh1 marker and express immune marker Ptprc (cd45)
epithelial_cells <- setdiff(intersect(colnames(s)[s$celltype == 'Epithelial cells'], 
                                      names(which(M['Cdh1',] > 0))),
                            names(which(M['Ptprc',] > 0)))
s <- s[,c(epithelial_cells, immune_cells)]

s$samples <- as.factor(s$sample)

# for each condition, compute cell-cell communication separately 
results <- sapply(simplify = F, unique(s$condition), function(x) {
  pool <- colnames(s)[s$condition == x]
  message(date(), " => creating cellchat object for condition ",x)
  cellChat <- createCellChat(object = s[,pool], group.by = "celltype2", assay = "RNA")
  CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  cellChat@DB <- CellChatDB 
  cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database
  message(date(), " => identifying overexpressed genes and interactions")
  cellChat <- identifyOverExpressedGenes(cellChat)
  cellChat <- identifyOverExpressedInteractions(cellChat)
  
  message(date(), " => computing communication probabilities")
  cellChat <- computeCommunProb(cellChat, type = "triMean", nboot = 20)
  cellChat <- filterCommunication(cellChat, min.cells = 10)
  
  message(date(), " => computing communication probabilities at pathway level")
  cellChat <- computeCommunProbPathway(cellChat)
  
  message(date(), " => aggregating networks")
  cellChat <- aggregateNet(cellChat)
  return(cellChat)
})

# use Secreted Signaling
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
# cellChat@DB <- CellChatDB.use

message(date(), " => saving cellchat objects")
saveRDS(results, file = 'cellChat_by_condition.RDS')

message(date(), " => finished cellchat analysis")

