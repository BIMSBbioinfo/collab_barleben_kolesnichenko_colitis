# import single-cell data from cellranger output 
# import reads that were found to contain egfp fragments; 
# find their mates 

library(ShortRead)
library(data.table)
library(Biostrings)

args <- commandArgs(trailingOnly = T)
readsDir <- args[1]
egfp_dir <- args[2]  

subset_fastq_by_ids <- function(fq, ids) {
  y <- 1000000
  f <- ShortRead::FastqStreamer(fq, n = y)
  base <- 0
  matches <- c()
  while (length(sub <- ShortRead::yield(f))) {
    m <- match(ids, paste(ShortRead::id(sub)))
    m <- setdiff(m, NA)
    if(length(m) > 0) {
      matches <- c(matches, sub[m])
    }
    base <- base + y
    if(sum(lengths(matches)) == length(ids)) {
      cat("Found all matches\n")
      print(matches)
      break
    }
    print(base)
  }
  close(f)
  return(matches)
}

process_matches <- function(matches) {
  matches <- do.call(rbind, lapply(matches, function(x) {
    data.table::data.table('id' = paste(ShortRead::id(x)), 
                           'read' = paste(ShortRead::sread(x)), 
                           'barcode' = paste(Biostrings::subseq(ShortRead::sread(x), 1, 16)))
  }))
  return(matches)
}


get_egfp_ids <- function(f) {
  ids <- unlist(lapply(strsplit(readLines(f)[-c(1:3)], "\t"), function(x) x[1]))
  # get egfp read ids; convert to R1 
  ids <- gsub(" 2:N", " 1:N", ids)
  return(ids)
}

reads <- list('PoolA' = list('fq' = file.path(readsDir, 'P3531_SP255_003_PoolA_S3_L001_R1_001.fastq.gz'), 
                            'egfp_reads' = file.path(egfp_dir, 'outA.sam')), 
              'PoolB' = list('fq' = file.path(readsDir, 'P3531_SP255_006_PoolB_S4_L001_R1_001.fastq.gz'), 
                            'egfp_reads' = file.path(egfp_dir, 'outB.sam')), 
              'PoolC' = list('fq' = file.path(readsDir, 'P3531_SP255_009_PoolC_S5_L001_R1_001.fastq.gz'),
                            'egfp_reads' =  file.path(egfp_dir, 'outC.sam')))

M <- lapply(names(reads), function(pool) {
  cat("processing, ",reads[[pool]]$fq,"\n")
  ids <- get_egfp_ids(reads[[pool]]$egfp_reads)
  matches <- subset_fastq_by_ids(reads[[pool]]$fq, ids)
  egfp_reads2cells <- process_matches(matches)
  write.table(egfp_reads2cells, file = paste0('egfp_reads2cells_', pool, ".tsv"), quote = F, row.names = F,
              sep = '\t')
  return(matches)
})

message("Finished processing reads")
