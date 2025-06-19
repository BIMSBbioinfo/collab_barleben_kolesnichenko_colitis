(use-modules (gnu)
             (guix profiles)
             (guix packages)
             (guix utils)
             (guix download)
             (guix build-system gnu)
             (guix licenses)
             ;; Import modules for Python and R packages
             (gnu packages bioinformatics)
)

;; Define the list of packages
(define my-packages
  (list
    ;; Bioinformatics tools
    (specification->package "r-seurat")
    (specification->package "r-shortread")
    (specification->package "r")
    (specification->package "r-data-table")
    (specification->package "r-pbapply")
    (specification->package "r-biostrings")
    (specification->package "r-aucell")
    (specification->package "r-singler")
    (specification->package "r-celldex")
    (specification->package "r-data-table")
    (specification->package "r-biocparallel")
    (specification->package "r-ggpubr")
    (specification->package "r-ggplot2")
    (specification->package "samtools")
    (specification->package "bbmap")
    (specification->package "fastp")
    (specification->package "bedops")
    ))

;; Convert the list of packages to a manifest
(packages->manifest my-packages)

