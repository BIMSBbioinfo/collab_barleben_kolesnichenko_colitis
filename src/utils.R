# get a basis matrix that is the mean value per factor of variable
get_basis_matrix <- function(M, factors) {
  B <- pbapply::pbapply(M, 2, function(x) {
    tapply(X = x, factors, FUN = mean)
  })
  return(B)
}

compute_PCA <- function(exp, topN = 50) {
  top <- head(names(sort(apply(exp, 1, sd), decreasing = T)), 5000)
  M <- t(exp[top,])
  pca <- stats::prcomp(M, rank. = topN)
  return(pca)
}

get_pca <- function(exp, colData, topN = 50) {
  top <- head(names(sort(apply(exp, 1, sd), decreasing = T)), 5000)
  M <- t(exp[top,])
  fit <- stats::prcomp(M, rank. = topN)
  df <- as.data.frame(fit[['x']])
  df <- merge(df, colData, by = 'row.names')
  var_exp <- round(diag(cov(fit$x))/sum(diag(cov(fit$x))) * 100, 1)
  var_exp <- sapply(names(var_exp), function(x) {
    paste0(x, " (", var_exp[[x]], "%)")
  })
  return(list('dat' = df, 'var_exp' = var_exp))
}

read_msigdb <- function(f) {
  geneSets <- lapply(readLines(f), function(x) {
    unlist(strsplit(x, "\t"))
  })
  #add names
  names(geneSets) <- unlist(lapply(geneSets, function(x) x[1]))
  #remove first two
  geneSets <- lapply(geneSets, function(x) x[-c(1:2)])
  return(geneSets)
}


# given a gene id, count table, factors, make a box plot 
plot_gex <- function(gene, M, factors) {
  df <- data.frame(t(M[gene,,drop=F]), group = factors) 
  colnames(df)[1] <- 'expression' 
  ggpubr::ggboxplot(df, x = 'group',y = 'expression', fill = 'group', add = 'jitter') +
    stat_compare_means()
}

score_gene_set <- function(rankData, genes_up, gene_down = NULL) {
  scores <- singscore::simpleScore(rankData, upSet = genes_up)
  return(scores$TotalScore)
}

compute_anova <- function(M, factors, return_val = 'pvalue') {
  factors <- as.factor(factors)
  if(length(levels(factors)) > 1) {
    dt <- do.call(rbind, pbapply::pblapply(colnames(M), function(x) {
      val <- M[,x]
      df <- data.frame("value" = as.numeric(val), "labels" = factors)
      rownames(df) <- names(val)
      #linear model 
      mod <- lm(value ~ labels, data = df)
      t <- anova(mod)
      effect_size <- effectsize::eta_squared(t)#[['etasq']]
      #stats::TukeyHSD(t)
      return(data.table::data.table("pval" = t$`Pr(>F)`[1], "effect_size"  = effect_size, "variable" = x))
    }))
    dt$padj <- p.adjust(dt$pval, method = 'BH')
    return(dt)
  } else {
    return(NULL)
  }
}
