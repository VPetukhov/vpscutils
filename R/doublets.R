#' @importFrom magrittr %<>% %$% %>%
NULL

#' @export
PredictGeneLimits <- function(cm, verbose=T, min.quant=0.15, max.quant=0.85) {
  cm@x <- (cm@x > 0) * 1
  n.genes <- Matrix::colSums(cm)
  quants <- seq(min.quant, max.quant, length.out=50)
  q.vals <- quantile(n.genes, quants)
  model <- lm(q.vals ~ poly(quants, 3))
  if (verbose) {
    message("Max abs error: ", round(max(abs(residuals(model))), 3), "; max rel error: ", round(max(abs(residuals(model) / q.vals)), 7))
  }

  return(predict(model, data.frame(quants=c(0, 1))))
}

#' @export
GetScrubletScores <- function(mat, min.molecules.per.gene=10) {
  tf.in <- tempfile()
  tf.out <- tempfile()
  dt <- mat[Matrix::rowSums(mat)>=min.molecules.per.gene,] %>% Matrix::t() %>% as.matrix() %>% data.table::data.table()
  data.table::fwrite(dt, file=tf.in)

  cmd <- paste0("/d0-mendel/home/viktor_petukhov/local/anaconda3/bin/python3.7 -c 'import sys; import pandas; import scrublet; ",
                "df = pandas.read_csv(\"", tf.in, "\"); scrub = scrublet.Scrublet(df); ",
                "doublet_scores, predicted_doublets = scrub.scrub_doublets();",
                "pandas.DataFrame(dict(score=doublet_scores, is_doublet=predicted_doublets)).to_csv(\"",tf.out,"\");'",
                sep='')
  system(cmd, intern=F)
  x <- data.table::fread(tf.out,sep=',')[,2:3] %>% as.data.frame() %>% as.list() %>%
    lapply(`names<-`, colnames(mat))

  return(x)
}

#' findBestClusterLinComb
#' @description compare all pair of clusters to explain given `clust` by their linear combination.
FindBestClusterLinComb <- function(clust, counts, clusters, plot=T) {
  clust.means <- conos:::collapseCellsByType(counts, clusters)
  clust.combs <- clusters %>% as.factor() %>% levels() %>% .[. != clust] %>% combn(2)

  lms <- pbapply::pbapply(clust.combs, 2, function(x) {
    lm(clust.means[clust, ] ~ clust.means[x[1], ] + clust.means[x[2], ] - 1 )
  })

  comb.id <- sapply(lms, function(lm) var(lm$residuals) / var(clust.means[clust, ])) %>% which.min()
  best.lm <- lms[[comb.id]]
  cl1 <- clust.combs[1, comb.id]
  cl2 <- clust.combs[2, comb.id]

  if (plot) {
    ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x=best.lm$coefficients[1] * clust.means[cl1, ] + best.lm$coefficients[2] * clust.means[cl2, ], y=clust.means[clust,]), size=0.1) +
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
      ggplot2::labs(x="lin comb", y="observed expression")
  }

  return(list(cl1=cl1, cl2=cl2, lm=best.lm))
}
