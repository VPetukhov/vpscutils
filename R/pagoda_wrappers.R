#' @export
GetPagoda <- function (cm, n.cores = 30, clustering.type = "leiden", embeding.type = "tSNE", verbose=TRUE,
                       n.pcs=100, distance="cosine", trim=5, n.odgenes=1000, graph.k=30,
                       od.genes=NULL, clustering.resolution=2, build.graph=TRUE, var.scale=TRUE, min.dist=0.5, spread=1.5, ...) {
  r <- pagoda2::Pagoda2$new(cm, trim=trim, n.cores=n.cores, verbose=verbose, ...)

  if (var.scale) {
    r$adjustVariance(plot = F, do.par = F, gam.k = 10, verbose = verbose)
  }

  r$calculatePcaReduction(nPcs = n.pcs, n.odgenes = n.odgenes, odgenes=od.genes, maxit = 1000, verbose=verbose, var.scale=var.scale)

  if (!build.graph)
    return(r)

  r$makeKnnGraph(k = graph.k, type = "PCA", center = T, distance = distance,
                 weight.type = "none", verbose = verbose)

  for (ct in clustering.type) {
    switch (ct,
            infomap = r$getKnnClusters(method = igraph::infomap.community, type = "PCA", name = "infomap"),
            multilevel = r$getKnnClusters(method = igraph::multilevel.community, type = "PCA", name = "multilevel"),
            leiden = r$getKnnClusters(method = conos::leiden.community, type = "PCA", name = "leiden", resolution=clustering.resolution),
            stop("Unknown clustering type: ", ct)
    )
  }

  for (et in embeding.type) {
    r$getEmbedding(type = "PCA", embeddingType = et, distance=distance, min_dist=min.dist, spread=spread)
  }

  return(r)
}

#' @export
GetPagodaWebApp <- function(p2, clusters, organism=NULL, additional.metadata=list(), verbose=T, go.sets=NULL, go.env=NULL, test.pathways=TRUE) {
  if (is.null(go.env)) {
    if (is.null(organism))
      stop("Either organism or go.env must be provided")

    if (verbose) cat("Generate go environment\n")
    go.env <- pagoda2::p2.generate.go(p2, organism=organism)
  }

  if (is.null(go.sets)) {
    if (verbose) cat("Generate genesets\n")

    go.sets <- ExtractGoSets(go.env, verbose=verbose)
  }

  if (verbose) cat("Generate de geneset\n")
  de.sets <- pagoda2::get.de.geneset(p2, groups = clusters, prefix = 'de_')

  if (test.pathways) {
    if (verbose) cat("Test pathway overdispersion\n")
    p2$testPathwayOverdispersion(setenv = go.env, verbose = verbose, correlation.distance.threshold = 0.8,
                                 recalculate.pca = F, min.pathway.size = 50, max.pathway.size = 1000)
  }

  if (verbose) cat("Create app\n")
  p2.web <- p2 %>%
    pagoda2::make.p2.app(
      dendrogramCellGroups = as.factor(clusters),
      geneSets = c(go.sets, de.sets),
      additionalMetadata=additional.metadata,
      show.clusters = T);

  if (verbose) cat("All done!\n")

  return(p2.web)
}

#' @export
ExtractGoSets <- function(go.env, verbose=F) {
  names(go.env) %>% setNames(., .) %>% sccore:::plapply(function(x)
    list(properties = list(locked=T, genesetname=x, shortdescription=GO.db::GOTERM[[x]]@Term),
         genes = c(go.env[[x]])), progress=verbose)
}

#' @export
Pagoda2FromConos <- function(con, embedding.type="tSNE", annotation=NULL, ...) {
  p2 <- con$getJointCountMatrix(raw=T) %>% Matrix::t() %>% GetPagoda(build.graph=F, ...)
  p2$graphs$conos <- con$graph

  if (!is.null(con$embedding)) {
    p2$embeddings$PCA[[embedding.type]] <- con$embedding
  }

  p2$clusters %<>% c(con$clusters %>% setNames(paste0("conos_", names(.))) %>%
                       lapply(`[[`, "groups") %>% lapply(as.factor))

  p2$clusters$dataset <- con$getDatasetPerCell()

  if (!is.null(annotation)) {
    p2$clusters$annotation <- as.factor(annotation[names(p2$clusters$dataset)])
  }

  return(p2)
}

#' @export
ConvertMetadataToPagoda2Format <- function(...) {
  metadata.list <- list(...)# lapply(list(...), as.factor)
  metadata.list <- metadata.list[!sapply(metadata.list, is.null)] %>% lapply(as.factor)
  return(mapply(pagoda2::p2.metadata.from.factor, metadata.list, names(metadata.list), SIMPLIFY=F))
}

#' @export
GetConos <- function (p2s, annotations = NULL, k = 30, k.self = 10, cluster.resolution=2,
                      k.self.weight = 0.5, embedding.method = "UMAP", matching.mask = NULL,
                      n.cores=1, space = "PCA", min.dist=0.5, spread=1.5, ...)
{
  con <- conos::Conos$new(p2s, n.cores = n.cores)
  con$buildGraph(base.groups=annotations, ...)
  con$findCommunities(method=leidenAlg::leiden.community, resolution=cluster.resolution)
  con$embedGraph(method=embedding.method, verbose=FALSE, min.dist=min.dist, spread=spread)
  return(con)
}
