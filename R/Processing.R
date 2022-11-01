#' Single-Cell Processing
#'
#' @import ArchR
#'
#' @param project ArchR project
#' @param force=T, if TRUE overwrite, ELSE add a {name}_new
#'
#' @rdname Cellcycling
#' @export
Cellcycling <- function(project, force=T){

  all.genes <- project@cellMetadata$genes
  assertthat::assert_that(length(all.genes) > 0)

  s.features <- paste0("^", Seurat::cc.genes$s.genes, "$", collapse = "|")
  s.features <- all.genes[grepl(s.features, all.genes, ignore.case = TRUE)]

  g2m.features <- paste0("^", Seurat::cc.genes$g2m.genes, "$", collapse = "|")
  g2m.features <- all.genes[grepl(g2m.features, all.genes, ignore.case = TRUE)]

  object <- PartialSeurat(project, useMatrix="GeneExpressionMatrix", assay="counts", meta.data=F)

  name <- 'Cell.Cycle'
  features <- list('S.Score' = s.features, 'G2M.Score' = g2m.features)
  ctrl <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  object.cc <- Seurat::AddModuleScore(
    object = object,
    features = features,
    name = name,
    ctrl = ctrl
  )
  cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), value = TRUE)
  cc.scores <- object.cc[[cc.columns]]
  rm(object.cc)
  rm(object)
  gc()
  assignments <- apply(
    X = cc.scores,
    MARGIN = 1,
    FUN = function(scores, first = 'S', second = 'G2M', null = 'G1') {
      if (all(scores < 0)) {
        return(null)
      } else {
        if (length(which(x = scores == max(scores))) > 1) {
          return('Undecided')
        } else {
          return(c(first, second)[which(x = scores == max(scores))])
        }
      }
    }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
  colnames(x = cc.scores) <- c('rownames', 'S.Score', 'G2M.Score', 'Phase')
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c('S.Score', 'G2M.Score', 'Phase')]
  cc.scores$G1.Score = 1 - cc.scores$S.Score - cc.scores$G2M.Score
  assertthat::assert_that(all(rownames(cc.scores) == rownames(project@cellColData)))

  if(!force){
      colnames(cc.scores) <- sapply(colnames(cc.scores), function(x)
                                   ifelse(x %in% colnames(project@cellColData), paste0(x, "_new"),  x)
      )
  }
  project@cellColData <- c(project@cellColData, DataFrame(cc.scores))
  gc()
  return(project)
}

