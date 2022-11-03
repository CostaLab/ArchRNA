
## Create a small seuratobject from ArchRProject
##TODO:
## select features from count matrix
## select assay during the reading of getMatrixFromProject
PartialSeurat <- function(project,
                          reducedDims=NULL,
                          embeddings=NULL,
                          useMatrix="GeneExpressionMatrix",
                          assay="data",
                          features=NULL,
                          meta.data=T){

  includes = c()
  if(!is.null(useMatrix)){
    assertthat::assert_that(useMatrix %in% ArchR::getAvailableMatrices(project))
    seMtx <- ArchR::getMatrixFromProject(project, useMatrix=useMatrix)
    if(is.null(assay)){
      warning(glue::glue("assay is NULL, use first assay {assay}!"))
      assay=names(assays(seMtx))[1]
    }
    assertthat::assert_that(assay %in% names(assays(seMtx)))
    counts <- assays(seMtx)[[assay]]
    genes <- elementMetadata(seMtx)$name
    rownames(counts) <- genes
    includes <- c(includes, "useMatrix")
  }else{
    ## there's no data useful in the Matrix
    ## random a tiny Matrix just for Seurat Creation
    counts <- matrix(rnbinom(nrow(project@cellColData)*10, mu = 4, size = 1), nrow=10 )
    colnames(counts) <- rownames(proj@cellColData)
    rownames(counts) <- letters[1:10]
  }
  if(!is.null(reducedDims)){
    assertthat::assert_that(reducedDims %in% names(project@reducedDims))
    reduction = project@reducedDims[[reducedDims]][[1]]
    includes <- c(includes, "reducedDims")
  }
  if(!is.null(embeddings)){
    assertthat::assert_that(embeddings %in% names(project@embeddings))
    embedd = project@embeddings[[embeddings]][[1]]
    includes <- c(includes, "embeddings")
  }
  meta_data = NULL
  if(meta.data){
    meta_data <- as.data.frame(project@cellColData)
    includes <- c(includes, "meta.data")
  }

  object <- Seurat::CreateSeuratObject(counts=counts,
                                       meta.data=meta_data,
                                       min.cells = 0,
                                       min.features = 0,
                                       verbose=F)

  if("reducedDims" %in% includes){
    object[[reducedDims]] <- Seurat::CreateDimReducObject(reduction, key=glue::glue("{reducedDims}_"))
  }
  if("embeddings" %in% includes){
    object[[embeddings]] <- Seurat::CreateDimReducObject(as.matrix(embedd), key=glue::glue("{embeddings}_"))
  }
  gc()
  return(object)
}

## Create a small seurat assay from ArchRProject
PartialSeuratAssay <- function(project){
  object <- 1
  return(object)
}

## convert Seurat to ArchRProject
Seurat2ArchRProject <- function(object){
  project <- 1
  return(project)
}
