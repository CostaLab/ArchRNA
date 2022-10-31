DimPlotS <- function(project, group.by=NULL, reduction=NULL, ...){
  if(reduction %in% names(project@reducedDims)){
    object <- PartialSeurat(project, useMatrix=NULL, reducedDims=reduction)
  }else if(reduction %in% names(project@embeddings)){
    object <- PartialSeurat(project, useMatrix=NULL, embeddings=reduction)
  }else{
    stop("No existing reduction")
  }
  p <- Seurat::DimPlot(object, group.by=group.by, reduction=reduction, ...)
  return(p)
}

DoHeatmapS <- function(project, useMatrix=NULL, assay=NULL, reduction=NULL, ...){
  if(is.null(useMatrix)){
    useMatrix = getAvailableMatrices(project)[1]
    warning(glue::glue("useMatrix is NULL, use first matrix {useMatrix}!"))
  }
  object <- PartialSeurat(project, useMatrix=useMatrix)
  object <- Seurat::ScaleData(object, features=rownames(object))

  p <- Seurat::DoHeatmap(object,  assay="RNA", slot="scale.data", ...)
  return(p)
}

FeaturePlotS <- function(project, useMatrix=NULL, assay=NULL, group.by=NULL, reduction=NULL, ...){
  if(is.null(useMatrix)){
    useMatrix = getAvailableMatrices(project)[1]
    warning(glue::glue("useMatrix is NULL, use first matrix {useMatrix}!"))
  }
  if(reduction %in% names(project@reducedDims)){
    object <- PartialSeurat(project, useMatrix=useMatrix, assay=assay, reducedDims=reduction)
  }else if(reduction %in% names(project@embeddings)){
    object <- PartialSeurat(project, useMatrix=useMatrix, assay=assay, embeddings=reduction)
  }else{
    stop("No existing reduction")
  }
  p <- Seurat::FeaturePlot(object, reduction=reduction, slot=assay, ...)
  return(p)
}
RidgePlotS <-function(project, useMatrix=NULL, assay=NULL, group.by=NULL, ...){
  if(is.null(useMatrix)){
    useMatrix = getAvailableMatrices(project)[1]
    warning(glue::glue("useMatrix is NULL, use first matrix {useMatrix}!"))
  }
  object <- PartialSeurat(project, useMatrix=useMatrix, assay=assay)
  p <- Seurat::RidgePlot(object,assay="RNA", slot=assay, ...)
  return(p)
}
DotPlotS <- function(project, useMatrix=NULL, assay=NULL, group.by=NULL, ...){
  if(is.null(useMatrix)){
    useMatrix = getAvailableMatrices(project)[1]
    warning(glue::glue("useMatrix is NULL, use first matrix {useMatrix}!"))
  }
  object <- PartialSeurat(project, useMatrix=useMatrix, assay=assay)
  p <- Seurat::DotPlot(object, assay="RNA", ...)
  return(p)

}
VlnPlotS <- function(project, useMatrix=NULL, assay=NULL, group.by=NULL, ...){
  if(is.null(useMatrix)){
    useMatrix = getAvailableMatrices(project)[1]
    warning(glue::glue("useMatrix is NULL, use first matrix {useMatrix}!"))
  }
  object <- PartialSeurat(project, useMatrix=useMatrix, assay=assay)
  p <- Seurat::VlnPlot(object, assay="RNA", slot=assay, ...)
  return(p)
}

