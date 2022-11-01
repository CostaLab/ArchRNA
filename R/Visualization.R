#' @export
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


#' @export
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


#' @export
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

#' @export
RidgePlotS <-function(project, useMatrix=NULL, assay=NULL, group.by=NULL, ...){
  if(is.null(useMatrix)){
    useMatrix = getAvailableMatrices(project)[1]
    warning(glue::glue("useMatrix is NULL, use first matrix {useMatrix}!"))
  }
  object <- PartialSeurat(project, useMatrix=useMatrix, assay=assay)
  p <- Seurat::RidgePlot(object,assay="RNA", slot=assay, ...)
  return(p)
}

#' @export
DotPlotS <- function(project, useMatrix=NULL, assay=NULL, group.by=NULL, ...){
  if(is.null(useMatrix)){
    useMatrix = getAvailableMatrices(project)[1]
    warning(glue::glue("useMatrix is NULL, use first matrix {useMatrix}!"))
  }
  object <- PartialSeurat(project, useMatrix=useMatrix, assay=assay)
  p <- Seurat::DotPlot(object, assay="RNA", ...)
  return(p)

}

#' @export
VlnPlotS <- function(project, useMatrix=NULL, assay=NULL, group.by=NULL, ...){

  ### if vln in meta.data, there's no need to retrieve count matrix
  if(is.null(useMatrix)){
    useMatrix = getAvailableMatrices(project)[1]
    warning(glue::glue("useMatrix is NULL, use first matrix {useMatrix}!"))
  }
  object <- PartialSeurat(project, useMatrix=useMatrix, assay=assay)
  p <- Seurat::VlnPlot(object, assay="RNA", slot=assay, ...)
  return(p)
}


#' @export
scProportionPlotS <- function(project,
                              clusterName="Cluster_0.5",
                              condition="Sample",
                              pair=c("a", "b"),
                              logFile = ArchR::createLogFile("scProportionPlotS")
){

  prop_test <- sc_utils(project)
  #prop_test@meta_data$sample_identity = prop_test@meta_data[, condition]
  #prop_test@meta_data$cluster_identity = prop_test@meta_data[, clusterName]

  prop_test@meta_data <-  prop_test@meta_data[eval(as.name(condition)) %in% pair ]

  ArchR:::.logMessage("run permutation test  ", logfile=logFile, verbose=T)
  prop_test <- permutation_test(prop_test,
                                cluster_identity = clusterName,
                                sample_1 = pair[1],
                                sample_2 = pair[2],
                                sample_identity = condition)

  ArchR:::.logMessage("finished permutation test  ", logfile=logFile, verbose=T)
  dff = data.table::as.data.table(prop_test@results$permutation)
  p <- permutation_plot(dff) + ggtitle(glue::glue("{pair[1]} vs {pair[2]}, positive means more {pair[2]}"))
  p
}
#Volcano

VolcanoPlotS <- function(){

}

##genebar
FlipBarPlotS <- function(){

}

GroupedBarPlotS <- function(){
  #each group is
  #proportions in each sample or condition in each cluster
  #
}

PropBarPlotS <- function(){

}

##
DensityPlotS<- function(project){
library(Nebulosa)
}

##
VlnMatrixPlotS <- function(){

}

##
PiePlotS <- function(){

}




## each Dot is a piePlots: Li's paper
PieMatrixPlotS <- function(){
}

## high variable genes



##clusteree


#' @export
clustreeS <- function(project, prefix="", suffix="",...){
  require(clustree, quietly=T)
  coldata = as.data.frame(project@cellColData)
  clustree(coldata, prefix=prefix, suffix=suffix, ...)
}



## statistics
#name 	nCount.Mean 	nCount.Median 	nFeature.Mean 	nFeature.Median 	pctMt.Mean 	pctMt.Median 	pctRb.Mean 	pctRb.Median 	Cells
#A_MxCre 	5717.538 	5698.0 	1985.620 	1985 	3.268716 	3.213078 	38.79467 	38.99316 	999
#B_MxCre 	5689.311 	5629.0 	1986.131 	1986 	2.717589 	2.679856 	38.03635 	38.15975 	793
#C_Csnk 	5853.518 	5810.5 	1981.780 	1981 	3.306313 	3.329845 	39.06403 	39.13145 	628
#D_Csnk 	5737.073 	5688.0 	1981.088 	1979 	3.377221 	3.351727 	38.33313 	38.67722 	478

