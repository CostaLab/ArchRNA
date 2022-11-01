

#' Single-Cell Utilities Constructor
#'
#' @import methods
#' @import ArchR
#' @import data.table
#'
#' @param project ArchR project with meta-data
#' @param resolutions vector
#'
#' @rdname addBatchClusters
#' @export
addBatchClusters <- function(project, resolutions=0.5, ...){

  for(resolution in resolutions){
    project <- addClusters(project,resolution=resolution, name=glue::glue("{cluster_prefix}_{resolution}"), ...)
  }
  project
}
