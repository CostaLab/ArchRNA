

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
addBatchClusters <- function(project, resolutions=0.5, cluster_prefix, ...){

  for(resolution in resolutions){
    project <- addClusters(project,resolution=resolution, name=glue::glue("{cluster_prefix}{resolution}"), ...)
  }
  project
}
