#' @return Dimension reduction
#' @param object Seurat or ArchR object
#' @rdname getDimRed
#' @export
getDimRed <- function(object, ...){
 UseMethod(generic = "getDimRed", object = object)
}

#' @return The object with dimension reduction
#' @param object Seurat or ArchR object
#' @rdname setDimRed
#' @export
setDimRed <- function(object, ...){
 UseMethod(generic = "setDimRed", object = object)
}

#' @return CellColDF
#' @param object Seurat or ArchR object
#' @rdname getCellCol
#' @export
getCellCol<- function(object, ...){
 UseMethod(generic = "getCellCol", object = object)
}

#' @return OBJ
#' @param object Seurat or ArchR object
#' @rdname setCellCol
#' @export
setCellCol <- function(object, ...){
 UseMethod(generic = "setCellCol", object = object)
}

#' @return Matrix
#' @param object Seurat or ArchR object
#' @rdname getMatrix
#' @export
getMatrix <- function(object, ...){
 UseMethod(generic = "getMatrix", object = object)
}

#' @return OBJ
#' @param object Seurat or ArchR object
#' @rdname setMatrix
#' @export
setMatrix  <- function(object, ...){
 UseMethod(generic = "setMatrix", object = object)
}

