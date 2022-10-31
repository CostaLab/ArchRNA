####################################################################
# Import 10X Data
####################################################################

#' Import Feature Matrix from 10x Feature HDF5 file.
#'
#' This function will import the feature matrix from a 10x feature hdf5 file.
#'
#' @param input A character of paths to 10x feature hdf5 file(s) or 10X sparse matrix directory. These will traditionally have a suffix similar to "filtered_feature_bc_matrix.h5".

#' @param input_type A character of input types 10XH5 or 10X
#' @param names A character of sample names associated with each input file.
#' @param strictMatch Only relevant when multiple input files are used. A boolean that indictes whether rows (genes) that do not match perfectly in the matrices
#' should be removed (`strictMatch = TRUE`) or coerced (`strictMatch = FALSE`). CellRanger seems to occassionally use different ensembl ids for the same gene across
#' different samples. If you are comfortable tolerating such mismatches, you can coerce all matrices to fit together, in which case the gene metadata present in
#' the first listed sample will be applied to all matrices for that particular gene entry. Regardless of what value is used for `strictMatch`, this function
#' cannot tolerate mismatched gene names, only mismatched metadata for the same gene.
#' @param verbose Only relevant when multiple input files are used. A boolean that indicates whether messaging about mismatches should be verbose (`TRUE`) or minimal (`FALSE`)
#' @param featureType The name of the feature to extract from the 10x feature file.
#' See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices for more information.
#' @export
import10xFeatureMatrix_ <- function(
  input = NULL,
  input_type = NULL, ## 10XH5 or 10X
  names = NULL,
  strictMatch = TRUE,
  verbose = TRUE,
  featureType = "Gene Expression"
){

  ArchR:::.validInput(input = input, name = "input", valid = c("character"))
  ArchR:::.validInput(input = names, name = "names", valid = c("character"))
  ArchR:::.validInput(input = strictMatch, name = "strictMatch", valid = c("boolean"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = featureType, name = "featureType", valid = c("character"))

  if (!all(file.exists(input))) {
    stop("Not all input file paths exist!")
  }
  if (input_type %ni% c("10XH5", "10X")){
    stop("not supported input type, options 10XH5, 10X!")
  }
  featureMats <- lapply(seq_along(input), function(y) {
    message("Importing Feature Matrix ", y, " of ", length(input))
    if(input_type == "10XH5"){
      .import10xh5FM_(featureMatrix = input[y], featureType = featureType,
              name = names[y])
    }else if(input_type == "10X"){
      .import10xFM_(featureMatrix = input[y], featureType = featureType,
              name = names[y])
    }

  })

  #message("Re-ordering RNA matricies for consistency.")
  #for(j in 1:length(featureMats)) {
  #  featureMats[[j]] <- sort.GenomicRanges(sortSeqlevels(featureMats[[j]]), ignore.strand = TRUE)
  #}

  #if more than one filtered feature barcode matrix is supplied, then merge the RSE objects
  if (length(featureMats) > 1) {
    message("Merging individual RNA objects...")
    #make the first matrix the base matrix and merge all others into it
    rse_final <- featureMats[[1]]

    rowsToRemove <- c() #rows that have previously been removed from rse_final

    #for each additional feature matrix (starting with the second), look for mismatches with rse_final and merge accordingly
    for (i in 2:length(featureMats)) {
      mismatchWarning <- TRUE #a boolean to prevent output of the warning message many times and only output it once

      message(sprintf("\nMerging %s", names[i]))

      if (!identical(rownames(rse_final), rownames(featureMats[[i]]))) {
        stop("Error - rownames (genes) of individual RNA objects are not equivalent.")
      }
      if (!identical(colnames(rowData(rse_final)), colnames(rowData(featureMats[[i]])))) {
        stop("Error - rowData (gene metadata) of individual RNA objects have different columns. This is highly unusual and merging has been aborted.")
      }
      if (!identical(names(assays(rse_final)), names(assays(featureMats[[i]])))) {
        stop("Error - available assays of individual RNA objects are not equivalent. Each object is expected to only have one assay named 'counts'.")
      }

      #check each column in rowData to check for mismatches that should be thrown as warnings
      #occasionally, it seems like 10x is annotating different ensembl IDs to the same gene which seems like a bad way to go
      #this is a bit heavy-handed but it seems like the safest thing to do is report any mismatch rather than merge blindly

      for (x in 1:ncol(rowData(rse_final))) {
        if (!identical(rowData(rse_final)[,x], rowData(featureMats[[i]])[,x])) {
          if(mismatchWarning) {
            message(sprintf("Warning! Some values within column \"%s\" of the rowData (gene metadata) of your objects do not precisely match!", colnames(rowData(rse_final))[x]))
            message("This is often caused by slight variations in Ensembl IDs and gene locations used by cellranger across different samples. ArchR will ignore these mismatches and allow merging to proceed but you should check to make sure that these are ok for your data.\n")
            mismatchWarning <- FALSE
          }

          #detect all of the mismatches betwenn rse_final and the current featureMat
          mismatch <- which(rowData(rse_final)[,x] != rowData(featureMats[[i]])[,x])
          #for each detected mismatch, handle the mismatch according to the value of strictMatch
          for (y in 1:length(mismatch)) {
            if (verbose) {
              message(sprintf("Mismatch in column \"%s\" row %s for %s: %s does not exactly match %s!", colnames(rowData(rse_final))[x], mismatch[y], names[i], rowData(rse_final)[mismatch[y],x], rowData(featureMats[[i]])[mismatch[y],x]))
            }
            if (strictMatch) {
              if (verbose) {
                message("strictMatch = TRUE so the corresponding gene entry with mismatching information will be removed.")
              }
              rowsToRemove <- unique(c(rowsToRemove, mismatch[y]))
              #temporarily force the data to match so that merging can occur easily. Mismatched rows will be removed later
              rowData(featureMats[[i]])[mismatch[y],] <- rowData(rse_final)[mismatch[y],]
              rowRanges(featureMats[[i]])[mismatch[y]] <- rowRanges(rse_final)[mismatch[y]]
            } else {
              if (verbose) {
                message("strictMatch = FALSE so mismatching information will be coerced to match the first sample provided.")
              }
              rowData(featureMats[[i]])[mismatch[y],] <- rowData(rse_final)[mismatch[y],]
              rowRanges(featureMats[[i]])[mismatch[y]] <- rowRanges(rse_final)[mismatch[y]]
            }
          }
        }
      }

      rse_final <- SummarizedExperiment::cbind(rse_final, featureMats[[i]])
    }
    if (strictMatch) {
      if(length(rowsToRemove) > 0) {
        rse_final <- rse_final[-rowsToRemove,]
      }
    }
    return(rse_final)
  }
  else {
    return(featureMats[[1]])
  }
}


.import10xh5FM_ <- function(featureMatrix = NULL, featureType = NULL, name = NULL){

  o <- h5closeAll()
  barcodes <- h5read(featureMatrix, "/matrix/barcodes")
  data <- h5read(featureMatrix, "/matrix/data")
  indices <- h5read(featureMatrix, "/matrix/indices")
  indptr <- h5read(featureMatrix, "/matrix/indptr")
  shape <- h5read(featureMatrix, "/matrix/shape")

  spMat <- sparseMatrix(
    i = indices,
    p = indptr,
    x = data,
    dims = shape,
    index1 = FALSE
  )

  colnames(spMat) <- paste0(name, "#", barcodes)

  features <- h5read(featureMatrix, "/matrix/features")
  features <- lapply(seq_along(features), function(x){
    if(length(features[[x]]) == nrow(spMat)){
      if(object.size(features[[x]]) > object.size(Rle(features[[x]]))){
        df <- DataFrame(x = Rle(features[[x]]))
      }else{
        df <- DataFrame(x = features[[x]])
      }
      colnames(df) <- names(features)[x]
      df
    }else{
      NULL
    }
  })
  features <- Reduce("cbind",features[!unlist(lapply(features,is.null))])

  se <- SummarizedExperiment(assays = SimpleList(counts = spMat), rowData = features)

  rownames(se) <- features$name

  if("feature_type" %in% colnames(rowData(se))){
    if(!is.null(featureType)){
      idx <- BiocGenerics::which(rowData(se)$feature_type %bcin% featureType)
      if(length(idx) == 0){
        stop("Error featureType not within provided features!")
      }
      se <- se[idx]
    }
  }

  if("interval" %in% colnames(rowData(se))){
    idxNA <- which(rowData(se)$interval=="NA")
    if(length(idxNA) > 0){
      se <- se[-idxNA, ]
    }
    rr <- GRanges(paste0(rowData(se)$interval))
    mcols(rr) <- rowData(se)
    se <- SummarizedExperiment(assays = SimpleList(counts = assay(se)), rowRanges = rr)
  }

  idxDup <- which(rownames(se) %in% rownames(se[duplicated(rownames(se))]))
  names(idxDup) <- rownames(se)[idxDup]
  if(length(idxDup) > 0){
    dupOrder <- idxDup[order(Matrix::rowSums(assay(se[idxDup])),decreasing=TRUE)]
    dupOrder <- dupOrder[!duplicated(names(dupOrder))]
    se <- se[-as.vector(idxDup[idxDup %ni% dupOrder])]
  }

  gc()

  se
}

.import10xFM_ <- function(featureMatrix = NULL, featureType = "Gene Expression", name = NULL){


  spMat <- Seurat::Read10X(featureMatrix)
  colnames(spMat) <- paste0(name, "#", colnames(spMat))
  features = DataFrame(feature_type="Gene Expression", name=rownames(spMat))
  se <- SummarizedExperiment(assays = SimpleList(counts = spMat), rowData = features)

  rownames(se) <- features$name

  #print(str(colnames(rowData(se))))
  if("feature_type" %in% colnames(rowData(se))){
    if(!is.null(featureType)){
      idx <- BiocGenerics::which(rowData(se)$feature_type %bcin% featureType)
      if(length(idx) == 0){
        stop("Error featureType not within provided features!")
      }
      se <- se[idx]
    }
  }

  idxDup <- which(rownames(se) %in% rownames(se[duplicated(rownames(se))]))
  names(idxDup) <- rownames(se)[idxDup]
  if(length(idxDup) > 0){
    dupOrder <- idxDup[order(Matrix::rowSums(assay(se[idxDup])),decreasing=TRUE)]
    dupOrder <- dupOrder[!duplicated(names(dupOrder))]
    se <- se[-as.vector(idxDup[idxDup %ni% dupOrder])]
  }

  gc()
  se
}


.importSparseFM_ <- function(spMat = NULL, name = NULL, featureType = "Gene Expression", addHashtag=TRUE){
  if(addHashtag){
    colnames(spMat) <- paste0(name, "#", colnames(spMat))
  }
  features = DataFrame(feature_type="Gene Expression", name=rownames(spMat))
  se <- SummarizedExperiment(assays = SimpleList(counts = spMat), rowData = features)

  rownames(se) <- features$name

  if("feature_type" %in% colnames(rowData(se))){
    if(!is.null(featureType)){
      idx <- BiocGenerics::which(rowData(se)$feature_type %bcin% featureType)
      if(length(idx) == 0){
        stop("Error featureType not within provided features!")
      }
      se <- se[idx]
    }
  }

  idxDup <- which(rownames(se) %in% rownames(se[duplicated(rownames(se))]))
  names(idxDup) <- rownames(se)[idxDup]
  if(length(idxDup) > 0){
    dupOrder <- idxDup[order(Matrix::rowSums(assay(se[idxDup])),decreasing=TRUE)]
    dupOrder <- dupOrder[!duplicated(names(dupOrder))]
    se <- se[-as.vector(idxDup[idxDup %ni% dupOrder])]
  }

  gc()
  se
}

