
getMatrixFeatures <- function(
  ArchRProj = NULL,
  matrixName = NULL,
  seqname = NULL, ## if null
  features = NULL,
  threads = getArchRThreads(),
  logFile = NULL
  ){

  o <- h5closeAll()

  ArchR:::.logMessage("Getting Matrix Values...", verbose = TRUE, logFile = logFile)

  featureDF <- ArchR:::.getFeatureDF(head(getArrowFiles(ArchRProj), 2), matrixName)
  ArchR:::.logThis(featureDF, "FeatureDF", logFile = logFile)

  matrixClass <- h5read(getArrowFiles(ArchRProj)[1], paste0(matrixName, "/Info/Class"))

  if(matrixClass == "Sparse.Assays.Matrix"){
      if(is.null(seqname)){
      ArchR:::.logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available formats for input!", logFile = logFile)
      stop("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available formats for input!")
    }
  }

  if(!is.null(seqname)){
    idx <- lapply(seq_along(features), function(x){
      ix <- intersect(which(tolower(features[x]) == tolower(featureDF$name)), BiocGenerics::which(tolower(seqname) == tolower(featureDF$seqnames)))
      if(length(ix)==0){
        ArchR:::.logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", features[x]), logFile = logFile)
      }
      ix
    }) %>% unlist

  }else{

    idx <- lapply(seq_along(features), function(x){
      ix <- which(tolower(features[x]) == tolower(featureDF$name))[1]
      if(length(ix)==0){
        ArchR:::.logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", features[x]), logFile = logFile)
      }
      ix
    }) %>% unlist

  }
  ArchR:::.logThis(idx, "idx", logFile = logFile)

  if(any(is.na(idx))){
    ArchR:::.logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", paste0(features[which(is.na(idx))], collapse=",")), logFile = logFile)
  }

  featureDF <- featureDF[idx, ,drop=FALSE]
  ArchR:::.logThis(featureDF, "FeatureDF-Subset", logFile = logFile)

  #Get Values for FeatureName
  cellNamesList <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj)$Sample)

  values <- ArchR:::.safelapply(seq_along(cellNamesList), function(x){
    #message(x, " ", appendLF = FALSE)
    valuesx <- tryCatch({
      o <- h5closeAll()
      ArrowFile <- getSampleColData(ArchRProj)[names(cellNamesList)[x],"ArrowFiles"]
      valuesx <- ArchR:::.getMatFromArrow(
          ArrowFile = ArrowFile,
          featureDF = featureDF,
          binarize = FALSE,
          useMatrix = matrixName,
          cellNames = cellNamesList[[x]],
          threads = 1
        )
      colnames(valuesx) <- cellNamesList[[x]]
      valuesx
    }, error = function(e){
      errorList <- list(
        x = x,
        ArrowFile = ArrowFile,
        ArchRProj = ArchRProj,
        cellNames = ArchRProj$cellNames,
        cellNamesList = cellNamesList,
        featureDF = featureDF
      )
      ArchR:::.logError(e, fn = ".getMatFromArrow", info = "", errorList = errorList, logFile = logFile)
    })
    valuesx
  }, threads = threads) %>% Reduce("cbind", .)
  values <- values[, ArchRProj$cellNames, drop = FALSE]
  message("")
  gc()
  ArchR:::.logThis(values, "Feature-Matrix", logFile = logFile)

  if(!inherits(values, "matrix")){
    values <- matrix(as.matrix(values), ncol = nCells(ArchRProj))
    colnames(values) <- ArchRProj$cellNames
  }

  #Values Summary

  rownames(values) <- features

  return(values)

}


