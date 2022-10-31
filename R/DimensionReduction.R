addPCA <- function(project, useMatrix, assay){
  #    total.variance <- sum(RowVar(x = object))
  #  if (approx) {
  #    npcs <- min(npcs, nrow(x = object) - 1)
  #    pca.results <- irlba(A = t(x = object), nv = npcs, ...)
  #    feature.loadings <- pca.results$v
  #    sdev <- pca.results$d/sqrt(max(1, ncol(object) - 1))
  #    if (weight.by.var) {
  #      cell.embeddings <- pca.results$u %*% diag(pca.results$d)
  #    } else {
  #      cell.embeddings <- pca.results$u
  #    }
  #  } else {
  #    npcs <- min(npcs, nrow(x = object))
  #    pca.results <- prcomp(x = t(object), rank. = npcs, ...)
  #    feature.loadings <- pca.results$rotation
  #    sdev <- pca.results$sdev
  #    if (weight.by.var) {
  #      cell.embeddings <- pca.results$x
  #    } else {
  #      cell.embeddings <- pca.results$x / (pca.results$sdev[1:npcs] * sqrt(x = ncol(x = object) - 1))
  #    }
  #  }
  #}
  #rownames(x = feature.loadings) <- rownames(x = object)
  #colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
  #rownames(x = cell.embeddings) <- colnames(x = object)
  #colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  #reduction.data <- CreateDimReducObject(
  #  embeddings = cell.embeddings,
  #  loadings = feature.loadings,
  #  assay = assay,
  #  stdev = sdev,
  #  key = reduction.key,
  #  misc = list(total.variance = total.variance)
  #)
  gc()
   return(0)
}

addICA <- function(project, useMatrix, assay){
  return(0)
}

addDiffusion <- function(project, useMatrix, assay){
  return(0)
}




##########################################################################################
# LSI Dimensionality Reduction Methods
##########################################################################################
## adjust >=2 assays in one matrix

#' Add an Iterative LSI-based dimensionality reduction to an ArchRProject
#'
#' This function will compute an iterative LSI dimensionality reduction on an ArchRProject.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix The name of the data matrix to retrieve from the ArrowFiles associated with the `ArchRProject`. Valid options are
#' "TileMatrix" or "PeakMatrix".
#' @param name The name to use for storage of the IterativeLSI dimensionality reduction in the `ArchRProject` as a `reducedDims` object.
#' @param iterations The number of LSI iterations to perform.
#' @param clusterParams A list of Additional parameters to be passed to `addClusters()` for clustering within each iteration.
#' These params can be constant across each iteration, or specified for each iteration individually. Thus each param must be of
#' length == 1 or the total number of `iterations` - 1. PLEASE NOTE - We have updated these params to `resolution=2` and `maxClusters=6`! To use previous settings use `resolution=0.2` and `maxClusters=NULL`.
#' @param firstSelection First iteration selection method for features to use for LSI. Either "Top" for the top accessible/average or "Var" for the top variable features.
#' "Top" should be used for all scATAC-seq data (binary) while "Var" should be used for all scRNA/other-seq data types (non-binary).
#' @param depthCol A column in the `ArchRProject` that represents the coverage (scATAC = unique fragments, scRNA = unique molecular identifiers) per cell.
#' These values are used to minimize the related biases in the reduction related. For scATAC we recommend "nFrags" and for scRNA we recommend "Gex_nUMI".
#' @param varFeatures The number of N variable features to use for LSI. The top N features will be used based on the `selectionMethod`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param LSIMethod A number or string indicating the order of operations in the TF-IDF normalization.
#' Possible values are: 1 or "tf-logidf", 2 or "log(tf-idf)", and 3 or "logtf-logidf".
#' @param scaleDims A boolean that indicates whether to z-score the reduced dimensions for each cell. This is useful forminimizing the contribution
#' of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since
#' it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the `reducedDims` were
#' originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param binarize A boolean value indicating whether the matrix should be binarized before running LSI. This is often desired when working with insertion counts.
#' @param outlierQuantiles Two numerical values (between 0 and 1) that describe the lower and upper quantiles of bias (number of acessible regions per cell, determined
#' by `nFrags` or `colSums`) to filter cells prior to LSI. For example a value of c(0.02, 0.98) results in the cells in the bottom 2 percent and upper 98 percent to be
#' filtered prior to LSI. These cells are then projected back in the LSI subspace. This prevents spurious 'islands' that are identified based on being extremely biased.
#' These quantiles are also used for sub-sampled LSI when determining which cells are used.
#' @param filterBias A boolean indicating whether to drop bias clusters when computing clusters during iterativeLSI.
#' @param sampleCellsPre An integer specifying the number of cells to sample in iterations prior to the last in order to perform a sub-sampled LSI and
#' sub-sampled clustering. This greatly reduced memory usage and increases speed for early iterations.
#' @param projectCellsPre A boolean indicating whether to reproject all cells into the sub-sampled LSI (see `sampleCellsPre`). Setting this to `FALSE`
#' allows for using the sub-sampled LSI for clustering and variance identification. If `TRUE` the cells are all projected into the sub-sampled LSI
#' and used for cluster and variance identification.
#' @param sampleCellsFinal An integer specifying the number of cells to sample in order to perform a sub-sampled LSI in final iteration.
#' @param selectionMethod The selection method to be used for identifying the top variable features. Valid options are "var" for
#' log-variability or "vmr" for variance-to-mean ratio.
#' @param scaleTo Each column in the matrix designated by `useMatrix` will be normalized to a column sum designated by `scaleTo` prior to
#' variance calculation and TF-IDF normalization.
#' @param totalFeatures The number of features to consider for use in LSI after ranking the features by the total number of insertions.
#' These features are the only ones used throught the variance identification and LSI. These are an equivalent when using a `TileMatrix` to a defined peakSet.
#' @param filterQuantile A number [0,1] that indicates the quantile above which features should be removed based on insertion counts prior
#' @param excludeChr A string of chromosomes to exclude for iterativeLSI procedure.
#' to the first iteration of the iterative LSI paradigm. For example, if `filterQuantile = 0.99`, any features above the 99th percentile in
#' insertion counts will be ignored for the first LSI iteration.
#' @param saveIterations A boolean value indicating whether the results of each LSI iterations should be saved as compressed `.rds` files in
#' the designated `outDir`.
#' @param UMAPParams The list of parameters to pass to the UMAP function if "UMAP" if `saveIterations=TRUE`. See the function `uwot::umap()`.
#' @param nPlot If `saveIterations=TRUE`, how many cells to sample make a UMAP and plot for each iteration.
#' @param outDir The output directory for saving LSI iterations if desired. Default is in the `outputDirectory` of the `ArchRProject`.
#' @param threads The number of threads to be used for parallel computing.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can
#' reproduce results downstream.
#' @param verbose A boolean value that determines whether standard output includes verbose sections.
#' @param force A boolean value that indicates whether or not to overwrite relevant data in the `ArchRProject` object.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addRNAIterativeLSI <- function(
  ArchRProj = NULL,
  useMatrix = "GeneExpressionMatrix",
  assay = "data",
  name = "RNAIterativeLSI",
  iterations = 2,
  clusterParams = list(
      resolution = 0.2,
      sampleCells = 10000,
      maxClusters = 6,
      n.start = 10
  ),
  firstSelection = "variable",
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  dimsToUse = 1:30,
  LSIMethod = 2,
  scaleDims = TRUE,
  corCutOff = 0.75,
  binarize = FALSE,
  outlierQuantiles = c(0.02, 0.98),
  filterBias = TRUE,
  sampleCellsPre = 10000,
  projectCellsPre = FALSE,
  sampleCellsFinal = NULL,
  selectionMethod = "var",
  scaleTo = 10000,
  totalFeatures = 500000,
  filterQuantile = 0.995,
  excludeChr = c(),
  saveIterations = F,
  UMAPParams = list(
    n_neighbors = 40,
    min_dist = 0.4,
    metric = "cosine",
    verbose = FALSE,
    fast_sgd = TRUE
  ),
  nPlot = 10000,
  outDir = getOutputDirectory(ArchRProj),
  threads = getArchRThreads(),
  seed = 1,
  verbose = TRUE,
  force = FALSE,
  logFile = createLogFile("addIterativeLSI")
  ){

  if(verbose) message("Checking Inputs...")
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  ArchR:::.validInput(input = name, name = "name", valid = c("character"))
  ArchR:::.validInput(input = iterations, name = "iterations", valid = c("integer"))
  ArchR:::.validInput(input = clusterParams, name = "clusterParams", valid = c("list"))
  ArchR:::.validInput(input = varFeatures, name = "varFeatures", valid = c("integer"))
  ArchR:::.validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer"))
  ArchR:::.validInput(input = LSIMethod, name = "LSIMethod", valid = c("integer", "character"))
  ArchR:::.validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  ArchR:::.validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
  ArchR:::.validInput(input = binarize, name = "binarize", valid = c("boolean"))
  ArchR:::.validInput(input = outlierQuantiles, name = "outlierQuantiles", valid = c("numeric", "null"))
  ArchR:::.validInput(input = filterBias, name = "filterBias", valid = c("boolean"))
  ArchR:::.validInput(input = sampleCellsPre, name = "sampleCellsPre", valid = c("integer", "null"))
  ArchR:::.validInput(input = sampleCellsFinal, name = "sampleCellsFinal", valid = c("integer", "null"))
  ArchR:::.validInput(input = selectionMethod, name = "selectionMethod", valid = c("character"))
  ArchR:::.validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  ArchR:::.validInput(input = totalFeatures, name = "totalFeatures", valid = c("integer"))
  ArchR:::.validInput(input = filterQuantile, name = "filterQuantile", valid = c("numeric"))
  ArchR:::.validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  ArchR:::.validInput(input = saveIterations, name = "saveIterations", valid = c("boolean"))
  ArchR:::.validInput(input = UMAPParams, name = "UMAPParams", valid = c("list"))
  ArchR:::.validInput(input = nPlot, name = "nPlot", valid = c("integer"))
  ArchR:::.validInput(input = outDir, name = "outDir", valid = c("character"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
  ArchR:::.validInput(input = seed, name = "seed", valid = c("integer"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))


  if(varFeatures < 1000){
    stop("Please provide more than 1000 varFeatures!")
  }

  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "IterativeLSI Input-Parameters", logFile=logFile)

  ArchR:::.requirePackage("Matrix", source = "cran")
  tstart <- Sys.time()

  if(!is.null(ArchRProj@reducedDims[[name]])){
    if(!force){
      stop("Error name in reducedDims Already Exists! Set force = TRUE or pick a different name!")
    }
  }

  #Set Seed
  set.seed(seed)
  outDir <- file.path(outDir, name)
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

  #All the Cell Names
  cellNames <- rownames(getCellColData(ArchRProj))
  if(!is.null(sampleCellsPre)){
    if(length(cellNames) < sampleCellsPre){
      sampleCellsPre <- NULL
    }
  }
  if(!is.null(sampleCellsFinal)){
    if(length(cellNames) < sampleCellsFinal){
      sampleCellsFinal <- NULL
    }
  }

  #MatrixFiles
  ArrowFiles <- getSampleColData(ArchRProj)[,"ArrowFiles"]

  #Check if Matrix is supported and check type
  if(tolower(useMatrix) == "tilematrix"){
    useMatrix <- "TileMatrix"
    tileSizes <- lapply(ArrowFiles, function(x){
      h5read(x, "TileMatrix/Info/Params/")$tileSize[1]
    }) %>% unlist
    if(length(unique(tileSizes)) != 1){
      stop("Error not all TileMatrices are the same tileSize!")
    }
    tileSize <- unique(tileSizes)
  }else if(tolower(useMatrix) == "peakmatrix"){
    useMatrix <- "PeakMatrix"
    tileSize <- NA
  }else{
    tileSize <- NA
  }

  units <- unique(unlist(lapply(ArrowFiles, function(x) h5read(x, paste0(useMatrix, "/Info/Units")))))
  if(assay %ni% units){
    stop(glue::glue("There's no assay of this matrix! options{units}"))
  }
  if(grepl("log",units,ignore.case=TRUE)){
    stop("Cannot use log transformed values for iterativeLSI!")
  }

  tstart <- Sys.time()
  ############################################################################################################################
  # Organize Information for LSI
  ############################################################################################################################
  chrToRun <- ArchR:::.availableSeqnames(ArrowFiles, subGroup = useMatrix)
  if(assay %in% chrToRun){
    chrToRun = assay
  }else{
    stop(glue::glue("There's no assay of this matrix! options{chrToRun}"))
  }

  if(tolower(firstSelection) == "top"){

    if(!binarize){
      stop("Please binarize data if using top selection for first iteration! Set binarize = TRUE!")
    }

    #Compute Row Sums Across All Samples
    ArchR:::.logDiffTime("Computing Total Across All Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    if(useMatrix == "TileMatrix"){
      totalAcc <- ArchR:::.getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = FALSE)
      totalAcc$start <- (totalAcc$idx - 1) * tileSize
    }else{
      totalAcc <- ArchR:::.getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = TRUE)
    }

    #Filter Chromosomes
    if(length(excludeChr) > 0){
      totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% excludeChr), , drop = FALSE]
    }

    #Identify the top features to be used here
    ArchR:::.logDiffTime("Computing Top Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    nFeature <- varFeatures[1]
    rmTop <- floor((1-filterQuantile) * totalFeatures)
    topIdx <- head(order(totalAcc$rowSums, decreasing=TRUE), nFeature + rmTop)[-seq_len(rmTop)]
    topFeatures <- totalAcc[sort(topIdx),]

    gc()

  }else if(tolower(firstSelection) %in% c("var", "variable")){

    if(binarize){
      stop("Please do not binarize data if using variable selection for first iteration! Set binarize = FALSE!")
    }

    if(units %in% "BinarizedCounts"){
      stop("Cannot do variable selection with BinarizedCounts. Set firstSelection = Top!")
    }

    #Compute Row Sums Across All Samples
    ArchR:::.logDiffTime("Computing Variability Across All Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    if(useMatrix == "TileMatrix"){
      totalAcc <- ArchR:::.getRowVars(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, useLog2 = TRUE)
      totalAcc$start <- (totalAcc$idx - 1) * tileSize
    }else{
      totalAcc <- ArchR:::.getRowVars(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, useLog2 = TRUE)
    }

    #Filter Chromosomes
    if(length(excludeChr) > 0){
      totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% excludeChr), , drop = FALSE]
    }

    #Identify the top features to be used here
    ArchR:::.logDiffTime("Computing Variable Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    nFeature <- varFeatures[1]
    if(nFeature > 0.5 * nrow(totalAcc)){
      stop("nFeature for variable selection must be at leat 1/2 the total features!")
    }
    topIdx <- head(order(totalAcc$combinedVars, decreasing=TRUE), nFeature)
    topFeatures <- totalAcc[sort(topIdx),]

    gc()

  }else{

    stop("firstSelect method must be Top or Var/Variable!")

  }

  cellDepth <- tryCatch({
      df <- getCellColData(ArchRProj = ArchRProj, select = depthCol)
      v <- df[,1]
      names(v) <- rownames(df)
      v
    }, error = function(e){
      tryCatch({
        ArchR:::.getColSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun)
      }, error = function(y){
        stop("Could not determine depth from depthCol or colSums!")
      })
    }
  )
  cellDepth <- log10(cellDepth + 1)

  ############################################################################################################################
  # LSI Iteration 1
  ############################################################################################################################
  ArchR:::.logDiffTime(paste0("Running LSI (1 of ",iterations,") on Top Features"), tstart, addHeader = TRUE, verbose = verbose, logFile = logFile)
  j <- 1

  if(!is.null(clusterParams$sampleCells)){
    if(!is.na(clusterParams$sampleCells[j])){
      sampleJ <- clusterParams$sampleCells[j]
    }else if(!is.na(clusterParams$sampleCells[1])){
      sampleJ <- clusterParams$sampleCells[1]
    }else{
      sampleJ <- sampleCellsPre
    }
  }else{
    sampleJ <- sampleCellsPre
  }

  outLSI <- ArchR:::.LSIPartialMatrix(
    ArrowFiles = ArrowFiles,
    featureDF = topFeatures,
    cellNames = cellNames,
    cellDepth = cellDepth,
    useMatrix = useMatrix,
    sampleNames = getCellColData(ArchRProj)$Sample,
    LSIMethod = LSIMethod,
    scaleTo = scaleTo,
    dimsToUse = dimsToUse,
    binarize = binarize,
    outlierQuantiles = outlierQuantiles,
    sampleCells = if(j != iterations) sampleCellsPre else sampleCellsFinal,
    projectAll = j == iterations | projectCellsPre | sampleJ > sampleCellsPre,
    threads = threads,
    useIndex = FALSE,
    seed = seed,
    tstart = tstart,
    verbose = verbose,
    logFile = logFile
  )
  outLSI$scaleDims <- scaleDims
  outLSI$useMatrix <- useMatrix
  outLSI$tileSize <- tileSize
  gc()
  ArchR:::.logThis(outLSI, paste0("outLSI-",j), logFile = logFile)

  if(iterations == 1){
    ArchR:::.logDiffTime("Finished Running addRNAIterativeLSI", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    ArchRProj@reducedDims[[name]] <- outLSI
    #.endLogging(logFile = logFile)
    return(ArchRProj)
  }

  #########################
  # Identify LSI Clusters
  #########################
  clusterDF <- ArchR:::.LSICluster(
    outLSI = outLSI,
    filterBias = filterBias,
    cellNames = cellNames,
    cellDepth = cellDepth,
    dimsToUse = dimsToUse,
    scaleDims = scaleDims,
    corCutOff = corCutOff,
    clusterParams = clusterParams,
    j = j,
    verbose = verbose,
    tstart = tstart,
    logFile = logFile
  )
  clusters <- clusterDF$clusters
  nClust <- length(unique(clusters))
  ArchR:::.logDiffTime(sprintf("Identified %s Clusters", nClust), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  ArchR:::.logThis(clusterDF, paste0("clusterDF-",j), logFile = logFile)

  #########################
  # Save LSI Iteration
  #########################
  if(saveIterations){
    ArchR:::.logDiffTime("Saving LSI Iteration", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    ArchR:::.saveIteration(outLSI=outLSI, clusters=clusters, scaleDims = scaleDims,
      dimsToUse = dimsToUse, corCutOff = corCutOff, outDir = outDir,
      nPlot=nPlot, UMAPParams=UMAPParams, ArchRProj=ArchRProj, j = j, threads = threads, logFile = logFile)
  }

  ############################################################################################################################
  # LSI Iteration 2+
  ############################################################################################################################
  variableFeatures <- topFeatures

  while(j < iterations){

    #Jth iteration
    j <- j + 1

    #########################
    # Identify Features for LSI Iteration
    #########################
    variableFeatures <- ArchR:::.identifyVarFeatures(
      outLSI = outLSI,
      clusterDF = clusterDF,
      ArrowFiles = ArrowFiles,
      useMatrix = useMatrix,
      prevFeatures = variableFeatures,
      scaleTo = scaleTo,
      totalAcc = totalAcc,
      totalFeatures = totalFeatures,
      firstSelection = firstSelection,
      selectionMethod = selectionMethod,
      varFeatures = varFeatures,
      tstart = tstart,
      threads = threads,
      verbose = verbose,
      logFile = logFile
    )

    #########################
    # LSI
    #########################
    ArchR:::.logDiffTime(sprintf("Running LSI (%s of %s) on Variable Features", j, iterations), tstart, addHeader = TRUE, verbose = verbose, logFile = logFile)
    if(!is.null(clusterParams$sampleCells)){
      if(!is.na(clusterParams$sampleCells[j])){
        sampleJ <- clusterParams$sampleCells[j]
      }else if(!is.na(clusterParams$sampleCells[1])){
        sampleJ <- clusterParams$sampleCells[1]
      }else{
        sampleJ <- sampleCellsPre
      }
    }else{
      sampleJ <- sampleCellsPre
    }

    #Compute Partial Matrix LSI
    outLSI <- ArchR:::.LSIPartialMatrix(
      ArrowFiles = ArrowFiles,
      featureDF = variableFeatures,
      useMatrix = useMatrix,
      cellNames = cellNames,
      cellDepth = cellDepth,
      sampleNames = getCellColData(ArchRProj)$Sample,
      LSIMethod = LSIMethod,
      scaleTo = scaleTo,
      dimsToUse = dimsToUse,
      binarize = binarize,
      outlierQuantiles = outlierQuantiles,
      sampleCells = if(j != iterations) sampleCellsPre else sampleCellsFinal,
      projectAll = j == iterations | projectCellsPre | sampleJ > sampleCellsPre,
      threads = threads,
      useIndex = FALSE,
      seed = seed,
      tstart = tstart,
      verbose = verbose,
      logFile = logFile
    )
    outLSI$scaleDims <- scaleDims
    outLSI$useMatrix <- useMatrix
    outLSI$tileSize <- tileSize
    ArchR:::.logThis(outLSI, paste0("outLSI-",j), logFile = logFile)

    if(j != iterations){

      #########################
      # Identify LSI Clusters
      #########################
      clusterDF <- ArchR:::.LSICluster(
        outLSI = outLSI,
        dimsToUse = dimsToUse,
        scaleDims = scaleDims,
        corCutOff = corCutOff,
        filterBias = filterBias,
        cellNames = cellNames,
        cellDepth = cellDepth,
        j = j,
        clusterParams = clusterParams,
        verbose = verbose,
        tstart = tstart,
        logFile = logFile
      )
      clusters <- clusterDF$clusters
      nClust <- length(unique(clusters))
      ArchR:::.logDiffTime(sprintf("Identified %s Clusters", nClust), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      ArchR:::.logThis(clusterDF, paste0("clusterDF-",j), logFile = logFile)

      #########################
      # Save LSI Iteration
      #########################
      if(saveIterations){
        ArchR:::.logDiffTime("Saving LSI Iteration", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        ArchR:::.saveIteration(outLSI=outLSI, clusters=clusters, scaleDims = scaleDims,
          dimsToUse = dimsToUse, corCutOff = corCutOff, outDir = outDir,
          nPlot=nPlot, UMAPParams=UMAPParams, ArchRProj=ArchRProj, j = j, threads = threads, logFile = logFile)
      }

    }

    gc()

  }

  #Organize Output
  ArchR:::.logDiffTime("Finished Running addRNAIterativeLSI", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  ArchRProj@reducedDims[[name]] <- outLSI

  return(ArchRProj)

}
