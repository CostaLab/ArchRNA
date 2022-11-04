####################################################################
# Trajectory Analysis Methods
####################################################################


#' Visualize a Trajectory from ArchR Project
#'
#' This function will plot a trajectory that was created onto an embedding.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param embedding The name of the embedding to use to visualize the given `trajectory`. See `addEmbedding()` for more information.
#' @param trajectory The column name in `cellColData` that refers the trajectory to be plotted. See `addTrajectory()` for more information.
#' @param colorBy A string indicating whether points in the plot should be colored by a column in `cellColData` ("cellColData")
#' or by a data matrix in the associated ArrowFiles (i.e. "GeneScoreMatrix", "MotifMatrix", "PeakMatrix").
#' @param name The name of the column in `cellColData` or the featureName/rowname of the data matrix to be used for plotting.
#' For example if colorBy is "cellColData" then `name` refers to a column name in the cellcoldata (see `getCellcoldata()`). If `colorBy`
#' is "GeneScoreMatrix" then `name` refers to a gene name which can be listed by `getFeatures(ArchRProj, useMatrix = "GeneScoreMatrix")`.
#' @param log2Norm A boolean value indicating whether a log2 transformation should be performed on the values from `colorBy`.
#' @param imputeWeights The weights to be used for imputing numerical values for each cell as a linear combination of other cells'
#' values. See `addImputationWeights()` and `getImutationWeights()` for more information.
#' @param pal The name of a custom palette from `ArchRPalettes` to use for coloring cells.
#' @param size A number indicating the size of the points to plot if `plotAs` is set to "points".
#' @param rastr A boolean value that indicates whether the plot should be rasterized. This does not rasterize lines and labels,
#' just the internal portions of the plot.
#' @param quantCut If this is not `NULL`, a quantile cut is performed to threshold the top and bottom of the distribution of numerical values.
#' This prevents skewed color scales caused by strong outliers. The format of this should be c(x,y) where x is the lower threshold and y is
#' the upper threshold. For example, quantileCut = c(0.025,0.975) will take the 2.5th percentile and 97.5 percentile of values and
#' set values below/above to the value of the 2.5th and 97.5th percentile values respectively.
#' @param quantHex The numeric xth quantile of all dots within each individual hexagon will determine the numerical value for
#' coloring to be displayed. This occurs when (i) `plotAs` is set to "hex" or (ii) `plotAs` is set to `NULL` and the values of `colorBy` are numeric.
#' @param discreteSet The name of a discrete palette from `ArchRPalettes` for visualizing `colorBy` in the embedding if a discrete color set is desired.
#' @param continuousSet The name of a continuous palette from `ArchRPalettes` for visualizing `colorBy` in the embedding if a continuous color set is desired.
#' @param randomize A boolean value that indicates whether to randomize points prior to plotting to prevent cells from one cluster
#' being present at the front of the plot.
#' @param keepAxis A boolean value that indicates whether the x and y axis ticks and labels should be plotted.
#' @param baseSize The base font size to use in the plot.
#' @param addArrow A boolean value that indicates whether to add a smoothed arrow in the embedding based on the aligned trajectory.
#' @param plotAs A string that indicates whether points ("points") should be plotted or a hexplot ("hex") should be plotted. By default
#' if `colorBy` is numeric, then `plotAs` is set to "hex".
#' @param smoothWindow An integer value indicating the smoothing window for creating inferred Arrow overlay on to embedding.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @param ... Additional parameters to pass to `ggPoint()` or `ggHex()`.
#' @export
plotTrajectory_ <- function(
  ArchRProj = NULL,
  embedding = "UMAP",
  trajectory = "Trajectory",
  colorBy = "colData",
  colorMatrixAssay="data", ## only works when colored by Matrix
  name = "Trajectory",
  log2Norm = NULL,
  imputeWeights = if(!grepl("coldata",tolower(colorBy[1]))) getImputeWeights(ArchRProj),
  pal = NULL,
  size = 0.2,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  quantHex = 0.5,
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 6,
  addArrow = TRUE,
  plotAs = NULL,
  smoothWindow = 5,
  logFile = createLogFile("plotTrajectory"),
  ...
  ){

  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = embedding, name = "reducedDims", valid = c("character"))
  ArchR:::.validInput(input = trajectory, name = "trajectory", valid = c("character"))
  ArchR:::.validInput(input = colorBy, name = "colorBy", valid = c("character"))
  ArchR:::.validInput(input = name, name = "name", valid = c("character"))
  ArchR:::.validInput(input = log2Norm, name = "log2Norm", valid = c("boolean", "null"))
  ArchR:::.validInput(input = imputeWeights, name = "imputeWeights", valid = c("list", "null"))
  ArchR:::.validInput(input = pal, name = "pal", valid = c("character", "null"))
  ArchR:::.validInput(input = size, name = "size", valid = c("numeric"))
  ArchR:::.validInput(input = rastr, name = "rastr", valid = c("boolean"))
  ArchR:::.validInput(input = quantCut, name = "quantCut", valid = c("numeric"))
  ArchR:::.validInput(input = quantHex, name = "quantHex", valid = c("numeric"))
  ArchR:::.validInput(input = discreteSet, name = "discreteSet", valid = c("character", "null"))
  ArchR:::.validInput(input = continuousSet, name = "continuousSet", valid = c("character", "null"))
  ArchR:::.validInput(input = randomize, name = "randomize", valid = c("boolean"))
  ArchR:::.validInput(input = keepAxis, name = "keepAxis", valid = c("boolean"))
  ArchR:::.validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  ArchR:::.validInput(input = addArrow, name = "addArrow", valid = c("boolean"))
  ArchR:::.validInput(input = plotAs, name = "plotAs", valid = c("character", "null"))
  ArchR:::.validInput(input = smoothWindow, name = "smoothWindow", valid = c("integer"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))

  ArchR:::.requirePackage("ggplot2", source = "cran")

  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "Input-Parameters", logFile=logFile)

  #Make Sure ColorBy is valid!
  if(length(colorBy) > 1){
    stop("colorBy must be of length 1!")
  }
  allColorBy <-  c("colData", "cellColData", ArchR:::.availableArrays(getArrowFiles(ArchRProj)))
  if(tolower(colorBy) %ni% tolower(allColorBy)){
    stop("colorBy (",colorBy,") must be one of the following :\n", paste0(allColorBy, sep=", "))
  }
  colorBy <- allColorBy[match(tolower(colorBy), tolower(allColorBy))]

  ##############################
  # Plot Helpers
  ##############################
  .summarizeHex <- function(x = NULL){
    quantile(x, quantHex, na.rm = TRUE)
  }

  ##############################
  # Get Trajectory
  ##############################
  dfT <- getCellColData(ArchRProj, select = trajectory)
  idxRemove <- which(is.na(dfT[,1]))
  ArchR:::.logThis(dfT, "dfT", logFile = logFile)

  ##############################
  # Get Embedding
  ##############################
  df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)
  ArchR:::.logThis(df, "embedding", logFile = logFile)
  dfT <- cbind(df, dfT[rownames(df),])
  colnames(dfT) <- c("x", "y", "PseudoTime")

  #Parameters
  plotParams <- list(...)
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- paste0(embedding, " of ", stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,1])
  plotParams$baseSize <- baseSize

  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){

    plotParams$color <- as.vector(getCellColData(ArchRProj, select = name, drop = FALSE)[rownames(df), 1])
    plotParams$discrete <- ArchR:::.isDiscrete(plotParams$color)
    plotParams$continuousSet <- "horizonExtra"
    plotParams$discreteSet <- "stallion"
    plotParams$title <- paste(plotParams$title, " colored by\ncolData : ", name)
    if(is.null(plotAs)){
      plotAs <- "hexplot"
    }

  }else{

    units <- tryCatch({
        ArchR:::.h5read(getArrowFiles(ArchRProj)[1], paste0(colorBy, "/Info/Units"))[1]
      },error=function(e){
        "values"
    })

    plotParams$continuousSet <- "solarExtra"

    if(is.null(log2Norm) & tolower(colorBy) == "genescorematrix"){
      log2Norm <- TRUE
      plotParams$continuousSet <- "horizonExtra"
    }

    if(is.null(log2Norm)){
      log2Norm <- FALSE
    }

    plotParams$color <- ArchR:::.getMatrixValues(ArchRProj, name = paste0(colorMatrixAssay, ":", name), matrixName = colorBy, log2Norm = FALSE)

    if(!all(rownames(df) %in% colnames(plotParams$color))){
      ArchR:::.logMessage("Not all cells in embedding are present in feature matrix. This may be due to using a custom embedding.", logFile = logFile)
      stop("Not all cells in embedding are present in feature matrix. This may be due to using a custom embedding.")
    }

    plotParams$color <- plotParams$color[, rownames(df), drop = FALSE]

    ArchR:::.logThis(plotParams$color, "colorMat-Before-Impute", logFile = logFile)

    plotParams$discrete <- FALSE
    plotParams$title <- sprintf("%s colored by\n%s : %s", plotParams$title, colorBy, name)
    if(is.null(plotAs)){
      plotAs <- "hexplot"
    }
    if(plotAs=="hexplot"){
      plotParams$fun <- .summarizeHex
    }

  }

  #Additional Params!
  plotParams$xlabel <- gsub("_", " ",stringr::str_split(colnames(df)[1],pattern="#",simplify=TRUE)[,2])
  plotParams$ylabel <- gsub("_", " ",stringr::str_split(colnames(df)[2],pattern="#",simplify=TRUE)[,2])

  if(!is.null(continuousSet)){
    plotParams$continuousSet <- continuousSet
  }
  if(!is.null(continuousSet)){
    plotParams$discreteSet <- discreteSet
  }
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize

  if(plotParams$discrete){
    plotParams$color <- paste0(plotParams$color)
  }

  if(!plotParams$discrete){

    if(!is.null(imputeWeights)){
      message("Imputing Matrix")
      mat <- matrix(as.vector(plotParams$color), nrow = 1)
      colnames(mat) <- rownames(df)
      plotParams$color <- imputeMatrix(mat = mat, imputeWeights = imputeWeights, logFile = logFile)
      ArchR:::.logThis(plotParams$color, "colorMat-After-Impute", logFile = logFile)
    }

    plotParams$color <- as.vector(plotParams$color)

    if(name != trajectory){
      plotParams$color <- ArchR:::.quantileCut(plotParams$color, min(quantCut), max(quantCut))
    }

    if(!is.null(log2Norm)){
      if(log2Norm){
        plotParams$color <- log2(plotParams$color + 1)
        plotParams$colorTitle <- paste0("Log2(",units," + 1)")
      }else{
        plotParams$colorTitle <- units
      }
    }

    plotParams$color[idxRemove] <- NA
    plotParams$pal <- paletteContinuous(set = plotParams$continuousSet)

    if(tolower(plotAs) == "hex" | tolower(plotAs) == "hexplot"){
      plotParams$addPoints <- TRUE
      if(is.null(plotParams$bins)){
        plotParams$bins <- 150
      }
      message("Plotting")
      ArchR:::.logThis(plotParams, name = "PlotParams", logFile = logFile)
      out <- do.call(ggHex, plotParams)
    }else{
      message("Plotting")
      ArchR:::.logThis(plotParams, name = "PlotParams", logFile = logFile)
      out <- do.call(ggPoint, plotParams)
    }

  }else{
    message("Plotting")
    ArchR:::.logThis(plotParams, name = "PlotParams", logFile = logFile)
    out <- do.call(ggPoint, plotParams)
  }

  if(!keepAxis){
    out <- out + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  ArchR:::.logMessage("Plotting Trajectory", logFile = logFile)

  #Prep Trajectory Vector
  dfT$value <- plotParams$color
  dfT <- dfT[order(dfT$PseudoTime), ]
  dfT <- dfT[!is.na(dfT$PseudoTime), ]

  ArchR:::.logThis(dfT, "TrajectoryDF", logFile = logFile)

  #Plot Pseudo-Time
  out2 <- ggPoint(
    x = dfT$PseudoTime,
    y = dfT$value,
    color = dfT$PseudoTime,
    discrete = FALSE,
    xlabel = "PseudoTime",
    ylabel = name,
    pal = plotParams$pal,
    ratioYX = 0.5,
    rastr = TRUE
  ) + geom_smooth(color = "black")

  attr(out2, "ratioYX") <- 0.5


  if(addArrow){

    ArchR:::.logMessage("Adding Inferred Arrow Trajectory to Plot", logFile = logFile)

    dfArrow <-  split(dfT, floor(dfT$PseudoTime / 1.01)) %>%
      lapply(colMeans) %>% Reduce("rbind",.) %>% data.frame
    dfArrow$x <- ArchR:::.centerRollMean(dfArrow$x, smoothWindow)
    dfArrow$y <- ArchR:::.centerRollMean(dfArrow$y, smoothWindow)
    dfArrow <- rbind(dfArrow, dfT[nrow(dfT), ,drop = FALSE])

    out <- out + geom_path(
            data = dfArrow, aes(x, y, color=NULL), size= 1,
            arrow = arrow(type = "open", length = unit(0.1, "inches"))
          )
  }

  ArchR:::.endLogging(logFile = logFile)

  list(out, out2)

}

