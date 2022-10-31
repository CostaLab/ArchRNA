.get10XType <- function(files){
  if(all(tools::file_ext(files)=="h5" | tools::file_ext(files)=="H5")){
    return("10XH5")
  }else if(all(tools::file_ext(files)=="")){
    return("10X")
  }
  return("UNKNOWN")
}

RnaArchRProject <- function(inputFiles,
                            sampleNames,
                            outputDirectory=NULL,
                            file_type = c("10X", "10XH5"),
                            matrixName="GeneExpressionMatrix"

){
  if(is.null(outputDirectory)){
    ArchR:::.logError("outputDirectory is NULL, please specify!")
  }
  for(i in seq_along(inputFiles)){
    typename = .get10XType(inputFiles)
    if(typename %ni% file_type){
      stop("Error: Unknown 10 input type, options, 10X or 10XH5!")
    }
    seRNA <- import10xFeatureMatrix_(inputFiles[i], typename, sampleNames[i])
    createRNAarrow(seRNA, sampleNames[i])
    addGeneExpressionMatrix_(input = paste0(sampleNames[i],".arrow"),
                            seRNA = seRNA,
                            chromSizes = getChromSizes(),
                            matrixName=matrixName
    )
  }

  proj <- ArchRProject(ArrowFiles = paste0(sampleNames, ".arrow"),
                       outputDirectory = outputDirectory,
                       copyArrows = TRUE)
  gc()
  return(proj)

}

Matrix2ArchRArrow <- function(mat,
                              name=NULL,
                              meta.data=NULL,
                              matrixName="GeneExpressionMatrix",
                              addHashtag=TRUE
){


  if(!addHashtag){#use exisiting name
    name = unique(stringr::str_split(colnames(mat), pattern = "#", simplify=TRUE)[,1])
    allname = stringr::str_split(colnames(mat), pattern = "#", simplify=TRUE)[,1]
  }else{
    allname = rep(name, ncol(mat))
  }
  spMat <- tryCatch({
    spMat <- as(mat, "dgCMatrix")
    }, error = function(e){
    stop("Erorr: Input mat is not valid matrix!")
  })
  if(is.null(rownames(spMat))){
    stop("Erorr: Rownames of the matrix is NULL!")
  }
  if(is.null(rownames(spMat))){
    stop("Erorr: Colnames of the matrix is NULL!")
  }
  if(!is.null(meta.data)){
    if(!all(rownames(meta.data) == colnames(spMat))){
      stop("Error: rownames of meta data is not consistent with rownames of mat")
    }
    ArchR:::.validInput(input = meta.data, name = "meta.data", valid = c("data.frame", "DataFrame"))
  }

  for(x in seq_along(name)){
    idx = which(allname %in% name[x])
    seRNA = .importSparseFM_(spMat[, idx, drop=F], name=name[x])
      if(!is.null(meta.data)){
        if(addHashtag){
          rownames(meta.data) <- paste0(name[x], "#", rownames(meta.data))
        }
        colData(seRNA) <- DataFrame(meta.data[idx, ,drop=F])
      }
      createRNAarrow(seRNA, name[x])
      addGeneExpressionMatrix_(input = paste0(name[x],".arrow"),
                               seRNA = seRNA,
                               chromSizes = getChromSizes(),
                               matrixName = matrixName
      )
  }
  gc()
  return(paste0(name,".arrow"))
}

Matrix2ArchRProject <- function(mat,
                                name=NULL,
                                meta.data=NULL,
                                outputDirectory=NULL,
                                matrixName="GeneExpressionMatrix",
                                addHashtag=TRUE
){

  if(is.null(outputDirectory)){
    ArchR:::.logError("outputDirectory is NULL, please specify!")
  }
  Matrix2ArchRArrow(mat=mat,
                    name=name,
                    meta.data=meta.data,
                    matrixName=matrixName,
                    addHashtag=addHashtag)

  proj <- ArchRProject(ArrowFiles = paste0(name, ".arrow"),
                       outputDirectory = outputDirectory,
                       copyArrows = TRUE)
  gc()
  return(proj)
}


#Seurat2ArchRProject <- function(SeuratObject,
#                                outputDirectory=NULL

.loadRangeseRNA <- function(h5, name){
    seRNA <- import10xFeatureMatrix_(
      input = c(h5),
      file_type = "10XH5",
      names = c(name)
    )
    if(is.null(rowRanges(seRNA))){
      geneAnnotation <- getArchRGenome(geneAnnotation = TRUE)
      if(is.null(geneAnnotation)){
        stop("getGenomeAnnotation : ArchRPRoj is NULL and there is no genome set with addArchRGenome!")
      }
      gtf <- geneAnnotation$gene
      shared_genes <- intersect(gtf$symbol, rownames(seRNA))
      seRNA <- seRNA[shared_genes, ]
      gtf <- gtf[gtf$symbol %in% shared_genes]

      gtfMatch <- gtf[na.omit(match(rownames(seRNA), gtf$symbol))]
      names(gtfMatch) <- rownames(seRNA)
      rowRanges(seRNA) <- gtfMatch
    }
    return(seRNA)
}

createRNAarrow <- function(seRNA=NULL,
                           name = NULL
                          ){

  bcs <- as.character(colnames(seRNA)) ## need split the bcs by the sample name and add to the
  bcs = stringr::str_split(bcs, pattern = "#", simplify=TRUE)[,2]
  ## for now only missing the metadata rownames
  outArrow=paste0(name, ".arrow")
  o <- ArchR:::.suppressAll(file.remove(outArrow))
  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = outArrow, name = "ArchRVersion")
  o <- h5createGroup(outArrow, paste0("Metadata"))
  o <- h5write(obj = bcs, file = outArrow, name = "Metadata/CellNames")
  o <- h5write(obj = name, file = outArrow, name = "Metadata/Sample")
  o <- h5write(obj = paste0(Sys.Date()), file = outArrow, name = "Metadata/Date")
  o <- h5closeAll()
  gc()
  return(outArrow)
}

addGeneExpressionMatrix_ <- function(
  input = NULL,
  seRNA = NULL,
  chromSizes = getChromSizes(input),
  matrixName = "GeneExpressionMatrix",
  excludeChr = c("chrM", "chrY"),
  scaleTo = 10000,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  strictMatch = FALSE,
  force = TRUE,
  logFile = createLogFile("addGeneExpressionMatrix")
  ){

  if(inherits(input, "ArchRProject")){
    ArrowFiles <- getArrowFiles(input)
    allCells <- rownames(getCellColData(input))
    outDir <- getOutputDirectory(input)
    if(is.null(chromSizes)){
      chromSizes <- getChromSizes(input)
    }
  }else if(inherits(input, "character")){
    outDir <- ""
    ArrowFiles <- input
    allCells <- NULL
    if(is.null(chromSizes)){
      chromSizes <- getChromSizes()
    }
  }else{
    stop("Error Unrecognized Input!")
  }
  if(!all(file.exists(ArrowFiles))){
    stop("Error Input Arrow Files do not all exist!")
  }

  cellsInArrows <- unlist(lapply(ArrowFiles, ArchR:::.availableCells), use.names=FALSE)
  if(!is.null(allCells)){
    cellsInArrows <- allCells
  }


  splitCells <- split(cellsInArrows, stringr::str_split(cellsInArrows, pattern = "#", simplify=TRUE)[,1])
  overlapPerSample <- unlist(lapply(splitCells, function(x) sum(x %in% colnames(seRNA))))
  ArchR:::.logMessage("Overlap Per Sample w/ scATAC : ", paste(paste(names(overlapPerSample), round(overlapPerSample,3), sep = "="), collapse=","), logFile = logFile, verbose = TRUE)

  #Get QC Info
  assay(seRNA) <- Matrix::Matrix(assay(seRNA), sparse=TRUE)
  nUMI <- Matrix::colSums(assay(seRNA))
  mb <- assay(seRNA)
  mb@x[mb@x > 0] <- 1
  nGenes <- Matrix::colSums(mb)
  rm(mb)
  MitoRatio <- Matrix::colSums(assay(seRNA)[grep("^MT", rownames(assay(seRNA))),]) / nUMI
  RiboRatio <- Matrix::colSums(assay(seRNA)[grep("^RP", rownames(assay(seRNA))),]) / nUMI
  qcInfo <- DataFrame(nUMI = nUMI, nGenes = nGenes, MitoRatio = MitoRatio, RiboRatio = RiboRatio)
  colnames(qcInfo) <- paste0("Gex_", colnames(qcInfo))

  #Dedup
  idxDup <- which(rownames(seRNA) %in% rownames(seRNA[duplicated(rownames(seRNA))]))
  names(idxDup) <- rownames(seRNA)[idxDup]
  if(length(idxDup) > 0){
    dupOrder <- idxDup[order(Matrix::rowSums(assay(seRNA[idxDup])),decreasing=TRUE)]
    dupOrder <- dupOrder[!duplicated(names(dupOrder))]
    seRNA <- seRNA[-as.vector(idxDup[idxDup %ni% dupOrder])]
  }

  #Add Index To RNA Ranges
  features <- rowRanges(seRNA)
  features <- Reduce("c",features)

   #Add args to list
  args <- mget(names(formals()), sys.frame(sys.nframe()))#as.list(match.call())
  #args <- list()
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$X <- seq_along(ArrowFiles)
  args$FUN <- addGeneExpressionMat_
  args$matrixName <- matrixName
  args$registryDir <- file.path(outDir, "addGeneExpressionMatRegistry")
  args$qcInfo <- qcInfo
  args$seRNA <- seRNA

  #Remove Input from args
  args$input <- NULL
  args$chromSizes <- NULL
  args$strictMatch <- NULL

  #Run With Parallel or lapply
  outList <- ArchR:::.batchlapply(args)


    #Return Output
    if(inherits(input, "ArchRProject")){

      qcInfo <- qcInfo[rownames(qcInfo) %in% input$cellNames, ]

      for(i in seq_len(ncol(qcInfo))){
    input <- addCellColData(
      ArchRProj = input,
      data = as.vector(qcInfo[,i]),
      name = paste0(colnames(qcInfo)[i]),
      cells = paste0(rownames(qcInfo)),
      force = force
    )
      }

      return(input)

    }else{

      return(unlist(outList))
    }
}

addGeneExpressionMat_ <- function(
  i = NULL,
  ArrowFiles = NULL,
  matrixName = "GeneExpressionMatrix",
  seRNA = NULL,
  qcInfo = NULL,
  excludeChr = NULL,
  scaleTo = NULL,
  cellNames = NULL,
  allCells = NULL,
  tstart = NULL,
  subThreads = 1,
  force = FALSE,
  verbose = TRUE,
  logFile = NULL
  ){

  ArrowFile <- ArrowFiles[i]
  sampleName <- ArchR:::.sampleName(ArrowFile)

  #Check
  o <- h5closeAll()
  o <- ArchR:::.createArrowGroup(ArrowFile = ArrowFile, group = matrixName, force = force, logFile = logFile)

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #Get all cell ids before constructing matrix
  if(is.null(cellNames)){
    cellNames <- ArchR:::.availableCells(ArrowFile)
  }

  if(!is.null(allCells)){
    cellNames <- cellNames[cellNames %in% allCells]
  }

  #Identify Overlapping Cells
  cellNames <- cellNames[cellNames %in% colnames(seRNA)]
  seRNA <- seRNA[, cellNames]

  dfParams <- data.frame( #####################################
    slot = "data",
    scaleTo = scaleTo,
    exclude = excludeChr,
    stringsAsFactors = FALSE
  )

  featureDF1 <- data.frame(
  seqnames = "data",
  idx = seq(1, nrow(seRNA)),
  name = rownames(seRNA),
  stringsAsFactors = FALSE
  )
  featureDF2 <- data.frame(
  seqnames = "counts",
  idx = seq(1, nrow(seRNA)),
  name = rownames(seRNA),
  stringsAsFactors = FALSE
  )

  featureDF <- rbind(featureDF1, featureDF2)
  #featureDF <- featureDF1




  ######################################
  # Initialize SP Mat Group
  ######################################

  o <- ArchR:::.initializeMat(
    ArrowFile = ArrowFile,
    Group = matrixName,
    Class = "Assays",
    Units = c("data", "counts"),
    cellNames = colnames(seRNA),
    #params = "GeneExpressionMatrix",
    params = dfParams,
    featureDF = featureDF,
    force = TRUE
  )



  ######################################
  # Normalize and Insert Log2 Normalized Counts
  ######################################

  #assay(seRNA) <- ArchR:::.normalizeCols(assay(seRNA), scaleTo = scaleTo)
  assays(seRNA)[["data"]] = ArchR:::.normalizeCols(assay(seRNA), scaleTo = scaleTo)

  o <- tryCatch({
    o <- h5closeAll()
    matz <- assays(seRNA)[["data"]]
    matc <- assays(seRNA)[["counts"]]
    ArchR:::.logDiffTime(sprintf("Adding data to GeneExpressionMatrix !"), tstart, verbose = verbose, logFile = logFile)

      #Write sparseMatrix to Arrow File!
      o <- addMatToArrow_(
        mat = matz,
        ArrowFile = ArrowFile,
        Group = paste0(matrixName, "/data"),
        binarize = FALSE,
        addColSums = TRUE,
        addRowSums = TRUE,
        addRowVarsLog2 = TRUE #add for integration analyses
      )

    ArchR:::.logDiffTime(sprintf("Adding counts to GeneExpressionMatrix !"), tstart, verbose = verbose, logFile = logFile)
      o <- addMatToArrow_(
        mat = matc,
        ArrowFile = ArrowFile,
        Group = paste0(matrixName, "/counts"),
        binarize = FALSE,
        addColSums = TRUE,
        addRowSums = TRUE,
        addRowVarsLog2 = FALSE #add for integration analyses
      )
      gc()
    }, error = function(e){

      errorList <- list(
        ArrowFile = ArrowFile,
        chr = chr,
        mat = if(exists("matz", inherits = FALSE)) matz else "matz"
      )

      ArchR:::.logError(e, fn = ".addGeneExpressionMat AddToArrow", info = sampleName, errorList = errorList, logFile = logFile)

    })

  #Add Info To Arrow Files
  allCells <- ArchR:::.availableCells(ArrowFile, passQC = FALSE)

  qcInfoi <- qcInfo[rownames(qcInfo) %in% colnames(seRNA), ]

  for(i in seq_len(ncol(qcInfo))){

    infoi <- rep(-1, length(allCells))
    names(infoi) <- allCells
    infoi[rownames(qcInfoi)] <- qcInfoi[,i]

    o <- h5closeAll()
    h5write(infoi, file = ArrowFile, paste0("Metadata/", colnames(qcInfoi)[i]))
    o <- h5closeAll()

  }
  ArrowFile
}

addMatToArrow_ <- function(
  mat = NULL,
  ArrowFile = NULL,
  Group = NULL,
  binarize = FALSE,
  addRowSums = FALSE,
  addColSums = FALSE,
  addRowMeans = FALSE,
  addRowVars = FALSE,
  addRowVarsLog2 = FALSE,
  logFile = NULL
  ){

  stopifnot(inherits(mat, "dgCMatrix"))

  #Create Group
  o <- h5closeAll()
  o <- h5createGroup(ArrowFile, Group)

  #Convert Columns to Rle
  j <- Rle(findInterval(seq(mat@x)-1,mat@p[-1]) + 1)

  #Info
  lengthRle <- length(j@lengths)
  lengthI <- length(mat@i)

  #Create Data Set
  o <- ArchR:::.suppressAll(h5createDataset(ArrowFile, paste0(Group,"/i"), storage.mode = "integer",
    dims = c(lengthI, 1), level = 0))

  o <- ArchR:::.suppressAll(h5createDataset(ArrowFile, paste0(Group,"/jLengths"), storage.mode = "integer",
    dims = c(lengthRle, 1), level = 0))

  o <- ArchR:::.suppressAll(h5createDataset(ArrowFile, paste0(Group,"/jValues"), storage.mode = "integer",
    dims = c(lengthRle, 1), level = 0))

  #Write Data Set
  o <- ArchR:::.suppressAll(h5write(obj = mat@i + 1, file = ArrowFile, name = paste0(Group,"/i")))
  o <- ArchR:::.suppressAll(h5write(obj = j@lengths, file = ArrowFile, name = paste0(Group,"/jLengths")))
  o <- ArchR:::.suppressAll(h5write(obj = j@values, file = ArrowFile, name = paste0(Group,"/jValues")))

  #If binary dont store x
  if(!binarize){

    o <- ArchR:::.suppressAll(h5createDataset(ArrowFile, paste0(Group, "/x"), storage.mode = "double",
      dims = c(lengthI, 1), level = 0))

    o <- ArchR:::.suppressAll(h5write(obj = mat@x, file = ArrowFile, name = paste0(Group, "/x")))

  }else{

    mat@x[mat@x > 0] <- 1

  }

  if(addColSums){
    cS <- Matrix::colSums(mat)
    o <-ArchR:::.suppressAll(h5createDataset(ArrowFile, paste0(Group, "/colSums"), storage.mode = "double",
      dims = c(ncol(mat), 1), level = 0))
    o <- ArchR:::.suppressAll(h5write(obj = cS, file = ArrowFile, name = paste0(Group, "/colSums")))

  }

  if(addRowSums){
    rS <- Matrix::rowSums(mat)
    o <-ArchR:::.suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowSums"), storage.mode = "double",
      dims = c(nrow(mat), 1), level = 0))
    o <- ArchR:::.suppressAll(h5write(obj = rS, file = ArrowFile, name = paste0(Group, "/rowSums")))

  }

  if(addRowMeans){
    rM <- Matrix::rowMeans(mat)
    o <- ArchR:::.suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowMeans"), storage.mode = "double",
      dims = c(nrow(mat), 1), level = 0))
    o <- ArchR:::.suppressAll(h5write(obj = rM, file = ArrowFile, name = paste0(Group, "/rowMeans")))

  }

  if(addRowVars){
    if(!addRowMeans){
      rM <- Matrix::rowMeans(mat)
    }
    rV <- ArchR:::computeSparseRowVariances(mat@i + 1, mat@x, rM, n = ncol(mat))
    o <- ArchR:::.suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowVars"), storage.mode = "double",
      dims = c(nrow(mat), 1), level = 0))
    o <- ArchR:::.suppressAll(h5write(obj = rV, file = ArrowFile, name = paste0(Group, "/rowVars")))

  }

  if(addRowVarsLog2){

    mat@x <- log2(mat@x + 1) #log-normalize
    rM <- Matrix::rowMeans(mat)
    idx <- seq_along(rM)
    idxSplit <- ArchR:::.splitEvery(idx, 2000)

    #Make sure not too much memory so split into 2000 gene chunks
    #Check this doesnt cause memory mapping issues!
    rV <- lapply(seq_along(idxSplit), function(x){
      idxX <- idxSplit[[x]]
      matx <- mat[idxX, , drop = FALSE]
      ArchR:::computeSparseRowVariances(matx@i + 1, matx@x, rM[idxX], n = ncol(mat))
    }) %>% unlist

    #Have to write rowMeansLog2 as well
    o <- ArchR:::.suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowMeansLog2"), storage.mode = "double",
      dims = c(nrow(mat), 1), level = 0))
    o <- ArchR:::.suppressAll(h5write(obj = rM, file = ArrowFile, name = paste0(Group, "/rowMeansLog2")))

    #Write rowVarsLog2
    o <- ArchR:::.suppressAll(h5createDataset(ArrowFile, paste0(Group, "/rowVarsLog2"), storage.mode = "double",
      dims = c(nrow(mat), 1), level = 0))
    o <- ArchR:::.suppressAll(h5write(obj = rV, file = ArrowFile, name = paste0(Group, "/rowVarsLog2")))
  }

  #Clean Up Memorys
  rm(j,mat)
  gc()

  o <- h5closeAll()
  gc()
  return(0)

}
