progenyClusters <- function(project){

  return(project)
}

 #https://saezlab.github.io/decoupleR/articles/pw_sc.html

generate_scrna_progeny <- function(scrna){
  ret_code = 0
  assertthat::assert_that(SPECIES == "Human" | SPECIES == "Mouse")
  Idents(scrna) <- DEFUALT_CLUSTER_NAME
  scrna <- progeny::progeny(scrna, scale=FALSE, organism=SPECIES, top=500, perm=1,return_assay = TRUE)

  da <- DefaultAssay(scrna)
  DefaultAssay(scrna) <-'progeny'
  pws <- rownames(scrna@assays$progeny)
  res <- list()
  Idents(scrna) <- DEFUALT_CLUSTER_NAME
  if (is.factor(scrna@meta.data[, DEFUALT_CLUSTER_NAME])){
    scrna@meta.data[, DEFUALT_CLUSTER_NAME] <- droplevels(scrna@meta.data[, DEFUALT_CLUSTER_NAME])
  }else{
    c_names <-unique(scrna@meta.data[, DEFUALT_CLUSTER_NAME])
    help_sort_func <- ifelse(all.is.numeric(c_names), as.numeric, function(x){x})
    scrna@meta.data[, DEFUALT_CLUSTER_NAME] <- factor(scrna@meta.data[, DEFUALT_CLUSTER_NAME],
                                                      levels=sort(help_sort_func(c_names)))
  }

  for(i in levels(scrna@meta.data[, DEFUALT_CLUSTER_NAME])){
    g <- as.character(scrna@meta.data[, DEFUALT_CLUSTER_NAME])
    g[!(g==i)] <- "others"
    g <- factor(g, levels=c(i, "others"))
    res[[i]] = scran::findMarkers(as.matrix(scrna@assays$progeny@data), g)[[1]]
    res[[i]] <- as.data.frame(res[[i]])
    r <- sapply(pws, function(pw) rcompanion::wilcoxonR(as.vector(scrna@assays$progeny@data[pw,]), g))
    nms <- sapply(stringr::str_split(names(r), "\\."), function(x)x[1])
    names(r) <- nms
    res[[i]][nms, "r"] <- r
    res[[i]] <- res[[i]][nms, ]
  }

  for (cl in names(res)) {
    res[[cl]]$pathway <- rownames(res[[cl]])
    res[[cl]]$CellType <- cl
    colnames(res[[cl]]) <-  c("Top","p.value","FDR", "summary.logFC","logFC","r","pathway","CellType")
  }
  res_df <- do.call("rbind", res)
  res_df$tag <- sapply(res_df$FDR, function(pval) {
      if(pval< 0.001) {
      txt <- "***"
      } else if (pval < 0.01) {
      txt <- "**"
      } else if (pval < 0.05) {
      txt <- "*"
      } else {
      txt <- ""
      }
      return(txt)
  })
  scrna@tools[[glue("progeny_{DEFUALT_CLUSTER_NAME}")]] <- res_df
  DefaultAssay(scrna) <- da
  return(list(scrna, ret_code))
}


generate_scrna_progeny_stage <- function(scrna){

  ret_code = 0
  if("progeny" %ni% names(scrna@assays)){
    scrna <- progeny::progeny(scrna, scale=FALSE, organism=SPECIES, top=500, perm=1,return_assay = TRUE)
  }
  conds <- scrna@tools$meta_order$stage
  m <- combn(conds, 2)
  n = length(m)/2
  lst <- vector("list", n)
  for (i in 1:n){
    lst[[i]] <- m[1:2, i]
  }
  if (is.factor(scrna@meta.data[, DEFUALT_CLUSTER_NAME])){
    scrna@meta.data[, DEFUALT_CLUSTER_NAME] <- droplevels(scrna@meta.data[, DEFUALT_CLUSTER_NAME])
  }else{
    c_names <-unique(scrna@meta.data[, DEFUALT_CLUSTER_NAME])
    help_sort_func <- ifelse(all.is.numeric(c_names), as.numeric, function(x){x})
    scrna@meta.data[, DEFUALT_CLUSTER_NAME] <- factor(scrna@meta.data[, DEFUALT_CLUSTER_NAME],
                                                      levels=sort(help_sort_func(c_names)))
  }

  pws <- rownames(scrna@assays$progeny)
  vs_df_list <- list()
  for(apair in lst){
    vs1 <- apair[1]
    vs2 <- apair[2]

    Idents(scrna) <- 'stage'
    ###
    res <- list()
    for(i in levels(scrna@meta.data[, DEFUALT_CLUSTER_NAME])){
        cells1 <- which(scrna@meta.data[, DEFUALT_CLUSTER_NAME]==i & (scrna@meta.data$stage %in% vs1))
        cells2 <- which(scrna@meta.data[, DEFUALT_CLUSTER_NAME]==i & (scrna@meta.data$stage %in% vs2))
        if( length(cells1) < 2 | length(cells2) < 2){
           next
        }
        a_sub = subset(scrna, cells=c(cells1, cells2))
        g <- as.character(a_sub@meta.data$stage)
        g <- factor(g, levels=c(vs1, vs2))
        res[[i]] = scran::findMarkers(as.matrix(a_sub@assays$progeny@data), g)[[1]]
        res[[i]] <- as.data.frame(res[[i]])
        r <- sapply(pws, function(pw) rcompanion::wilcoxonR(as.vector(a_sub@assays$progeny@data[pw,]), g))
        nms <- sapply(stringr::str_split(names(r), "\\."), function(x)x[1])
        names(r) <- nms
        res[[i]][nms, "r"] <- r
        res[[i]] <- res[[i]][nms, ]

    }

    for (cl in names(res)) {
        res[[cl]]$pathway <- rownames(res[[cl]])
        res[[cl]]$CellType <- cl
        colnames(res[[cl]]) <-  c("Top","p.value","FDR", "summary.logFC","logFC","r","pathway","CellType")

    }

    res_df <- do.call("rbind", res)
#    res_df$FDR <- p.adjust(res_df$p_val, method="fdr")
    res_df$tag <- sapply(res_df$FDR, function(pval) {
        if(pval< 0.001) {
        txt <- "***"
        } else if (pval < 0.01) {
        txt <- "**"
        } else if (pval < 0.05) {
        txt <- "*"
        } else {
        txt <- ""
        }
        return(txt)
    })

    vs_df_list[[glue("{vs1}.vs.{vs2}")]] <- res_df

  }
  scrna@tools[[glue("progeny_stage_{DEFUALT_CLUSTER_NAME}")]] <- vs_df_list
  return(list(scrna, ret_code))
}

