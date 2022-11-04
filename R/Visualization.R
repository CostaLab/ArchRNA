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
FeaturePlotS <- function(project, features=NULL, useMatrix=NULL, assay=NULL, reduction=NULL, ...){
  if(is.null(useMatrix) & all(features %ni% colnames(project@cellColData))){
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
  p <- Seurat::FeaturePlot(object, reduction=reduction, features=features,  ...)
  return(p)
}

#' @export
RidgePlotS <-function(project, useMatrix=NULL, assay=NULL, ...){
  if(is.null(useMatrix)){
    useMatrix = getAvailableMatrices(project)[1]
    warning(glue::glue("useMatrix is NULL, use first matrix {useMatrix}!"))
  }
  object <- PartialSeurat(project, useMatrix=useMatrix, assay=assay)
  p <- Seurat::RidgePlot(object,assay="RNA", slot=assay, ...)
  return(p)
}

#' @export
DotPlotS <- function(project, useMatrix=NULL, assay=NULL, ...){
  if(is.null(useMatrix)){
    useMatrix = getAvailableMatrices(project)[1]
    warning(glue::glue("useMatrix is NULL, use first matrix {useMatrix}!"))
  }
  object <- PartialSeurat(project, useMatrix=useMatrix, assay=assay)
  p <- Seurat::DotPlot(object, assay="RNA", ...)
  return(p)
}

### can add mean or median

#' @export
VlnPlotS <- function(project, features=c(), useMatrix=NULL, assay=NULL, ...){

  ### if vln in meta.data, there's no need to retrieve count matrix
  if(is.null(useMatrix) & all(features %ni% colnames(project@cellColData))){
    useMatrix = getAvailableMatrices(project)[1]
    warning(glue::glue("useMatrix is NULL, use first matrix {useMatrix}!"))
  }
  object <- PartialSeurat(project, useMatrix=useMatrix, assay=assay)
  p <- Seurat::VlnPlot(object, features=features, assay="RNA", slot=assay, ...)
  rm(object)
  gc()
  return(p)
}

#' @export
scProportionPlotS <- function(project,
                              clusterName="Cluster_0.5",
                              condition="Sample",
                              pair=c("a", "b"),
                              logFile = ArchR::createLogFile("scProportionPlotS")
){

#TODO:
  #add checking the input parameters's validation

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

#' @export
FlipBarPlotS <- function(DF, top_n=10, name="", sort.val="desc", x.text.angle=90){
    #Cluster,
    #DF <- DF %>% dplyr::filter(Cluster == Cluster)
    dfm <- as.data.frame(DF) %>% dplyr::slice_max(Log2FC, n=top_n) # %>% dplyr::mutate(updn = ifelse(Log2FC>0, "up", "down"))
    p <- ggbarplot(dfm,
              x = "name",
              y = "Log2FC",
              fill = "lightblue",         # change fill color by mpg_level
              color = "white",            # Set bar border colors to white
              palette = "jco",            # jco journal color palett. see ?ggpar
              sort.val = sort.val,           # Sort the value in ascending order
              sort.by.groups = FALSE,     # Don't sort inside each group
              x.text.angle = 90,          # Rotate vertically x axis texts
              ylab = "Log2FC",
              xlab = FALSE,
              legend.title = glue::glue("{name}")
    )
    p
}

#' @export
GroupedBarPlotS <- function(project, Cluster="Cluster", condition="Sample", cells=NULL){
  #each group is
  #proportions in each sample or condition in each cluster
  #nCK124 = length(which(scrna$name == "CK124"))
  CK124_nCells <- table(scrna$clusters[scrna$name=="CK124"])
  clusters_seq = sort(unique(scrna$clusters))
  for(i in clusters_seq){
    if(is.na(CK124_nCells[as.character(i)])){
        CK124_nCells[as.character(i)] = 0
    }
  }
  CK124_nCells=CK124_nCells[order(as.numeric(names(CK124_nCells)))]

  #CK124_nCells['23'] = 0
  CK124_total_cells <- sum(CK124_nCells)
  CK124_Prop <- round(100.0*CK124_nCells/CK124_total_cells, 3)

  #CK124_Prop['23'] = 0
  CK124_tbl <- as.table(cbind(CK124_nCells, CK124_Prop))

  CK142_nCells <- table(scrna@meta.data$clusters[scrna$name=="CK142"])
  for(i in 1:clusters_seq){
    if(is.na(CK142_nCells[as.character(i)])){
      CK142_nCells[as.character(i)] = 0
    }
  }
  CK142_nCells=CK142_nCells[order(as.numeric(names(CK142_nCells)))]
  CK142_total_cells <- sum(CK142_nCells)
  CK142_Prop <- round(100.0*CK142_nCells/CK142_total_cells, 3)

  CK142_tbl <- as.table(cbind(CK142_nCells, CK142_Prop))

  tbl <- cbind(CK124_tbl, CK142_tbl)

  suppressPackageStartupMessages(library(knitr))
  suppressPackageStartupMessages(library(kableExtra))
  suppressPackageStartupMessages(library(formattable))



  #-------
  library(data.table)
  ptbl <- tbl[, c("CK124_Prop", "CK142_Prop")]/100
  colnames(ptbl) <- c("CK124","CK142")
  mtbl <- melt(t(ptbl))
  colnames(mtbl) <- c("name", "cluster", "prop")
  df = data.frame(mtbl)
  ggplot(df, aes(x = name, y=prop, fill=name)) +
                  facet_wrap(~cluster, scales = "free", ncol=4)  +
                  scale_y_continuous(labels=scales::percent) +
                  geom_bar(stat = "identity") +
                  scale_fill_manual(values=c("#e41a1c","#4daf4a"))
                  #scale_fill_manual(values=c("CK124" = "#e41a1c", "Trans9" = "#3
}



#' @export
PropBarPlotS <- function(project, Cluster="Cluster", condition="Sample"){

  dfm = as.data.frame.table(table(as.vector(project@cellColData[, condition]), as.vector(project@cellColData[, Cluster]))) %>%
                                                dplyr::mutate(proportion=Freq/sum(Freq)) %>%
                                                dplyr::select(-c(Freq))
  colnames(dfm) <- c(condition, Cluster, "proportion")

  ggplot(dfm) +
          aes(x=!!sym(Cluster), y=proportion,  fill = !!sym(condition)) +
          geom_bar(position = "fill", stat = "identity") +
          scale_y_continuous(labels=scales::percent)
}


##
DensityPlotS<- function(project){
library(Nebulosa)
}





RidgesMatrixPlotS <- function(project, features){
}

#' @export
VlnMatrixPlotS <- function(project, Cluster="Cluster", features=NULL, useMatrix="GeneExpressionMatrix", show_cluster=T,leftmost_width=1.4){
  ps <- VlnPlotS(project,  useMatrix=useMatrix, features = features, pt.size = 0, ncol = 4, group.by = Cluster, combine = FALSE)

  len = length(ps)
  p1 <- ps[[1]] + coord_flip() + NoLegend() + scale_fill_discrete(guide='none') +
          theme_bw() +
          theme(plot.title = element_text(angle = 90),
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                #axis.text.y=element_blank(),
                axis.title.y=element_blank(),
                #axis.ticks.y=element_blank(),
                plot.margin = unit(c(0, 0, 0, 0), "cm"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(colour = "black", size=1))
  ps <- lapply(ps[2:len], function(x) x + coord_flip() + NoLegend() + scale_fill_discrete(guide='none') +
          theme_bw() +
          theme(plot.title = element_text(angle = 90),
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.y=element_blank(),
                axis.title.y=element_blank(),
                #axis.ticks.y=element_blank(),
                plot.margin = unit(c(0, 0, 0, 0), "cm"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(colour = "black", size=1)))


  p2 = cowplot::plot_grid(plotlist = ps, nrow=1, align = "h")
  if(show_cluster){
    return(cowplot::plot_grid(plotlist=c(list(p1), ps), align='h', nrow=1, rel_widths=c(leftmost_width, rep(1, len-1))))
  }else{
    return(cowplot::plot_grid(plotlist=c(list(p1+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())), ps),
                                align='h',
                                nrow=1,
                                rel_widths=c(1, rep(1, len-1))))
  }

}

BoxPlotS <- function(){

}

VlnboxPlotS <- function(){

}

## https://github.com/crazyhottommy/scATACutils
## density plots for QC



##
# can have style like https://github.com/harbourlab/PieParty

#' @export
PiePlotS <- function(project, Cluster, condition=NULL, cols=ggsci::pal_igv()(51), round_n=2, ...){

  if (is.null(condition)){
      data = as.data.frame.table(table(project@cellColData[, Cluster])) %>% dplyr::mutate(proportion=round(100.0*Freq/sum(Freq), round_n))
      colnames(data) <- c(Cluster, "cells", "proportion")
      p <- ggpubr::ggpie(data, x="cells", label = paste0(data$proportion, '%'),
                 fill = Cluster, palette = cols, ...) + theme(legend.direction='vertical', legend.position='right')

  }else{
      plist <- list()
      conditions = unique(project@cellColData[, condition])
      for(i in seq_along(conditions)){
        cond = conditions[i]
        meta = project@cellColData
        idx = which(as.vector(meta[, condition]) == cond)
        data = as.data.frame.table(table(meta[idx, Cluster])) %>% dplyr::mutate(proportion=round(100.0*Freq/sum(Freq), round_n))
        colnames(data) <- c(Cluster, "cells", "proportion")
        px <-ggpubr::ggpie(data, x="cells", label = paste0(data$proportion, '%'),
                   fill = Cluster, palette = cols, ...) + ggtitle(cond)

        plist[[cond]] <- px+theme(legend.position='none')
      }
      pl =ggpubr::ggpie(data, x="cells", fill=Cluster, palette=cols) + theme(legend.direction='vertical')
      leg <- ggpubr::get_legend(pl)

      p <- cowplot::plot_grid(plotlist=plist, ncol=1)
      p <- p + ggpubr::as_ggplot(leg) + patchwork::plot_layout(widths = c(6, 1))
  }
  p
}




## each Dot is a piePlots: Li's paper
# x cluster, y features, pie proportions.
PieMatrixPlotS <- function(){
}


## diffusion can show edges with igraph

## kawaii layouts of a graph


## velocity compatible



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

