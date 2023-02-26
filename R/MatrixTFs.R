#' export
Dorothea <- function(species="human"){
  #species human or mouse
  net <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))
  mat <- as.matrix(data@assays$RNA@data)
  # Run wmean
  acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',
                    .mor='mor', times = 100, minsize = 5)

  # Extract norm_wmean and store it in tfswmean in pbmc
  data[['tfswmean']] <- acts %>%
    filter(statistic == 'norm_wmean') %>%
    pivot_wider(id_cols = 'source', names_from = 'condition',
                values_from = 'score') %>%
    column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)

  # Change assay
  DefaultAssay(object = data) <- "tfswmean"

  # Scale the data
  data <- ScaleData(data)
  data@assays$tfswmean@data <- data@assays$tfswmean@scale.data

  return(0)
}
# run dorothea
#https://github.com/saezlab/transcriptutorial/blob/master/scripts/04_TranscriptionFactor_activity_with_Dorothea.md


#https://saezlab.github.io/decoupleR/articles/tf_bk.html
#https://saezlab.github.io/decoupleR/articles/tf_sc.html


