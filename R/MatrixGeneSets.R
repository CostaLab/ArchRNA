genesets <- function(){

}


generate_scrna_MSigDB_geneset <- function(scrna){
    #http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.4/msigdb.v7.4.symbols.gmt
  ret_code = 0
  suppressPackageStartupMessages(require(GSEABase))
  tryCatch(
           {
             HS_Gmt <- getGmt(gzfile(MSigDB_GENESET_HUMAN_GMT_FILE))
             Gmt <- HS_Gmt
             if(SPECIES == "Mouse"){
               all_mappings <- read.csv(file="external/Human2Mouse_mappings.tsv", sep="\t")
               for(idx in 1:length(HS_Gmt)){
                  geneIds <- HS_Gmt[[idx]]@geneIds
                  mm_geneIds <-all_mappings[which(all_mappings$HGNC.symbol %in% geneIds), "MGI.symbol"]
                  Gmt[[idx]]@geneIds <- mm_geneIds
               }
             }
             assertthat::assert_that(all(MSigDB_Geneset_names %in% names(Gmt)))
             subGmt <- Gmt[MSigDB_Geneset_names]

             for(nm in names(subGmt)){
                 geneIds <- subGmt[[nm]]@geneIds
                 geneIds <- intersect(geneIds, rownames(scrna))
                 ## scrna@assays$RNA@data
                 scrna <- AddModuleScore(object = scrna, features = list(geneIds), name = nm, assay="RNA")
                 scrna@meta.data[, nm] <- scrna@meta.data[, paste0(nm, 1)]
                 scrna@meta.data[, paste0(nm, 1)] <- NULL
             }
             scrna@tools$genesets <- names(subGmt)
           },
           error=function(cond) {
             ret_code <<- -1
             logger.error(cond)
             logger.error(traceback())
           },
           finally={
             return(list(scrna, ret_code))
           })
  return(list(scrna, ret_code))
}

