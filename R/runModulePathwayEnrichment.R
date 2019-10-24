#' Run pathways analysis for each module from a given network
#'
#' Use GOseq to preform enrichment analysis on each module in a
#' given network
#' ITEMDataSet object defining module assignment
#' @param input.type the type of input id e.g. ensembl
#' @param organism.db organism of input
#' @param pathway.db database used by GOSeq to perform
#' enrichment testing
#' @param bias.data median length of transcripts for each gene in data, data frame mapping gene id to length
#' @return 
#' @examples
#' 
#' @export

runModulePathwayEnrichment <- function(ITEMDataSet, input.type="ensGene", organism.db="hg38", pathway.db="GO:BP", bias.data=NULL){

    tissue1.modules <- ITEMDataSet[[10]]
    name1 <- ITEMDataSet$ref
    tissue2.modules <- ITEMDataSet[[11]]
    name2 <- ITEMDataSet$comparison

    # create background sets
    tissue1.background <- tissue1.modules[tissue1.modules$Module != 0,]$Gene

    tissue2.background <- tissue2.modules[tissue2.modules$Module != 0,]$Gene

    # rename modules
    tissue1.modules$Module <- paste0(name1, ".", tissue1.modules$Module)
    tissue2.modules$Module <- paste0(name2, ".", tissue2.modules$Module)

    # bind rows here
    modules.to.test <- data.frame(rbind(tissue1.modules, tissue2.modules))

    # enrichment is performed amongst genes that have been placed in a module i.e. remove those
    # that are in module 0 from the background
    modules.to.test <- modules.to.test[grep("\\.0", modules.to.test$Module, invert=TRUE),]    

    # iterate over modules and return enrichment results
    enrichment.results <- list()
    modules <- unique(modules.to.test$Module)

    for (i in 1:length(modules)){
        module <- modules[i]
        flog.info(paste0("runnning pathway enrichment analysis for module ", module))
        if (name1 %in% module){
	    background <- tissue1.background
	}
	else{
	    background <- tissue2.background
        }
	foreground <- modules.to.test[modules.to.test$Module == module,]$Gene
        background <- setdiff(background, foreground)

        enrichment.result <- runPathwayEnrichment(background,
	                                          foreground,
						  input.type=input.type,
						  organism.db=organism.db,
						  pathway.db=pathway.db,
						  bias.data=bias.data)

        enrichment.results[[i]] <- enrichment.result
    }
    names(enrichment.results) <- modules
    return(enrichment.results)
}