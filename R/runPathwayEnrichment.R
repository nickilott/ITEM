# ' Run pathways analysis for each module from a given network
#'
#' Use GOseq to preform enrichment analysis on each module in a
#' given network
#' @param background background set of genes
#' @param foreground foreground set of genes
#' @param input.type the type of input id - in supportedGeneIDs() from goseq
#' @param organism.db organism of input - in supportedGenomes() from goseq 
#' @param pathway.db database used by GOSeq to perform
#' enrichment testing: default is GO:BP
#' @param bias.data median length of transcripts for each gene in data, data frame mapping gene id to length
#' @return 
#' @import goseq futile.logger
#' @examples
#' 
#' @export

runPathwayEnrichment <- function(background, foreground, input.type="ensGene", organism.db="hg38", pathway.db="GO:BP", bias.data=NULL){

    foreground <- as.character(foreground)
    background <- as.character(background)

    flog.info("creating dataset to test")
    de <- rep(1, length(foreground))
    not.de <- rep(0, length(background))
    to.test <- append(de, not.de)
    gene.ids <- append(foreground, background)

    names(to.test) <- gene.ids

    if (!(is.null(bias.data))){
        bias.data <- bias.data[gene.ids,]$length
	}

    flog.info("creating null distribution")
    pwf <- nullp(to.test, organism.db, input.type, bias.data=bias.data)

    flog.info("getting enrichment results")
    goseq.res <- goseq_res <- goseq(pwf, organism.db, input.type, test.cats=pathway.db)

    return(goseq.res)

}