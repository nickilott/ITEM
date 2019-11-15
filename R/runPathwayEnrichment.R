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
#' @param gene2cat data frame specifying gene to category mappings (different to goseq gene2cat specification)
#' It is user-defined and is a data.frame that also contains additional information that is used to produce output
#' compatible with ITEM. i.e it has 5 columns:
#' 1. geneset 2. gene 3. category 4. term (note that this can be the same as category depending on the gene set)
#' @return 
#' @import goseq futile.logger
#' @examples
#' 
#' @export

runPathwayEnrichment <- function(background, foreground, input.type="ensGene", organism.db="hg38", pathway.db="GO:BP", bias.data=NULL, gene2cat=NULL){

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
    pwf <- nullp(to.test, organism.db, input.type, bias.data=bias.data, plot.fit=FALSE)

    flog.info("getting enrichment results")
    if (!(is.null(gene2cat))){
        category.df <- gene2cat
        gene2cat <- data.frame(category=gene2cat[,3], Gene=gene2cat[,2])
    }
    goseq.res <- goseq(pwf, organism.db, input.type, test.cats=pathway.db, gene2cat=gene2cat)

    if (!("term" %in% colnames(goseq.res))){
        cat2term <- as.character(unlist(category.df[,4]))
        print(str(cat2term))
        names(cat2term) <- category.df[,3]
        goseq.res$term <- cat2term[goseq.res$category]
    }
    
    goseq.res$over_represented_fdr <- p.adjust(goseq.res$over_represented_pvalue)
    return(goseq.res)

}