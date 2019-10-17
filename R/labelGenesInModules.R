#' Assign module genes to tissues
#'
#' Assign genes that are in modules to the tissue that they are associated with
#' @param wgcna.net output from WCGNA
#' @param gene.annotation genes annotated to tissue from runDESeq2
#' @return data.frame of genes, their module annotations and tissue assignment
#' @examples
#' labelGenesInModules(wgcna.net, gene.annotation)
#' @export

labelGenesInModules <- function(wgcna.net, gene.annotation, residuals.matrix){

    module.colours <- labels2colors(wgcna.net$colors)
    module.labels <- data.frame(wgcna.net$colors)
    rownames(module.labels) <- rownames(residuals.matrix)

    rownames(gene.annotation) <- gene.annotation$gene

    gene.labels <- data.frame(Gene=rownames(module.labels),
                               Module=module.labels$wgcna.net.colors,
     			       Module.colour=module.colours,
     			       Tissue=gene.annotation[rownames(module.labels),]$status)
    
}