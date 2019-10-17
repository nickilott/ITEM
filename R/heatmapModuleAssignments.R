#' Plot a heatmap of residuals matrix with module assignment
#'
#' Heatmap module assignments using residualised matrix
#' @param residualsMatrix residuals of regressed tissue matrix
#' @param moduleLabels module labels for genes
#' @param metadata metadata to take sample tissue annotations from 
#' @param tissue.var variable in metadata corresponding to tissue annotation
#' @import pheatmap
#' @examples
#' heatmapTissueDE(residualsMatrix, moduleLabels, metadata, tissue.var="Tissue.location")
#' @export

heatmapModuleAssignments <- function(residualsMatrix, moduleLabels, metadata, tissue.var="Tissue.location", scale="row"){

    # remove genes in module 0
    moduleLabels <- moduleLabels[moduleLabels$Module != 0,]
    genes <- moduleLabels$Gene

    # subset input
    residualsMatrix.sub <- residualsMatrix[genes,]

    # define heatmap colours
    colours <- colorRampPalette(c("blue", "white", "red"))(75) 

    # define column annotations
    Tissue <- metadata[,tissue.var]
    sample.annotation <- data.frame("Tissue" = Tissue)
    rownames(sample.annotation) <- rownames(metadata)

    # define row annotations
    row.annotation <- data.frame("Module" = moduleLabels$Module.colour)
    rownames(row.annotation) <- genes

    # draw the heatmap
    pheatmap(residualsMatrix.sub,
             color=colours,
	     scale=scale,
	     annotation_col=sample.annotation,
	     annotation_row=row.annotation,
	     show_rownames=FALSE)
}
