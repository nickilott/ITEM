#' Run WCGNA
#'
#' This function regresses a normalised matrix for the effect of a covariate
#' and uses this as input into WGCNA.
#' @param normalised.counts normalised (cpm) counts
#' @param metadata metadata with rows in order of columns of counts matrix
#' @param regress whether to regress out a variable that contributes to variation
#' @param regress.on variable in metadata to regress on
#' @return WGCNA results
#' @import WGCNA
#' @examples
#' runWGCNA(normalised.counts, metadata, regress=TRUE, regress.on="Tissue.location")
#' @export

runWGCNA <- function(normalised.counts,
                     metadata,
                     regress=TRUE,
                     regress.on="Tissue.location",
		     wgcna.power=6,
		     wgcna.maxBlockSize=30000,
		     wgcna.minModuleSize=30,
		     wgcna.reassignThreshold=0,
		     wgcna.mergeCutHeight=0.25,
		     wgcna.numericLabels=TRUE,
		     wgcna.pamRespectsDendro=FALSE,
		     wgcna.saveTOMS=TRUE,
		     wgcna.saveTOMFileBase="test",
		     wgcna.verbose=3,
                     wgcna.TOMType="signed",
                     wgcna.type="signed"){

    if (regress == TRUE & regress.on %in% colnames(metadata)){
        input.counts <- regressOnTissue(normalised.counts, metadata, tissue.var=regress.on)}
    else{
        stop("regress.on variable name must match a column in metadata")}

    if (regress == FALSE){
        input.counts <- normalised.counts}

    net = blockwiseModules(t(input.counts),
                           maxBlockSize=wgcna.maxBlockSize,
                           power=wgcna.power,
			   TOMType=wgcna.TOMType,
			   minModuleSize=wgcna.minModuleSize,
                           reassignThreshold=wgcna.reassignThreshold,
			   mergeCutHeight=wgcna.mergeCutHeight,
			   numericLabels=wgcna.numericLabels,
                           pamRespectsDendro=wgcna.pamRespectsDendro,
                           saveTOMs=TRUE,
                           saveTOMFileBase="ITEMTest",
                           verbose=wgcna.verbose,
                           type=wgcna.type)

    return(list(input.counts, net))
}
    
