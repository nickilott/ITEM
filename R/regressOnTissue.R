#' Regress out effect of tissue
#' Uses linear model to regress out the effect of tissue on gene expression
#' @param row row of normalised counts matrix to perform the regressions
#' @param tissue variable containing a vector of tissues that correspond to row
#' @examples
#' residualise(row, tissues)
#' @export

residualise <- function(row, tissues, scale=TRUE){

	    # save standardised residuals after regressing
	    # out effect of tissue

	    linmod <- lm(unlist(row)~as.factor(tissues))
	    residuals <- residuals(linmod)
	    if (scale==TRUE){
	       residuals <- scale(residuals)}
	    else{
		residuals <- residuals
	    }
	    return(residuals)
	    }

#' Regress out effect of tissue
#' Uses linear model to regress out the effect of tissue on gene expression
#' @param normalised.counts normalised counts matrix (log2(cpm))
#' @param metadata data to take tissue variable
#' @param tissue.var variable in metadata that corresponds to tissue
#' @param scale boolean of whether to scale residuals after regression
#' @examples
#' regressOnTissue(normalised.counts, metadata, tissue.var="Tissue.location", scale=TRUE)
#' @export

regressOnTissue <- function(normalised.counts, metadata, tissue.var="Tissue.location", scale=TRUE){

    tissues <- metadata[,tissue.var]
    regressed.data <- data.frame(t(apply(normalised.counts, 1, function(x) residualise(x, tissues, scale=scale))))
    colnames(regressed.data) <- colnames(normalised.counts)
    return(regressed.data)
    }


