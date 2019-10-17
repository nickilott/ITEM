#' Run PCA
#'
#' Run PCA on matrix
#' @param df data frame of normalised counts
#' @export
#' @examples
#' runPCA(df)

runPCA <- function(df, scale=TRUE){

    pc <- prcomp(t(df), scale=scale)
    return (pc)
}

##################################################
##################################################
##################################################

#' Get variance explained
#'
#' Get variance explained for prcomp object
#' @param pc prcomp object
#' @param component string (component to get VE for)
#' @export
#' @examples
#' getVE(runPCA(df))

getVE <- function(pc, component="PC1"){

    pve <- summary(pc)$importance[,component][[2]]
    return (pve)
}
