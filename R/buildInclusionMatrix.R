#' Build inclusion matrix
#'
#' Build a matrix that includes modules as columns and genes as rows annotated
#' with a presence/absence score of 1/0. This can then be used for clustering
#' purposes
#' @param ITEMOverlapSet
#' @return data.frame
#' @import reshape2
#' @examples
#' buildInclusionMatrix(ITEMDataSet)
#' @export

buildInclusionMatrix <- function(ITEMDataSet){

    name1 <- ITEMDataSet$ref
    name2 <- ITEMDataSet$comparison

    # Get module assignments and rename
    modules1 <- ITEMDataSet[[10]]
    modules2 <- ITEMDataSet[[11]]

    modules1$Module <- paste0(name1, ".", modules1$Module)
    modules2$Module <- paste0(name2, ".", modules2$Module)

    modules1 <- modules1[,1:2]
    modules2 <- modules2[,1:2]

    # cast the data into matrix
    modules1.cast <- dcast(modules1, Gene~Module)
    modules2.cast <- dcast(modules2, Gene~Module)

    # Merge the data
    inclusion.matrix <- data.frame(cbind(modules1.cast, modules2.cast[2:ncol(modules2.cast)]))
    rownames(inclusion.matrix) <- inclusion.matrix[,1]
    inclusion.matrix <- inclusion.matrix[,2:ncol(inclusion.matrix)]
    inclusion.matrix[!(is.na(inclusion.matrix))] <- 1
    inclusion.matrix[is.na(inclusion.matrix)] <- 0

    return(inclusion.matrix)

}
