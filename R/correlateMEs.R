#' Correlate module eigengenes between tissue modules
#'
#' Correlate the module eigengene values for matched samples across tissues
#' i.e. modules called in each module separately
#' @param ITEMDataSet
#' @return matrix correlation matrix
#' @examples
#' correlateMEs(ITEMDataSet)
#' @export

correlateMEs <- function(ITEMDataSet, individual.var="Patient.ID"){

    mes1 <- ITEMDataSet[[7]]
    mes2 <- ITEMDataSet[[9]]
    mes1 <- mes1$MEs
    mes2 <- mes2$MEs

    # rownames come from metadata
    rownames(mes1) <- ITEMDataSet[[12]][,individual.var]
    rownames(mes2) <- ITEMDataSet[[13]][,individual.var]

    # Rename ME column names so that the tissues
    # are distinguishable
    colnames(mes1) <- paste0(ITEMDataSet$ref, ".", gsub("ME", "", colnames(mes1)))
    colnames(mes2) <- paste0(ITEMDataSet$comparison, ".", gsub("ME", "", colnames(mes2)))

    # get matched samples
    matched <- intersect(rownames(mes1), rownames(mes2))
    mes1 <- mes1[matched,]
    mes2 <- mes2[matched,]

    all.mes <- data.frame(cbind(mes1, mes2))
    rownames(all.mes) <- rownames(mes1)

    # do the correlation - still working out whether to perform the adjustment
    cor.mes <- corr.test(all.mes, adjust="BH")

    return(cor.mes)

}