#' Create ITEMDataSet class
#'
#' Create ITEMDataSet class
#' @param results.list list of results from ITEM()
#' @return ITEMDataSet object
#' @examples
#' ITEMDataSet(results.list)
#' @export

ITEMDataSet <- function(results.list){

    ref <- results.list[[1]]
    comparison <- results.list[[2]]
    net.name1 <- paste0(ref, ".net")
    net.name2 <- paste0(comparison, ".net")

    resids.matrix1 <- paste0(ref, ".residualsMatrix")
    resids.matrix2 <- paste0(comparison, ".residualsMatrix")

    labels1 <- paste0(ref, ".moduleLabels")
    labels2 <- paste0(comparison, ".moduleLabels")

    metadata1 <- paste0(ref, ".metadata")
    metadata2 <- paste0(comparison, ".metadata")

    names(results.list) <- c("ref",
                             "comparison",
                             "filteredCounts",
			     "DESeq2Results",
                             "diffGenes",
			     resids.matrix1,
			     net.name1,
			     resids.matrix2,
                             net.name2,
			     labels1,
			     labels2,
			     metadata1,
			     metadata2)
    class(results.list) <- "ITEMDataSet"
    return(results.list)
    }

####
####
summary <- function(obj){

    UseMethod("summary")
    }

#' Create ITEMDataSet summary method
#'
#' Create ITEMDataSet summary method
#' @param obj ITEMDataSet object
#' @return printed summary
#' @examples
#' summary(ITEMDataSet)
#' @export

summary.ITEMDataSet <- function(obj){

    total.genes <- nrow(obj$filteredCounts)
    ngenes.diff <- nrow(obj$DESeq2Results[obj$DESeq2Results$padj < 0.05 & !(is.na(obj$DESeq2Results$padj)),])

    nmodules1 <- obj[[10]]
    nmodules1 <- length(unique(nmodules1$Module[nmodules1$Module != 0]))

    nmodules2 <- obj[[11]]
    nmodules2 <- length(unique(nmodules2$Module[nmodules2$Module != 0]))


    net1 <- obj[[7]]
    net2 <- obj[[9]]

    cat("\n## DESeq2 results\n\n")
    cat(paste0("Total number of genes analysed = ", total.genes, "\n"))
    cat(paste0("Number of genes differentially expressed by tissue = ", ngenes.diff, "\n"))
    cat("\n\n")
    cat("## WGCNA results\n\n")
    cat("Number of modules for tissue", obj$ref, " = ", nmodules1, "\n")
    cat("Number of genes in modules\n")
    print(table(net1$colors))

    cat("\n\nNumber of modules for tissue", obj$comparison, " = ", nmodules2, "\n")
    cat("Number of genes in modules\n")
    print(table(net2$colors))

    }