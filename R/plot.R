
plot <- function(){
    UseMethod("plot")
    }


#' Plot
#'
#' Scatterplot the relationship between overlaps (enrichments) and eigengene correlations
#' between co-expression modules present across different tissues
#' @param ITEMOverlapSet
#' @return ggplot object
#' @import ggplot2
#' @examples
#' plot(ITEMOverlapSet)
#' @export

plot.ITEMOverlapSet <- function(ITEMOverlapSet, threshold.f=1, threshold.cor=0.6){

    overall <- ITEMOverlapSet$overall
    overall$foldEnrichment[overall$foldEnrichment == -Inf] <- 0
    overall$foldEnrichment[overall$foldEnrichment == Inf] <- max(overall$foldEnrichment[is.finite(overall$foldEnrichment)])

    overall$status <- ifelse(overall$cor.p < 0.05 & overall$overlap.p < 0.05, "Significant module overlap and cor", "NS")
    overall$status <- ifelse(overall$cor.p < 0.05 & overall$overlap.p > 0.05, "Significant module-module cor", overall$status)
    overall$status <- ifelse(overall$cor.p > 0.05 & overall$overlap.p < 0.05, "Significant module overlap", overall$status)

    p1 <- ggplot(overall, aes(x=foldEnrichment, y=cor, size=overlap, color=status))
    p1 <- p1 + geom_point() + geom_vline(xintercept=threshold.f, linetype="dashed", color="blue")
    p1 <- p1 + scale_color_manual(values=c("grey", "red3", "blue3", "orange")) + theme_bw()
    p1 <- p1 + xlab("Overlap fold enrichment") + ylab("Module-module correlation")
    p1 <- p1 + ggtitle("Inter-tissue module features") + xlim(0, max(overall$foldEnrichment) + 1)
    return(p1)
}


