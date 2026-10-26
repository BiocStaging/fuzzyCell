#' @keywords internal
"_PACKAGE"

#' @importFrom stats kmeans lm p.adjust predict as.formula reorder setNames
#' @importFrom stats var
#' @importFrom utils head
#' @importFrom methods is
#' @importFrom rlang .data :=
#' @importFrom SummarizedExperiment colData assay
#' @import ggplot2
NULL

# column names used via non-standard evaluation (dplyr / data.table) - silences
# R CMD check
utils::globalVariables(c(
    "xRange", "yRange", "n", "cellTypeI", "cellTypeJ", "wj", "edge", "d", "N", "i", "j", "s",
    "from", "to", "leaf", "group", "child", "Var1", "Var2", "Freq", "cellType", "broad",
    "imageID", "subject", "proportion", "niche", "lineage", "estimate", "statistic", "p.value",
    "colocalisation", "value", "weightedValue", "weight",
    "misclassType", "x", "y", "normAScaled", "normBScaled", "normScaled", "type", "parent", ".", "gene", "plotValue", "label", "proportionUnweighted", "enrichment", "proportionWeighted", "..relevantCols"
))
