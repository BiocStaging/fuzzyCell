#' Cell type weights from bootstrap confidences
#'
#' Builds the soft cell type membership of every cell, the weights that the
#' weighted analyses pass to spicyR and lisaClust: one column per leaf cell type
#' and per group of the tree. A cell of a bootstrapped group has, for each
#' subtype of the group, its bootstrap confidence in that subtype; every other
#' cell has weight 1 for its own cell type. A group's column is 1 for each of its
#' cells, since the bootstrap reclassifies cells only within their group.
#'
#' @inheritParams runColocalisationModels
#' @param weighted `TRUE` for the bootstrap confidences, `FALSE` for the
#'   original labels (0/1 weights) on the same columns.
#' @param cellType,cellID Names of the `colData(spe)` columns holding each
#'   cell's cell type label and a unique cell ID (matching the `cellID` column
#'   of `confList`).
#'
#' @return A numeric matrix with one row per cell of `spe` (in its order) and one
#'   named column per leaf cell type and group.
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' w = confidenceWeights(lungSubset, lungCellTypeTree, lungConfList)
#' head(w[, c("Arteriole", "Venous", "Endothelial lineage")])
#' @export
confidenceWeights = function(spe, cellTypeTree, confList, weighted = TRUE,
                             cellType = "cellType", cellID = "cellID") {
    cells = standardiseCells(getCellData(spe), list(cellID = cellID, cellType = cellType))
    weightMatrix(cells, cellTypeTree, confList, weighted)
}

# the weight matrix of a standardised per-cell table (cellID, cellType)
weightMatrix = function(cells, cellTypeTree, confList, weighted) {
    treeInfo = buildCellTypeMap(cellTypeTree)
    cellMap = treeInfo$cellMap
    leaves = cellMap$cellType
    groups = treeInfo$groupsWithChildren
    missingTypes = setdiff(unique(as.character(cells$cellType)), leaves)
    if (length(missingTypes) > 0) {
        stop("Cell types not in cellTypeTree: ", paste(utils::head(missingTypes, 10), collapse = ", "),
            call. = FALSE)
    }
    labels = c(leaves, setdiff(groups, leaves))
    w = matrix(0, nrow(cells), length(labels), dimnames = list(NULL, labels))
    ct = as.character(cells$cellType)
    w[cbind(seq_len(nrow(cells)), match(ct, labels))] = 1
    groupOf = cellMap$group[match(ct, cellMap$cellType)]
    for (g in groups) w[which(groupOf == g), g] = 1
    if (!weighted) return(w)
    for (g in intersect(names(confList), groups)) {
        conf = confList[[g]]
        idx = which(groupOf == g)
        rows = match(as.character(cells$cellID[idx]), as.character(conf$cellID))
        if (anyNA(rows)) {
            stop("confList[[\"", g, "\"]] has no confidence for ", sum(is.na(rows)), " of the group's cells ",
                "(matched by cell ID)", call. = FALSE)
        }
        subtypes = intersect(setdiff(names(conf), "cellID"), leaves)
        w[idx, cellMap$cellType[cellMap$group %in% g]] = 0
        for (s in subtypes) w[idx, s] = conf[[s]][rows]
    }
    w
}
