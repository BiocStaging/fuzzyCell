# The parent population of each pair's `to` label for relative co-localisation:
# its own group, or an override; pairs whose `to` has no parent other than
# itself are dropped. Returns a parentDf for spicyR::kontextualEngine().
relativeTriples = function(pairs, cellTypeTree, overrides) {
    cellMap = buildCellTypeMap(cellTypeTree)$cellMap
    leafToGroup = stats::setNames(cellMap$group, cellMap$cellType)
    parent = ifelse(pairs$to %in% names(overrides), unname(overrides[pairs$to]), unname(leafToGroup[pairs$to]))
    keep = !is.na(parent) & parent != pairs$to & pairs$from != pairs$to
    out = data.frame(from = pairs$from[keep], to = pairs$to[keep], parent_name = parent[keep],
        stringsAsFactors = FALSE)
    out$parent = lapply(out$parent_name, function(p) resolveTreeLeaves(cellTypeTree, p))
    out
}

# Per-image effects of a spicyR result in long format (imageID, from, to,
# colocalisation, and parent for relative co-localisation).
perImageEffects = function(res) {
    comp = res$comparisons
    out = do.call(rbind, lapply(seq_len(nrow(comp)), function(k) {
        d = data.frame(imageID = res$imageIDs, from = comp$from[k], to = comp$to[k],
            colocalisation = res$pairwiseAssoc[[comp$labels[k]]], stringsAsFactors = FALSE)
        if (!is.null(comp$parent)) d$parent = comp$parent[k]
        d
    }))
    out[!is.na(out$colocalisation), , drop = FALSE]
}

# The model table of a spicyR result for the requested pairs: estimate (the
# difference between the conditions, second level minus the first), its
# standard error, p-value and BH-adjusted p-value over the pairs kept.
modelTable = function(res, pairs = NULL) {
    tab = res$cellResults
    if (!is.null(pairs)) tab = tab[paste(tab$from, tab$to) %in% paste(pairs$from, pairs$to), , drop = FALSE]
    out = data.frame(from = tab$from, to = tab$to, stringsAsFactors = FALSE)
    if (!is.null(tab$parent)) out$parent = tab$parent
    out$estimate = tab$excess_difference
    out$se = tab$se
    out$statistic = tab$excess_difference / tab$se
    out$p.value = tab$p_value
    if (!is.null(tab$side)) out$side = tab$side
    out$p.adj = stats::p.adjust(out$p.value, "BH")
    rownames(out) = NULL
    out
}

#' Co-localisation analysis with soft cell type weights
#'
#' Tests, for every pair of leaf cell types and groups, whether the
#' co-localisation of the two differs between the two conditions of an outcome,
#' once with the original labels and once with the bootstrap confidences as soft
#' cell type weights. The tests are those of spicyR (`spicyR::spicy()`), and the
#' relative co-localisation is Statial's Kontextual test
#' (`spicyR::kontextualEngine()`).
#'
#' @details
#' Co-localisation of a pair `from` -> `to` is, by default, the extra fraction of
#' `to` cells with a `from` cell within radius `r`, beyond the `to` labels being a
#' random choice among the cells that are not `from` cells (`effect =
#' "allocation"`); `effect = "count"` gives the extra number of `from` cells within
#' `r` of each `to` cell. Images are combined within subjects and the conditions
#' compared by spicyR's frailty model with a cluster-robust variance, with
#' subjects as the units.
#'
#' With weights (see [confidenceWeights()]), a cell of a bootstrapped group
#' belongs to each subtype of the group with its bootstrap confidence. A `to` cell
#' counts with its weight, and its chance of having a `from` neighbour is computed
#' as if the labels of its neighbours were drawn from their confidences (for the
#' count, it is the summed `from` confidence of its neighbours). Under the null
#' hypothesis the `to` weights are shuffled over the cells, which gives the
#' expected value and variance of each image's statistic. With the original labels
#' (weights of 0 and 1) this is spicyR's own test.
#'
#' Relative co-localisation compares the co-localisation of `from` with `to` to
#' that of `from` with a parent population of `to`, so that attraction that only
#' arises because both cell types prefer the same region is not counted. By
#' default the parent of a `to` cell type is its own group; `parentOverrides` and
#' `autoOverrideParents` give broader parents. Pairs whose `to` has no parent
#' other than itself are not tested.
#'
#' A label is not paired with itself, nor with a group that contains it.
#'
#' @param spe A `SummarizedExperiment` (e.g. `SpatialExperiment`) with a
#'   cell ID, cell type, image and coordinates for every cell (see
#'   [getCellData()]).
#' @param cellTypeTree A cell type tree, see [cellTypeTree].
#' @param confList Confidence list from [bootstrapLineageConfidence()].
#' @param outcomeCol Name of the `colData(spe)` column with the per-image
#'   outcome. It must have exactly two levels and a single value per image.
#' @param subject Name of the `colData(spe)` column identifying the subject (for
#'   example the patient) each image comes from. Images of one subject are
#'   combined, and subjects are the units of the test. With `NULL` every image is
#'   its own subject.
#' @param r Radius within which cells count as neighbours, in the units of the
#'   coordinates; a single positive number.
#' @param effect `"allocation"` or `"count"`, see Details.
#' @param cores Number of threads.
#' @param parentOverrides Named character vector mapping a `to` label to the
#'   ancestor node (a group, wrapper or lineage in `cellTypeTree`) that
#'   relative co-localisation should compare against instead of the label's
#'   own group, for example `c("Fibroblast lineage" = "Mesenchymal")`.
#' @param autoOverrideParents Whether to add overrides found by
#'   [autoParentOverrides()]. Use `FALSE` to condition only on the parents you
#'   list in `parentOverrides`.
#' @param cellType,cellID,imageID Names of the `colData(spe)` columns holding
#'   each cell's cell type label, a unique cell ID (matching the `cellID`
#'   column of `confList`), and the image it is in.
#' @param x,y Names of the columns holding the cell coordinates. For a
#'   `SpatialExperiment` whose coordinates are in `spatialCoords()`, the
#'   defaults `"x"` and `"y"` refer to its first two coordinates.
#'
#' @return A list with the tests `modelWeighted`, `modelUnweighted`,
#'   `modelRelativeWeighted` and `modelRelativeUnweighted` (one row per pair:
#'   `estimate`, the difference for the second level of the outcome compared with
#'   the first, its `se`, `statistic`, `p.value` and BH-adjusted `p.adj`, and for
#'   co-localisation the `side`, `"attract"` or `"avoid"`), the per-image effects
#'   `weighted`, `unweighted`, `relativeWeighted` and `relativeUnweighted`
#'   (columns `imageID`, `from`, `to`, `colocalisation`, and `parent` for the
#'   relative version), and the spicyR results themselves (`spicy`, a list of
#'   `weighted`, `unweighted`, `relativeWeighted` and `relativeUnweighted`), for
#'   spicyR's plotting functions.
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 10,
#'     seed = 1)
#' res = runColocalisationModels(spe, heartCellTypeTree, confList,
#'     outcomeCol = "resolution_broad", subject = "subject")
#' head(res$modelWeighted)
#' }
#' @export
runColocalisationModels = function(spe, cellTypeTree, confList, outcomeCol, r = 50, subject = NULL,
                                   effect = c("allocation", "count"), cores = 1,
                                   parentOverrides = character(0), autoOverrideParents = TRUE,
                                   cellType = "cellType", cellID = "cellID",
                                   imageID = "imageID", x = "x", y = "y") {
    # check the arguments before doing any work
    checkSpe(spe)
    effect = match.arg(effect)
    checkRadii(r)
    if (length(r) != 1) {
        stop("`r` must be a single positive number (the radius)", call. = FALSE)
    }
    if (!is.logical(autoOverrideParents) || length(autoOverrideParents) != 1 ||
        is.na(autoOverrideParents)) {
        stop("`autoOverrideParents` must be TRUE or FALSE", call. = FALSE)
    }
    checkColumn(spe, outcomeCol, "outcomeCol")
    if (!is.null(subject)) checkColumn(spe, subject, "subject")

    cells = standardiseCells(getCellData(spe), list(cellID = cellID, cellType = cellType,
        imageID = imageID, x = x, y = y))
    checkOutcome(cells, outcomeCol, subject, oneImagePerSubject = FALSE)
    messageIfNoSubject(subject)
    if (!is.null(subject)) cells$.subject = cells[[subject]]

    message("Generating comparisons...")
    comparisons = generateComparisons(cellTypeTree, confList)
    pairs = dplyr::bind_rows(comparisons$groupGroup, comparisons$groupLeaf, comparisons$leafLeaf) |>
        dplyr::distinct(from, to)
    # a label paired with itself, or with a group that contains it, is not a pair of distinct cell populations
    covers = function(a, b) vapply(seq_along(a), function(k) b[k] %in% resolveTreeLeaves(cellTypeTree, a[k]), TRUE)
    pairs = pairs[pairs$from != pairs$to & !covers(pairs$from, pairs$to) & !covers(pairs$to, pairs$from), ,
        drop = FALSE]
    effectiveOverrides = if (autoOverrideParents) {
        auto = autoParentOverrides(cellTypeTree)
        c(auto[!names(auto) %in% names(parentOverrides)], parentOverrides)
    } else {
        parentOverrides
    }
    triples = relativeTriples(pairs, cellTypeTree, effectiveOverrides)

    weights = list(weighted = weightMatrix(cells, cellTypeTree, confList, TRUE),
        unweighted = weightMatrix(cells, cellTypeTree, confList, FALSE))
    sub = if (is.null(subject)) NULL else ".subject"
    labels = unique(c(pairs$from, pairs$to))

    res = list()
    for (v in names(weights)) {
        message("Testing ", v, " co-localisation...")
        res[[v]] = spicyR::spicy(cells, condition = outcomeCol, subject = sub, r = r, from = labels, to = labels,
            effect = effect, cores = cores, cellTypeWeights = weights[[v]])
        rel = paste0("relative", toupper(substring(v, 1, 1)), substring(v, 2))
        message("Testing ", v, " relative co-localisation...")
        res[[rel]] = if (nrow(triples)) spicyR::kontextualEngine(cells, triples, condition = outcomeCol,
            subject = sub, r = r, cellTypeWeights = weights[[v]])
    }

    model = function(x, p = NULL) if (is.null(x)) data.frame() else modelTable(x, p)
    perImage = function(x) if (is.null(x)) data.frame() else perImageEffects(x)
    list(
        weighted = perImage(res$weighted), unweighted = perImage(res$unweighted),
        modelWeighted = model(res$weighted, pairs), modelUnweighted = model(res$unweighted, pairs),
        relativeWeighted = perImage(res$relativeWeighted), relativeUnweighted = perImage(res$relativeUnweighted),
        modelRelativeWeighted = model(res$relativeWeighted),
        modelRelativeUnweighted = model(res$relativeUnweighted),
        spicy = res
    )
}
