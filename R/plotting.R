# Cell type labels, confidence scores and markers ----
# Spatial maps of one image: the original labels, the bootstrap confidence
# scores, and marker expression within a group. relabelLineage(),
# getSubtypeWeight() and getTwoLineageWeightDF() are internal helpers for the
# confidence maps.

# relabel everything in a lineage to a specific subtype
relabelLineage = function(speSub, group, newLabel) {
    if (is.null(speSub$misclassType)) speSub$misclassType = speSub$cellType
    speSub$misclassType[speSub$group == group] = newLabel
    speSub
}

getSubtypeWeight = function(speSub, confList, group, subtype) {
    n = if (is.data.frame(speSub)) nrow(speSub) else ncol(speSub)
    weight = rep(1, n)
    # subtype == group means "the whole lineage", not one bootstrapped subtype
    # within it -
    # group membership is a hard, deterministic classification (no confList
    # column for it,
    # since confList's columns are the group's subtypes), so those cells just
    # keep the
    # default hard weight of 1 instead of looking up a per-subtype confidence
    if (identical(group, subtype)) {
        return(weight)
    }
    inGroup = which(speSub$group == group)
    conf = confList[[group]]
    matchIdx = match(speSub$cellID[inGroup], conf$cellID)
    weight[inGroup] = conf[[subtype]][matchIdx]
    weight
}

# same as getSubtypeWeight(), but for two lineages being plotted side by side at
# once -
# groupA cells get their confidence of being subtypeA, groupB cells get their
# confidence
# of being subtypeB, everyone else gets the default hard weight of 1
getTwoLineageWeightDF = function(speSub, confList, groupA, subtypeA, groupB, subtypeB) {
    weightA = getSubtypeWeight(speSub, confList, groupA, subtypeA)
    weightB = getSubtypeWeight(speSub, confList, groupB, subtypeB)
    ifelse(speSub$group == groupA, weightA, ifelse(speSub$group == groupB, weightB, 1))
}

#' Spatial plots of cell types and confidence scores
#'
#' Plots the cells of one image, highlighting a cell type or group, or showing
#' bootstrap confidence scores.
#'
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with a cell
#'   ID, cell type, image and coordinates for every cell in the form given by
#'   [getCellData()] (`cellID`, `cellType`, `imageID`, `x` and `y`).
#' @param cellTypeTree A cell type tree, see [cellTypeTree]. Optional for
#'   `plotHardSpatialPair()`, where it is only needed to use a group name as a
#'   label.
#' @param confList Confidence list from [bootstrapLineageConfidence()].
#' @param imageID The image to plot. An error is raised if it is not an image in
#'   the data.
#' @param groupA,groupB,targetGroup Group names. An error is raised for a group
#'   that has no cells in the data.
#' @param subtypeA,subtypeB,targetSubtype Cell types within the groups, which
#'   must have confidence scores in `confList` (for `plotHardSpatialPair()`, a
#'   cell type or a group name in the data).
#' @param encoding How to show the confidence scores: `"alpha"` (the default),
#'   as transparency, or `"colour"`, as a colour gradient. Either way each of the
#'   two cell types has its own legend.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 10,
#'     seed = 1)
#' img = unique(spe$imageID)[1]
#' plotHardSpatialPair(spe, img, "CD8+ T cells", "Macrophages",
#'     heartCellTypeTree)
#' plotWeightedSpatial(spe, heartCellTypeTree, confList, img, "T cells",
#'     "CD8+ T cells")
#' plotOriginalLabelSpatial(spe, heartCellTypeTree, img, "T cells")
#' plotWeightedSpatialPair(spe, heartCellTypeTree, confList, img,
#'     "T cells", "CD8+ T cells", "Macrophage lineage", "Macrophages")
#' }
#' @param cellSize Size of the highlighted cells, the ones drawn in colour. The
#'   grey background cells always have size 0.7.
#' @name plotSpatial
#' @export
plotWeightedSpatialPair = function(spe, cellTypeTree, confList, imageID,
                                    groupA, subtypeA, groupB, subtypeB,
                                    encoding = c("alpha", "colour"), cellSize = 1) {
    checkCellSize(cellSize)
    checkSpe(spe)
    encoding = match.arg(encoding)

    cells = getCellData(spe)
    checkImage(cells, imageID)
    checkGroupsInData(cells, cellTypeTree, list(groupA = groupA, groupB = groupB))
    checkSubtype(confList, groupA, subtypeA, "subtypeA")
    checkSubtype(confList, groupB, subtypeB, "subtypeB")
    cellsSub = cells[cells$imageID == imageID, ]
    cellsSub = annotateCellTypeTree(withoutAnnotation(cellsSub), cellTypeTree)
    cellsSub = relabelLineage(cellsSub, groupA, subtypeA)
    cellsSub = relabelLineage(cellsSub, groupB, subtypeB)
    cellsSub$weight = getTwoLineageWeightDF(cellsSub, confList, groupA, subtypeA, groupB, subtypeB)

    p = ggplot() +
        geom_point(
            data = cellsSub |> dplyr::filter(!misclassType %in% c(subtypeA, subtypeB)),
            aes(x = x, y = y), colour = "grey80", size = 0.7, alpha = 0.5
        )

    if (encoding == "alpha") {
        # each cell type has its own alpha scale and legend (ggnewscale), keyed
        # in that cell type's colour
        p = p +
            geom_point(
                data = cellsSub |> dplyr::filter(misclassType == subtypeA),
                aes(x = x, y = y, alpha = weight), colour = "red", size = cellSize
            ) +
            scale_alpha_continuous(
                limits = c(0, 1),
                name = paste0(subtypeA, "\nconfidence")
            ) +
            ggnewscale::new_scale("alpha") +
            geom_point(
                data = cellsSub |> dplyr::filter(misclassType == subtypeB),
                aes(x = x, y = y, alpha = weight), colour = "blue", size = cellSize
            ) +
            scale_alpha_continuous(
                limits = c(0, 1),
                name = paste0(subtypeB, "\nconfidence")
            )
    } else {
        p = p +
            geom_point(
                data = cellsSub |> dplyr::filter(misclassType == subtypeA),
                aes(x = x, y = y, colour = weight), size = cellSize
            ) +
            scale_colour_gradient(
                low = "grey80", high = "red", limits = c(0, 1),
                name = paste0(subtypeA, "\nconfidence")
            ) +
            ggnewscale::new_scale_colour() +
            geom_point(
                data = cellsSub |> dplyr::filter(misclassType == subtypeB),
                aes(x = x, y = y, colour = weight), size = cellSize
            ) +
            scale_colour_gradient(
                low = "grey80", high = "blue", limits = c(0, 1),
                name = paste0(subtypeB, "\nconfidence")
            )
    }

    p +
        scale_y_reverse() +
        coord_equal() +
        labs(
            x = "", y = "",
            title = paste0(imageID, "\n", subtypeA, " (red) vs ", subtypeB, " (blue)")
        ) +
        theme_classic() +
        theme(legend.position = "right")
}

#' @rdname plotSpatial
#' @export
plotHardSpatialPair = function(spe, imageID, subtypeA, subtypeB, cellTypeTree = NULL, cellSize = 1) {
    checkCellSize(cellSize)
    checkSpe(spe)
    cells = getCellData(spe)
    checkImage(cells, imageID)
    checkLabelsInData(cells, cellTypeTree, list(subtypeA = subtypeA, subtypeB = subtypeB))

    # with a cell type tree the cells get a `group` column, so that a group name
    # (not only a cell type) can be used as a label; without one, a group name
    # matches nothing and its layer is empty
    if (!is.null(cellTypeTree)) cells = annotateCellTypeTree(withoutAnnotation(cells), cellTypeTree)
    cellsSub = cells[cells$imageID == imageID, ]

    # a label may be a specific leaf cellType, or an entire lineage/group -
    # matching by
    # group covers "plot this whole lineage as one colour" (e.g. all Fibroblast
    # lineage
    # cells), matching by cellType covers a single subtype as before
    matchLabel = function(label) {
        m = if (!is.null(cellsSub$group) && label %in% cellsSub$group) {
            cellsSub$group == label
        } else {
            cellsSub$cellType == label
        }
        m[is.na(m)] = FALSE # cells with an NA group shouldn't leak into any layer as NA rows
        m
    }
    matchA = matchLabel(subtypeA)
    matchB = matchLabel(subtypeB)

    ggplot() +
        geom_point(
            data = cellsSub[!matchA & !matchB, ],
            aes(x = x, y = y), colour = "grey80", size = 0.7, alpha = 0.5
        ) +
        geom_point(
            data = cellsSub[matchA, ],
            aes(x = x, y = y), colour = "red", size = cellSize
        ) +
        geom_point(
            data = cellsSub[matchB, ],
            aes(x = x, y = y), colour = "blue", size = cellSize
        ) +
        scale_y_reverse() +
        coord_equal() +
        labs(
            x = "", y = "",
            title = paste0(imageID, "\n", subtypeA, " (red) vs ", subtypeB, " (blue)")
        ) +
        theme_classic() +
        theme(legend.position = "right")
}

#' @rdname plotSpatial
#' @export
plotWeightedSpatial = function(spe, cellTypeTree, confList, imageID, targetGroup, targetSubtype, cellSize = 1) {
    checkCellSize(cellSize)
    checkSpe(spe)
    cells = getCellData(spe)
    checkImage(cells, imageID)
    checkGroupsInData(cells, cellTypeTree, list(targetGroup = targetGroup))
    checkSubtype(confList, targetGroup, targetSubtype, "targetSubtype")
    cells = annotateCellTypeTree(withoutAnnotation(cells), cellTypeTree)
    imgCells = cells |> dplyr::filter(imageID == !!imageID)
    imgCells$weight = getSubtypeWeight(imgCells, confList, targetGroup, targetSubtype)

    ggplot() +
        geom_point(
            data = imgCells |> dplyr::filter(is.na(group) | group != targetGroup),
            aes(x = x, y = y), colour = "grey80", size = 0.7
        ) +
        geom_point(
            data = imgCells |> dplyr::filter(group == targetGroup),
            aes(x = x, y = y, alpha = weight), colour = "red", size = cellSize
        ) +
        scale_alpha_continuous(range = c(0.05, 1), limits = c(0, 1)) +
        coord_equal() +
        scale_y_reverse() +
        labs(
            title = paste0(imageID, " - ", targetSubtype, " weight"), alpha = "Weight",
            x = "", y = ""
        ) +
        theme_classic()
}

#' @rdname plotSpatial
#' @export
plotOriginalLabelSpatial = function(spe, cellTypeTree, imageID, targetGroup, cellSize = 1) {
    checkCellSize(cellSize)
    checkSpe(spe)
    cells = getCellData(spe)
    checkImage(cells, imageID)
    checkGroupsInData(cells, cellTypeTree, list(targetGroup = targetGroup))
    cells = annotateCellTypeTree(withoutAnnotation(cells), cellTypeTree)
    imgCells = cells |> dplyr::filter(imageID == !!imageID)

    ggplot() +
        geom_point(
            data = imgCells |> dplyr::filter(is.na(group) | group != targetGroup),
            aes(x = x, y = y), colour = "grey80", size = 0.7
        ) +
        geom_point(
            data = imgCells |> dplyr::filter(group == targetGroup),
            aes(x = x, y = y, colour = cellType), size = cellSize
        ) +
        coord_equal() +
        scale_y_reverse() +
        labs(
            title = paste0(imageID, " - ", targetGroup, " lineage (original labels)"), colour = "Cell type",
            x = "", y = ""
        ) +
        theme_classic()
}

#' Spatial plots of marker expression within a group
#'
#' Plots one image, with the cells of a group shaded by the (library-size
#' normalised, rescaled) expression of one or two markers.
#'
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with a cell
#'   ID, cell type, image and coordinates for every cell in the form given by
#'   [getCellData()], and a counts assay.
#' @param imageID The image to plot.
#' @param group,groupA,groupB Group name(s) whose cells are shaded.
#' @param marker,markerA,markerB Gene name(s), which must be rows of the assay.
#' @param colour,colourA,colourB Colour(s) used for the shading.
#' @param cellTypeTree A cell type tree, used to find the cells of each group.
#'   It is needed unless `colData(spe)` already has a `group` column.
#' @param assayName Name of the assay with the counts. Expression is normalised
#'   by each cell's library size.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' img = unique(spe$imageID)[1]
#' plotOneMarkerSpatial(spe, img, "T cells", "CD3E",
#'     cellTypeTree = heartCellTypeTree)
#' plotTwoMarkerSpatial(spe, img, "T cells", "CD3E", "purple",
#'     "Macrophage lineage", "CD68", "green", cellTypeTree = heartCellTypeTree)
#' }
#' @param cellSize Size of the highlighted cells, the ones drawn in colour. The
#'   grey background cells always have size 0.7.
#' @name plotMarkerSpatial
#' @export
plotOneMarkerSpatial = function(spe, imageID, group, marker, colour = "purple",
                                 cellTypeTree = NULL, cellSize = 1, assayName = "counts") {
    checkCellSize(cellSize)
    checkSpe(spe)
    checkAssay(spe, assayName)
    cells = getCellData(spe)
    checkImage(cells, imageID)
    if (!is.null(cellTypeTree)) cells = annotateCellTypeTree(withoutAnnotation(cells), cellTypeTree)
    if (is.null(cells$group)) {
        stop("`cellTypeTree` is needed to find the cells of the group, because ",
            "colData(spe) has no `group` column", call. = FALSE)
    }
    checkGroupsInData(cells, cellTypeTree, list(group = group))
    countMat = SummarizedExperiment::assay(spe, assayName)
    checkGene(countMat, marker, "marker", assayName)

    librarySize = Matrix::colSums(countMat) # works for dense and sparse counts
    cells$norm = countMat[marker, ] / librarySize

    imgCells = cells |> dplyr::filter(imageID == !!imageID)
    inGroup = !is.na(imgCells$group) & imgCells$group == group

    # scale to [0, 1] using its own range WITHIN the target group, so the alpha
    # below uses
    # the full transparency range instead of being washed out by whatever raw
    # library-size-normalised scale the marker happens to sit at
    rng = range(imgCells$norm[inGroup], na.rm = TRUE)
    imgCells$normScaled = if (diff(rng) == 0) {
        rep(0, nrow(imgCells))
    } else {
        pmin(pmax((imgCells$norm - rng[1]) / diff(rng), 0), 1)
    }

    ggplot() +
        geom_point(
            data = imgCells |> dplyr::filter(!inGroup),
            aes(x = x, y = y), colour = "grey80", size = 0.7
        ) +
        geom_point(
            data = imgCells |> dplyr::filter(inGroup),
            aes(x = x, y = y, alpha = normScaled), colour = colour, size = cellSize
        ) +
        scale_alpha_continuous(range = c(0.05, 1), name = paste0(marker, " (norm.)")) +
        coord_equal() +
        scale_y_reverse() +
        labs(title = paste0(imageID, " - ", marker, " in ", group), x = "", y = "") +
        theme_classic()
}

#' @rdname plotMarkerSpatial
#' @export
plotTwoMarkerSpatial = function(spe, imageID,
                                 groupA, markerA, colourA = "purple",
                                 groupB, markerB, colourB = "green",
                                 cellTypeTree = NULL, cellSize = 1, assayName = "counts") {
    checkCellSize(cellSize)
    checkSpe(spe)
    checkAssay(spe, assayName)
    cells = getCellData(spe)
    checkImage(cells, imageID)
    if (!is.null(cellTypeTree)) cells = annotateCellTypeTree(withoutAnnotation(cells), cellTypeTree)
    if (is.null(cells$group)) {
        stop("`cellTypeTree` is needed to find the cells of the groups, because ",
            "colData(spe) has no `group` column", call. = FALSE)
    }
    checkGroupsInData(cells, cellTypeTree, list(groupA = groupA, groupB = groupB))
    countMat = SummarizedExperiment::assay(spe, assayName)
    checkGene(countMat, markerA, "markerA", assayName)
    checkGene(countMat, markerB, "markerB", assayName)
    librarySize = Matrix::colSums(countMat) # works for dense and sparse counts
    cells$normA = countMat[markerA, ] / librarySize
    cells$normB = countMat[markerB, ] / librarySize

    imgCells = cells |> dplyr::filter(imageID == !!imageID)
    inA = !is.na(imgCells$group) & imgCells$group == groupA
    inB = !is.na(imgCells$group) & imgCells$group == groupB

    # scale each marker to [0, 1] using its own range WITHIN its own target
    # group, so the
    # alpha below uses the full transparency range instead of being washed out
    # by whatever
    # raw library-size-normalised scale that marker happens to sit at
    scale01 = function(x, mask) {
        rng = range(x[mask], na.rm = TRUE)
        if (diff(rng) == 0) {
            return(rep(0, length(x)))
        }
        pmin(pmax((x - rng[1]) / diff(rng), 0), 1)
    }
    imgCells$normAScaled = scale01(imgCells$normA, inA)
    imgCells$normBScaled = scale01(imgCells$normB, inB)

    ggplot() +
        geom_point(
            data = imgCells |> dplyr::filter(!inA & !inB),
            aes(x = x, y = y), colour = "grey80", size = 0.7
        ) +
        geom_point(
            data = imgCells |> dplyr::filter(inA),
            aes(x = x, y = y, alpha = normAScaled), colour = colourA, size = cellSize
        ) +
        geom_point(
            data = imgCells |> dplyr::filter(inB),
            aes(x = x, y = y, alpha = normBScaled), colour = colourB, size = cellSize
        ) +
        scale_alpha_continuous(range = c(0.05, 1), name = "Expression (norm.)") +
        coord_equal() +
        scale_y_reverse() +
        labs(
            title = paste0(
                imageID, "\n", markerA, " in ", groupA, " (", colourA, ") / ",
                markerB, " in ", groupB, " (", colourB, ")"
            ),
            x = "", y = ""
        ) +
        theme_classic()
}
# Proportions ----
# Plot of the results of runProportionsPipeline() and
# runProportionsRelativeToParentPipeline().

#' Plot one cell type's proportion by outcome
#'
#' Shows, for one cell type (or group), each image's proportion measured with
#' the original labels (unweighted) and with the bootstrap confidence weights
#' (weighted), for each level of the outcome. For the results of
#' [runProportionsPipeline()] the proportion is of all cells in the image; for
#' the results of [runProportionsRelativeToParentPipeline()] it is of the cells
#' of the parent group. Each point is one image.
#'
#' @param results Output of [runProportionsPipeline()] or
#'   [runProportionsRelativeToParentPipeline()].
#' @param cellType The cell type to show: one of the `cellType` values in
#'   `results$proportions`. For [runProportionsPipeline()] results this can also
#'   be one of the bootstrapped groups.
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with
#'   `imageID` and the outcome in `colData`.
#' @param outcomeCol Name of the outcome column in `colData(spe)`.
#'
#' @return An interactive plot, a `girafe` widget that needs the ggiraph
#'   package.
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 10,
#'     seed = 1)
#' if (requireNamespace("ggiraph", quietly = TRUE)) {
#'     res = runProportionsPipeline(spe, heartCellTypeTree, confList,
#'         outcomeCol = "resolution_broad", subject = "subject")
#'     plotProportionBoxplot(res, "CD8+ T cells", spe, "resolution_broad")
#'
#'     rel = runProportionsRelativeToParentPipeline(spe, heartCellTypeTree,
#'         confList, outcomeCol = "resolution_broad", subject = "subject")
#'     plotProportionBoxplot(rel, "CD8+ T cells", spe, "resolution_broad")
#' }
#' }
#' @export
plotProportionBoxplot = function(results, cellType, spe, outcomeCol) {
    checkSpe(spe)
    p = proportionBoxplotGG(results, cellType, getCellData(spe), outcomeCol)
    ggiraph::girafe(ggobj = p, height_svg = 3.5)
}

# the ggplot behind plotProportionBoxplot(), before it is made interactive
proportionBoxplotGG = function(results, cellType, cells, outcomeCol) {
    needed = c("imageID", "cellType", "lineage", "proportionUnweighted", "proportionWeighted")
    if (!is.list(results) || !is.data.frame(results$proportions) ||
        !all(needed %in% names(results$proportions))) {
        stop("`results` must be the output of runProportionsPipeline() or ",
            "runProportionsRelativeToParentPipeline()", call. = FALSE)
    }
    available = unique(results$proportions$cellType)
    if (!is.character(cellType) || length(cellType) != 1 || is.na(cellType)) {
        stop("`cellType` must be a single string", call. = FALSE)
    }
    if (!cellType %in% available) {
        stop("`cellType = \"", cellType, "\"` is not in the results. Available: ",
            paste(utils::head(available, 10), collapse = ", "),
            if (length(available) > 10) ", ...", call. = FALSE)
    }
    checkColumn(cells, outcomeCol, "outcomeCol")
    checkColumn(cells, "imageID", "imageID")

    # the relative-to-parent pipeline has no separate weighted/unweighted tables
    isRelative = !"weighted" %in% names(results)
    dat = results$proportions[results$proportions$cellType == cellType, ]
    yLab = if (isRelative) {
        paste0("Proportion of ", unique(dat$lineage)[1], " cells")
    } else {
        "Proportion of all cells"
    }

    outcomeDat = dplyr::distinct(cells, imageID, .data[[outcomeCol]])
    plotDat = dat |>
        dplyr::select(imageID, value = proportionUnweighted, weightedValue = proportionWeighted) |>
        dplyr::left_join(outcomeDat, by = "imageID") |>
        tidyr::pivot_longer(cols = c(value, weightedValue), names_to = "type", values_to = "plotValue") |>
        dplyr::mutate(
            type = factor(dplyr::recode(type, value = "Unweighted", weightedValue = "Weighted"),
                levels = c("Unweighted", "Weighted")),
            !!outcomeCol := as.factor(.data[[outcomeCol]])
        ) |>
        dplyr::filter(!is.na(plotValue))

    ggplot(plotDat, aes(x = type, y = plotValue, colour = .data[[outcomeCol]])) +
        geom_boxplot(aes(group = interaction(.data[[outcomeCol]], type)),
            outlier.shape = NA, width = 0.35, position = position_dodge(0.45)
        ) +
        ggiraph::geom_point_interactive(
            aes(
                tooltip = imageID, data_id = imageID,
                group = interaction(.data[[outcomeCol]], type)
            ),
            position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.45),
            alpha = 0.6, size = 1.5
        ) +
        labs(x = NULL, y = yLab, colour = NULL, title = cellType) +
        theme_classic()
}

# Co-localisation ----
# Plot of the results of runColocalisationModels().

#' Plot one pair's co-localisation by outcome
#'
#' Shows, for one pair of cell types or groups, the per-image co-localisation
#' measured with the original labels (unweighted) and with the bootstrap
#' confidence weights (weighted), for each level of the outcome. Each point is
#' one image.
#'
#' @param results Output of [runColocalisationModels()].
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with
#'   `imageID` and the outcome in `colData`.
#' @param fromType,toType The pair to show. An error is raised if the pair is
#'   not in the results; relative co-localisation (`"kontextual"`) only has
#'   pairs whose `toType` has a parent other than itself.
#' @param outcomeCol Name of the outcome column in `colData(spe)`.
#' @param metric The measure to plot: `"lfunction"` (co-localisation, the
#'   cross-type L-function) or `"kontextual"` (co-localisation relative to a
#'   parent population, the Kontextual measure).
#'
#' @return An interactive plot, a `girafe` widget that needs the ggiraph
#'   package.
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 10,
#'     seed = 1)
#' if (requireNamespace("ggiraph", quietly = TRUE)) {
#'     res = runColocalisationModels(spe, heartCellTypeTree, confList,
#'         outcomeCol = "resolution_broad", subject = "subject")
#'     plotMetricBoxplot(res, spe, "CD8+ T cells", "Macrophages",
#'         outcomeCol = "resolution_broad")
#'     plotMetricBoxplot(res, spe, "CD8+ T cells", "Macrophages",
#'         outcomeCol = "resolution_broad", metric = "kontextual")
#' }
#' }
#' @export
plotMetricBoxplot = function(results, spe, fromType, toType, outcomeCol,
                              metric = c("lfunction", "kontextual")) {
    checkSpe(spe)
    metric = match.arg(metric)
    p = metricBoxplotGG(results, getCellData(spe), fromType, toType, outcomeCol, metric)
    ggiraph::girafe(ggobj = p, height_svg = 3.5)
}

# the ggplot behind plotMetricBoxplot(), before it is made interactive
metricBoxplotGG = function(results, cells, fromType, toType, outcomeCol, metric) {
    source = switch(metric,
        lfunction = list(
            unweighted = results$unweighted, weighted = results$weighted,
            valueCol = "colocalisation", yLab = "Co-localisation"
        ),
        kontextual = list(
            unweighted = results$relativeUnweighted, weighted = results$relativeWeighted,
            valueCol = "colocalisation", yLab = "Relative co-localisation"
        )
    )

    needed = c("imageID", "from", "to", "colocalisation")
    if (!is.data.frame(source$unweighted) || !is.data.frame(source$weighted) ||
        !all(needed %in% names(source$unweighted))) {
        stop("`results` must be the output of runColocalisationModels()", call. = FALSE)
    }
    for (arg in c("fromType", "toType")) {
        value = get(arg)
        if (!is.character(value) || length(value) != 1 || is.na(value)) {
            stop("`", arg, "` must be a single cell type or group name", call. = FALSE)
        }
    }
    if (!any(source$unweighted$from == fromType & source$unweighted$to == toType)) {
        stop("The pair '", fromType, "' -> '", toType, "' is not in the results",
            if (metric == "kontextual") paste0(" for metric = \"kontextual\". Relative ",
                "co-localisation only has pairs whose `toType` has a parent other ",
                "than itself (see ?runColocalisationModels)"), call. = FALSE)
    }
    checkColumn(cells, outcomeCol, "outcomeCol")

    outcomeDat = dplyr::distinct(cells, imageID, .data[[outcomeCol]])

    plotDat = dplyr::inner_join(
        source$unweighted |> dplyr::select(imageID, from, to, value = dplyr::all_of(source$valueCol)),
        source$weighted |> dplyr::select(imageID, from, to, weightedValue = dplyr::all_of(source$valueCol)),
        by = c("imageID", "from", "to")
    ) |>
        dplyr::filter(from == fromType, to == toType) |>
        dplyr::left_join(outcomeDat, by = "imageID") |>
        tidyr::pivot_longer(cols = c(value, weightedValue), names_to = "type", values_to = "plotValue") |>
        dplyr::mutate(type = dplyr::recode(type, value = "Unweighted", weightedValue = "Weighted"))

    p = ggplot(plotDat, aes(x = type, y = plotValue, colour = .data[[outcomeCol]])) +
        geom_boxplot(aes(group = interaction(.data[[outcomeCol]], type)),
            outlier.shape = NA, width = 0.35, position = position_dodge(0.45)
        ) +
        ggiraph::geom_point_interactive(
            aes(
                tooltip = imageID, data_id = imageID,
                group = interaction(.data[[outcomeCol]], type)
            ),
            position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.45),
            alpha = 0.6, size = 1.5
        ) +
        labs(
            x = NULL, y = source$yLab, colour = NULL,
            title = paste0(fromType, " \u2192 ", toType)
        ) +
        theme_classic()

    p
}

# Niches ----
# Plots of the niches from runNicheIdentification() and testNicheOutcome().

#' Niche plots
#'
#' `plotNicheSpatial()` shows the cells of one niche in one image.
#' `nicheCompositionMap()` shows the relative enrichment of each cell type in
#' each niche as a bubble plot or heatmap. `nicheBarplot()` shows it for a single
#' niche. `nicheProportionBoxplot()` compares a niche's per-image proportion
#' between outcome groups.
#'
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`). For the
#'   niche plots the niche assignments need to be in `colData`, for example
#'   `spe$nicheWeighted = niches$nicheWeighted[spe$cellID]`; the other plots
#'   need the cell type, image and coordinates in the form given by
#'   [getCellData()].
#' @param imageID The image to plot.
#' @param niche For `plotNicheSpatial()` and the proportion boxplot, the niche
#'   to show. For `nicheCompositionMap()` and `nicheBarplot()`, the name of the
#'   `colData` column holding niche labels.
#' @param nicheCol Column of `colData(spe)` with niche assignments.
#' @param cellType Column of `colData(spe)` with cell type labels.
#' @param type `"bubble"` or `"heatmap"`.
#' @param limit Lower and upper limit of the relative enrichment.
#' @param ... Further arguments passed to [pheatmap::pheatmap()].
#' @param target The niche to plot in `nicheBarplot()`.
#' @param proportions The `propWeighted` or `propUnweighted` element of the
#'   results of [testNicheOutcome()]: the table you pass decides whether the
#'   weighted or the unweighted niches are plotted.
#' @param outcomeCol Name of the per-image outcome column in `colData(spe)`.
#'
#' @return A `ggplot` object; for `nicheProportionBoxplot()` a `girafe`
#'   (interactive) widget; for `nicheCompositionMap(type = "heatmap")` the
#'   result of [pheatmap::pheatmap()].
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 10,
#'     seed = 1)
#' niches = runNicheIdentification(spe, heartCellTypeTree, confList,
#'     r = c(20, 50), k = 5)
#' spe$nicheWeighted = niches$nicheWeighted[spe$cellID]
#' plotNicheSpatial(spe, unique(spe$imageID)[1], "niche1")
#' nicheCompositionMap(spe, "nicheWeighted")
#' nicheBarplot(spe, "nicheWeighted", "niche1")
#' if (requireNamespace("ggiraph", quietly = TRUE)) {
#'     res = testNicheOutcome(spe, niches$nicheUnweighted,
#'         niches$nicheWeighted, outcomeCol = "resolution_broad",
#'         subject = "subject")
#'     nicheProportionBoxplot(res$propWeighted, "niche1", spe, "resolution_broad")
#'     nicheProportionBoxplot(res$propUnweighted, "niche1", spe,
#'         "resolution_broad")
#' }
#' }
#' @param cellSize Size of the highlighted cells in `plotNicheSpatial()`, the
#'   cells of the chosen niche. The grey background cells always have size 0.7.
#' @name nichePlots
#' @export
plotNicheSpatial = function(spe, imageID, niche, nicheCol = "nicheWeighted", cellSize = 1) {
    checkCellSize(cellSize)
    checkSpe(spe)
    cells = getCellData(spe)
    checkColumn(cells, nicheCol, "nicheCol")
    checkImage(cells, imageID)
    checkNiche(cells, niche, nicheCol)
    imgCells = cells |> dplyr::filter(imageID == !!imageID)

    ggplot() +
        geom_point(
            data = imgCells |> dplyr::filter(is.na(.data[[nicheCol]]) | .data[[nicheCol]] != niche),
            aes(x = x, y = y), colour = "grey80", size = 0.7
        ) +
        geom_point(
            data = imgCells |> dplyr::filter(.data[[nicheCol]] == niche),
            aes(x = x, y = y, colour = cellType), size = cellSize
        ) +
        coord_equal() +
        scale_y_reverse() +
        labs(title = paste0(imageID, " - ", niche), colour = "Cell type", x = "", y = "") +
        theme_classic()
}

#' @rdname nichePlots
#' @export
nicheCompositionMap = function(spe, niche, cellType = "cellType", type = c("bubble", "heatmap"),
                                limit = c(0.33, 3), ...) {
    type = match.arg(type)
    checkSpe(spe)
    cells = getCellData(spe)
    checkColumn(cells, cellType, "cellType")
    checkColumn(cells, niche, "niche")
    df = cells[, c(cellType, niche)]

    tab = table(df[[cellType]], df[[niche]])
    tab = tab / rowSums(tab) %*% t(colSums(tab)) * sum(tab)
    tabClipped = pmax(pmin(tab, limit[2]), limit[1])

    ph = pheatmap::pheatmap(tabClipped, cluster_cols = FALSE, silent = TRUE, ...)

    if (type == "heatmap") {
        return(pheatmap::pheatmap(tabClipped, cluster_cols = FALSE, ...))
    }

    tab |>
        as.data.frame() |>
        dplyr::mutate(
            cellType = factor(Var1, levels = levels(Var1)[ph$tree_row$order]),
            niche = Var2,
            enrichment = pmax(pmin(Freq, limit[2]), limit[1])
        ) |>
        ggplot(aes(x = niche, y = cellType, colour = enrichment, size = enrichment)) +
        geom_point() +
        scale_colour_gradient2(low = "#4575B4", mid = "grey90", high = "#D73027", midpoint = 1, guide = "legend") +
        theme_minimal() +
        labs(x = "Niche", y = "Cell type", colour = "Relative\nenrichment", size = "Relative\nenrichment") +
        theme(axis.text.x = element_text(angle = 90))
}

#' @rdname nichePlots
#' @export
nicheBarplot = function(spe, niche, target, cellType = "cellType", limit = c(0.33, 3)) {
    checkSpe(spe)
    cells = getCellData(spe)
    checkColumn(cells, cellType, "cellType")
    checkColumn(cells, niche, "niche")
    df = cells[, c(cellType, niche)]

    tab = table(df[[cellType]], df[[niche]])
    tab = tab / rowSums(tab) %*% t(colSums(tab)) * sum(tab) # same relative-enrichment formula as nicheCompositionMap
    if (!is.character(target) || length(target) != 1 || !target %in% colnames(tab)) {
        stop("`target` must be one of the niches in `", niche, "`: ",
            paste(colnames(tab), collapse = ", "), call. = FALSE)
    }

    tab[, target] |>
        as.data.frame() |>
        tibble::rownames_to_column("cellType") |>
        setNames(c("cellType", "enrichment")) |>
        dplyr::mutate(
            enrichment = pmax(pmin(enrichment, limit[2]), limit[1]),
            cellType = forcats::fct_reorder(cellType, enrichment)
        ) |>
        ggplot(aes(x = cellType, y = enrichment)) +
        geom_col(fill = "#2a78d6", width = 0.7) +
        geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
        coord_flip() +
        labs(x = "Cell type", y = "Relative enrichment", title = paste0(target, " (", niche, ")")) +
        theme_classic()
}

#' @rdname nichePlots
#' @export
nicheProportionBoxplot = function(proportions, niche, spe, outcomeCol) {
    checkSpe(spe)
    p = nicheProportionBoxplotGG(proportions, niche, getCellData(spe), outcomeCol)
    ggiraph::girafe(ggobj = p, height_svg = 3.5)
}

# the ggplot behind nicheProportionBoxplot(), before it is made interactive
nicheProportionBoxplotGG = function(proportions, niche, cells, outcomeCol) {
    if (!is.data.frame(proportions) ||
        !all(c("imageID", "niche", "proportion") %in% names(proportions))) {
        stop("`proportions` must be the `propWeighted` or `propUnweighted` ",
            "element of the results of testNicheOutcome()", call. = FALSE)
    }
    if (!is.character(niche) || length(niche) != 1 || !niche %in% proportions$niche) {
        stop("`niche` must be one of the niches in `proportions`: ",
            paste(sort(unique(proportions$niche)), collapse = ", "), call. = FALSE)
    }
    checkColumn(cells, outcomeCol, "outcomeCol")

    outcomeDat = dplyr::distinct(cells, imageID, .data[[outcomeCol]])
    plotDat = proportions |>
        dplyr::filter(niche == !!niche) |>
        dplyr::inner_join(outcomeDat, by = "imageID") |>
        dplyr::mutate(!!outcomeCol := as.factor(.data[[outcomeCol]]))

    # testNicheOutcome() records whether the table is weighted or unweighted
    type = attr(proportions, "nicheType")

    ggplot(plotDat, aes(x = .data[[outcomeCol]], y = proportion, colour = .data[[outcomeCol]])) +
        geom_boxplot(outlier.shape = NA, width = 0.5) +
        ggiraph::geom_point_interactive(aes(tooltip = imageID, data_id = imageID),
            position = position_jitter(width = 0.1), alpha = 0.7, size = 1.5
        ) +
        labs(
            x = NULL, y = paste0(niche, " proportion"), colour = NULL,
            title = paste0(niche, if (!is.null(type)) paste0(" (", type, ")"))
        ) +
        theme_classic()
}
