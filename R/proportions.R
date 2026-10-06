# Developer-level building blocks for runProportionsPipeline() (not exported).
# `cells` is a per-cell data frame with imageID, cellID and cellType columns, plus
# `group` for the group-level function.
#
# calcUnweightedProportions(): images x cell types proportions from the original
# (hard) labels.
calcUnweightedProportions = function(cells) {
    tab = table(sample = cells$imageID, clusters = cells$cellType)
    tab / rowSums(tab)
}

# Unweighted proportions of every group per image. Group membership is a hard,
# deterministic classification (unlike a subtype within a group, which carries
# bootstrap confidence), so there is no weighted counterpart to this one.
calcUnweightedGroupProportions = function(cells) {
    tab = table(sample = cells$imageID, clusters = cells$group)
    tab / rowSums(tab)
}

# A linear model of each column of `props` (images x cell types) against the
# outcome, one per cell type. `estimate` is the second level of the outcome
# minus the first. The outcome is checked first: two levels, one value per image
# and, if `subject` is given, one image per subject.
runProportionModels = function(props, cells, outcomeCol, subject = NULL) {
    checkOutcome(cells, outcomeCol, subject)
    outcomeDat = dplyr::distinct(cells, imageID, .data[[outcomeCol]])

    propDat = data.frame(imageID = rownames(props), unclass(props), check.names = FALSE) |>
        dplyr::left_join(outcomeDat, by = "imageID")

    lapply(colnames(props), function(cellType) {
        dat = propDat |> dplyr::filter(!is.na(.data[[cellType]]), !is.na(.data[[outcomeCol]]))

        # skip cell types with too little data to fit a meaningful model
        if (length(unique(dat[[outcomeCol]])) < 2 || nrow(dat) < 5) {
            return(data.frame(cellType = cellType, estimate = NA, statistic = NA, p.value = NA))
        }

        fit = tryCatch(lm(as.formula(paste0("`", cellType, "` ~ `", outcomeCol, "`")), data = dat),
            error = function(e) NULL
        )
        if (is.null(fit)) {
            return(data.frame(cellType = cellType, estimate = NA, statistic = NA, p.value = NA))
        }

        coef = summary(fit)$coefficients
        data.frame(
            cellType = cellType, estimate = coef[2, "Estimate"],
            statistic = coef[2, "t value"], p.value = coef[2, "Pr(>|t|)"]
        )
    }) |>
        dplyr::bind_rows() |>
        dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))
}

propsToLong = function(props, valueCol) {
    dplyr::rename(
        data.frame(
            imageID = rep(rownames(props), times = ncol(props)),
            cellType = rep(colnames(props), each = nrow(props)),
            value = as.vector(props),
            stringsAsFactors = FALSE
        ),
        !!valueCol := value
    )
}

# Weighted proportion of every subtype of every group in confList, per image:
# the sum of the confidences of the group's cells for that subtype, divided by
# all cells in the image.
runWeightedProportions = function(cells, cellTypeTree, confList, skipGroups = character(0)) {
    cells = annotateCellTypeTree(withoutAnnotation(cells), cellTypeTree)

    # these do not depend on the group or subtype, so they are computed once
    # instead of once per (group, subtype) pair
    trueTotal = table(cells$imageID)
    allImages = names(trueTotal)
    rowsByGroup = split(seq_len(nrow(cells)), cells$group)

    pairs = dplyr::bind_rows(lapply(setdiff(names(confList), skipGroups), function(group) {
        data.frame(
            group = group, subtype = setdiff(colnames(confList[[group]]), "cellID"),
            stringsAsFactors = FALSE
        )
    }))

    weightedList = lapply(seq_len(nrow(pairs)), function(i) {
        group = pairs$group[i]
        subtype = pairs$subtype[i]

        # every member of the group counts towards `subtype` with weight equal
        # to its bootstrap confidence of actually being that subtype
        inGroupRows = rowsByGroup[[group]]
        if (is.null(inGroupRows)) inGroupRows = integer(0)
        conf = confList[[group]]
        groupWeight = conf[[subtype]][match(cells$cellID[inGroupRows], conf$cellID)]

        subtypeSum = tapply(groupWeight, cells$imageID[inGroupRows], sum)

        full = setNames(rep(0, length(allImages)), allImages)
        full[names(subtypeSum)] = subtypeSum

        data.frame(
            imageID = allImages, lineage = group, cellType = subtype,
            proportionWeighted = unname(full / as.vector(trueTotal[allImages])),
            stringsAsFactors = FALSE
        )
    })

    dplyr::bind_rows(weightedList)
}

#' Cell type proportion analysis
#'
#' Computes per-image proportions of every cell type using original labels
#' (unweighted) and bootstrap confidences (weighted), and tests each against an
#' outcome with a linear model. `runProportionsRelativeToParentPipeline()`
#' instead expresses each subtype as a proportion of its parent group.
#'
#' @details
#' In `runProportionsPipeline()` the unweighted results cover every cell type
#' and also the whole groups that were bootstrapped. Only the subtypes of those
#' groups have confidences, so for every other cell type (cell types without a
#' group, and the group-level proportions) the weighted proportion is the
#' unweighted one: their labels are not uncertain. The weighted and unweighted
#' models therefore test the same set of cell types, and their adjusted
#' p-values are comparable.
#'
#' The `lineage` column of the results is the group a subtype belongs to (not
#' the broad lineage). For cell types without a group and for the group-level
#' proportions it is the cell type or group itself.
#'
#' @param spe A `SummarizedExperiment` (e.g. `SpatialExperiment`) with a cell
#'   ID, cell type and image for every cell in `colData`.
#' @param cellTypeTree A cell type tree, see [cellTypeTree].
#' @param confList Confidence list from [bootstrapLineageConfidence()]: a named
#'   list with one data frame per group, whose names must be groups of
#'   `cellTypeTree`.
#' @param outcomeCol Name of the per-image outcome column in `colData(spe)`. It
#'   must have exactly two levels (subset your data to two conditions first)
#'   and a single value per image: repeated measures and mixed models are not
#'   supported. An error is raised otherwise, before any computation.
#' @param skipGroups Names of groups of `cellTypeTree` to exclude, as a character
#'   vector (the default `character(0)` excludes none). An error is raised for
#'   a name that is not a group of the tree, so the same `skipGroups` can be used
#'   here and in [bootstrapLineageConfidence()].
#' @param subject Name of the `colData(spe)` column identifying the subject (for
#'   example the patient) each image comes from, or `NULL`. If given, every
#'   image must belong to one subject and every subject must contribute only
#'   one image, as in [runColocalisationModels()]. With `NULL` this cannot be checked,
#'   and a message says so.
#' @param cellType,cellID,imageID Names of the `colData(spe)` columns holding
#'   each cell's cell type label, a unique cell ID (matching the `cellID` column
#'   of `confList`), and the image it is in.
#'
#' @return A list with `proportions`, `weighted`, `unweighted` (long
#'   `data.frame`s; only `proportions` and not `weighted`/`unweighted` for
#'   the relative-to-parent version), and `modelWeighted` and
#'   `modelUnweighted` (per cell type model results, where `estimate` is the
#'   difference for the second level of the outcome compared with the first).
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 10,
#'     seed = 1)
#' res = runProportionsPipeline(spe, heartCellTypeTree, confList,
#'     outcomeCol = "resolution_broad", subject = "subject")
#' head(res$modelWeighted)
#' }
#' @export
runProportionsPipeline = function(spe, cellTypeTree, confList, outcomeCol, skipGroups = character(0),
                                  subject = NULL, cellType = "cellType", cellID = "cellID",
                                  imageID = "imageID") {
    # check the arguments before doing any work
    checkProportionArgs(spe, cellTypeTree, confList, skipGroups)
    checkColumn(spe, outcomeCol, "outcomeCol")
    if (!is.null(subject)) checkColumn(spe, subject, "subject")

    cells = standardiseCells(getCellData(spe), list(cellID = cellID, cellType = cellType,
        imageID = imageID))
    cells = annotateCellTypeTree(cells, cellTypeTree)
    cells$group[is.na(cells$group)] = cells$cellType[is.na(cells$group)]

    checkOutcome(cells, outcomeCol, subject)
    messageIfNoSubject(subject)

    message("Computing unweighted proportions once (shared across all subtypes)...")
    unweightedProps = calcUnweightedProportions(cells)

    # bundle in whole-group proportions for every bootstrapped group (one that
    # actually has confList subtypes). Other groups are left out: for an
    # ungrouped cell type the group equals the cell type, so its group-level
    # column would just duplicate the cell type one already present.
    bootstrappedGroups = setdiff(names(confList), skipGroups)
    groupProps = calcUnweightedGroupProportions(cells)
    groupProps = groupProps[rownames(unweightedProps), bootstrappedGroups, drop = FALSE]
    unweightedProps = cbind(unweightedProps, groupProps)

    unweightedLong = propsToLong(unweightedProps, "proportionUnweighted")
    unweightedModels = runProportionModels(unweightedProps, cells, outcomeCol, subject)

    message("Computing weighted proportions...")
    weighted = runWeightedProportions(cells, cellTypeTree, confList, skipGroups = skipGroups)

    # only subtypes of bootstrapped groups have confidences. Everything else
    # (ungrouped cell types and whole groups) has hard labels, so its weighted
    # proportion is the unweighted one. This puts both sets of models on the
    # same cell types, so their adjusted p-values are comparable.
    notWeighted = unweightedLong[!unweightedLong$cellType %in% weighted$cellType, ]
    weighted = dplyr::bind_rows(weighted, data.frame(
        imageID = notWeighted$imageID, lineage = notWeighted$cellType,
        cellType = notWeighted$cellType,
        proportionWeighted = notWeighted$proportionUnweighted,
        stringsAsFactors = FALSE
    ))

    # need to pivot for runProportionModels
    weightedWide = weighted |>
        tidyr::pivot_wider(id_cols = imageID, names_from = cellType, values_from = proportionWeighted) |>
        tibble::column_to_rownames("imageID") |>
        as.matrix()

    message("Fitting weighted models...")
    modelWeighted = runProportionModels(weightedWide, cells, outcomeCol, subject) |>
        dplyr::left_join(dplyr::distinct(weighted, lineage, cellType), by = "cellType")

    proportions = dplyr::inner_join(unweightedLong, weighted, by = c("imageID", "cellType"))

    list(
        proportions     = proportions,
        weighted        = weighted,
        unweighted      = unweightedLong,
        modelWeighted   = modelWeighted,
        modelUnweighted = unweightedModels
    )
}

# Proportion of each subtype within its parent group, per image, both unweighted
# (original labels) and weighted (bootstrap confidences). A long data.frame with
# the group in `lineage`; images with no cells of the group get NaN.
calcProportionsRelativeToParent = function(cells, cellTypeTree, confList, skipGroups = character(0)) {
    cells = annotateCellTypeTree(withoutAnnotation(cells), cellTypeTree)
    allImages = unique(cells$imageID)

    rowsByGroup = split(seq_len(nrow(cells)), cells$group)

    # map parents to child
    pairs = dplyr::bind_rows(lapply(setdiff(names(confList), skipGroups), function(group) {
        data.frame(
            group = group, subtype = setdiff(colnames(confList[[group]]), "cellID"),
            stringsAsFactors = FALSE
        )
    }))

    results = lapply(seq_len(nrow(pairs)), function(i) {
        group = pairs$group[i]
        subtype = pairs$subtype[i]

        inGroupRows = rowsByGroup[[group]]
        if (is.null(inGroupRows)) inGroupRows = integer(0)
        groupImages = cells$imageID[inGroupRows]

        # denominator: hard count of cells in the parent group, per image
        parentTotal = setNames(rep(0, length(allImages)), allImages)
        parentTab = table(groupImages)
        parentTotal[names(parentTab)] = parentTab

        # unweighted numerator: hard count of cells actually labelled `subtype`
        unweightedCount = setNames(rep(0, length(allImages)), allImages)
        hardTab = table(groupImages[cells$cellType[inGroupRows] == subtype])
        unweightedCount[names(hardTab)] = hardTab

        # weighted numerator
        conf = confList[[group]]
        weight = conf[[subtype]][match(cells$cellID[inGroupRows], conf$cellID)]
        weightedSum = setNames(rep(0, length(allImages)), allImages)
        weightedTab = tapply(weight, groupImages, sum)
        weightedSum[names(weightedTab)] = weightedTab

        data.frame(
            imageID = allImages,
            lineage = group,
            cellType = subtype,
            proportionUnweighted = unname(unweightedCount) / unname(parentTotal),
            proportionWeighted = unname(weightedSum) / unname(parentTotal),
            stringsAsFactors = FALSE
        )
    })

    dplyr::bind_rows(results)
}

#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 10,
#'     seed = 1)
#' rel = runProportionsRelativeToParentPipeline(spe, heartCellTypeTree,
#'     confList, outcomeCol = "resolution_broad", subject = "subject")
#' head(rel$modelWeighted)
#' }
#' @rdname runProportionsPipeline
#' @export
runProportionsRelativeToParentPipeline = function(spe, cellTypeTree, confList, outcomeCol,
                                                   skipGroups = character(0), subject = NULL,
                                                   cellType = "cellType", cellID = "cellID",
                                                   imageID = "imageID") {
    # check the arguments before doing any work
    checkProportionArgs(spe, cellTypeTree, confList, skipGroups)
    checkColumn(spe, outcomeCol, "outcomeCol")
    if (!is.null(subject)) checkColumn(spe, subject, "subject")

    cells = standardiseCells(getCellData(spe), list(cellID = cellID, cellType = cellType,
        imageID = imageID))
    cells = annotateCellTypeTree(cells, cellTypeTree)
    cells$group[is.na(cells$group)] = cells$cellType[is.na(cells$group)]

    checkOutcome(cells, outcomeCol, subject)
    messageIfNoSubject(subject)

    message("Computing proportions relative to parent...")
    proportions = calcProportionsRelativeToParent(cells, cellTypeTree, confList, skipGroups = skipGroups)

    unweightedWide = proportions |>
        tidyr::pivot_wider(id_cols = imageID, names_from = cellType, values_from = proportionUnweighted) |>
        tibble::column_to_rownames("imageID") |>
        as.matrix()

    weightedWide = proportions |>
        tidyr::pivot_wider(id_cols = imageID, names_from = cellType, values_from = proportionWeighted) |>
        tibble::column_to_rownames("imageID") |>
        as.matrix()

    message("Fitting unweighted models...")
    modelUnweighted = runProportionModels(unweightedWide, cells, outcomeCol, subject) |>
        dplyr::left_join(dplyr::distinct(proportions, lineage, cellType), by = "cellType")

    message("Fitting weighted models...")
    modelWeighted = runProportionModels(weightedWide, cells, outcomeCol, subject) |>
        dplyr::left_join(dplyr::distinct(proportions, lineage, cellType), by = "cellType")

    list(
        proportions = proportions,
        modelWeighted = modelWeighted,
        modelUnweighted = modelUnweighted
    )
}
