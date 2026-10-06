# Internal helpers shared across the package (not exported).

# Per-cell data.frame (one row per cell) from a SummarizedExperiment /
# SpatialExperiment.
# Pulls colData; for a SpatialExperiment whose coordinates live in
# spatialCoords() rather
# than colData, adds them as columns `x` and `y`.
#' Per-cell data frame from a SummarizedExperiment
#'
#' Returns `colData(spe)` as a `data.frame`. For a `SpatialExperiment` whose
#' coordinates are stored in `spatialCoords()` rather than in `colData()`,
#' the coordinates are added as columns `x` and `y`. The spatial functions in
#' scBoot expect `cellID`, `imageID`, `cellType`, `x` and `y` columns.
#'
#' @param spe A `SummarizedExperiment` or `SpatialExperiment`.
#'
#' @return A `data.frame` with one row per cell.
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' head(getCellData(spe))
#' }
#' @export
getCellData = function(spe) {
    cells = as.data.frame(SummarizedExperiment::colData(spe))
    if (methods::is(spe, "SpatialExperiment") && !all(c("x", "y") %in% names(cells))) {
        xy = SpatialExperiment::spatialCoords(spe)
        cells$x = xy[, 1]
        cells$y = xy[, 2]
    }
    cells
}

# generate all possible pairwise comparisons, one row per pair
generateComparisons = function(cellTypeTree, confList) {
    groups = names(confList)
    cellMap = buildCellTypeMap(cellTypeTree)$cellMap

    # ungrouped leaf nodes - their cell type label is the group
    cellMap$group[is.na(cellMap$group)] = cellMap$cellType[is.na(cellMap$group)]
    leaves = cellMap$cellType

    # group x group - every ordered pair, both directions
    groupGroup = expand.grid(from = groups, to = groups, stringsAsFactors = FALSE)

    # group x leaf - every ordered pair, both directions
    groupLeafBase = expand.grid(group = groups, leaf = leaves, stringsAsFactors = FALSE) |>
        dplyr::mutate(leafGroup = cellMap$group[match(leaf, cellMap$cellType)])

    groupLeaf = dplyr::bind_rows(
        dplyr::mutate(groupLeafBase, from = group, to = leaf),
        dplyr::mutate(groupLeafBase, from = leaf, to = group)
    )

    # leaf x leaf - every ordered pair, both directions
    leafLeaf = expand.grid(from = leaves, to = leaves, stringsAsFactors = FALSE)

    list(groupGroup = groupGroup, groupLeaf = groupLeaf, leafLeaf = leafLeaf)
}

# L-value calculation
Lfunction = function(numerator, r, area, sumWeightsChild1, sumWeightsChild2) {
    lambda2 = sumWeightsChild2 / area
    kvalue = (numerator / lambda2) * (1 / sumWeightsChild1)
    sqrt(kvalue / pi) - r
}

# log1p-normalizes and z-scores every row (cell) of a cells x genes count matrix
# in one
# vectorized pass. Cells with zero total count come out all-zero.
scaleCountMat = function(mat) {
    totals = rowSums(mat)
    norm = log1p(mat / pmax(totals, 1) * 1e4)
    norm[totals == 0, ] = 0

    means = rowMeans(norm)
    centred = norm - means
    sds = sqrt(rowSums(centred^2) / (ncol(mat) - 1))
    sds[sds == 0] = 1 # avoid div-by-zero for zero-variance rows
    centred / sds
}

# Edge-correction weights for every point in a ppp: the fraction of the
# radius-`maxD`
# disc around each point that lies inside the observation window (1 for interior
# points).
# Copied from Statial:::.borderEdge (Statial, GPL-3,
# https://github.com/SydneyBioX/Statial)
# so we don't depend on an unexported function.
borderEdge = function(X, maxD) {
    W = X$window
    bW = spatstat.geom::union.owin(
        spatstat.geom::border(W, maxD, outside = FALSE),
        spatstat.geom::border(W, 2, outside = TRUE)
    )
    inB = spatstat.geom::inside.owin(X$x, X$y, bW)
    e = rep(1, X$n)
    if (any(inB)) {
        circs = spatstat.geom::discs(X[inB], maxD, separate = TRUE)
        circs = spatstat.geom::solapply(circs, spatstat.geom::intersect.owin, X$window)
        areas = unlist(lapply(circs, spatstat.geom::area)) / (pi * maxD^2)
        e[inB] = areas
    }
    e
}

# names of the per-cell columns of a SummarizedExperiment (its colData) or of a
# data.frame / DataFrame
cellDataNames = function(x) {
    if (methods::is(x, "SummarizedExperiment")) {
        names(SummarizedExperiment::colData(x))
    } else {
        names(x)
    }
}

# checks that `value`, the argument called `arg`, is a single string naming a
# column of x (colData of a SummarizedExperiment, or a data.frame)
checkColumn = function(x, value, arg) {
    where = if (methods::is(x, "SummarizedExperiment")) "colData(spe)" else "the data"
    if (!is.character(value) || length(value) != 1 || is.na(value)) {
        stop("`", arg, "` must be a single string naming a column of ", where,
            call. = FALSE)
    }
    available = cellDataNames(x)
    if (!value %in% available) {
        stop("`", arg, " = \"", value, "\"` is not a column of ", where,
            ". Available columns: ", paste(utils::head(available, 10), collapse = ", "),
            if (length(available) > 10) ", ...", call. = FALSE)
    }
    invisible(value)
}

# checks that `value`, the argument called `arg`, is a single string naming an
# assay of the SummarizedExperiment `spe`
checkAssay = function(spe, value, arg = "assayName") {
    if (!is.character(value) || length(value) != 1 || is.na(value)) {
        stop("`", arg, "` must be a single string naming an assay of spe",
            call. = FALSE)
    }
    available = SummarizedExperiment::assayNames(spe)
    if (!value %in% available) {
        stop("`", arg, " = \"", value, "\"` is not an assay of spe. ",
            "Available assays: ", paste(available, collapse = ", "), call. = FALSE)
    }
    invisible(value)
}

# checks the per-image outcome used by the co-localisation models: the column
# exists, has a single value per image, and has exactly two levels. If `subject`
# is given, also checks that every image comes from a different subject.
checkOutcome = function(cells, outcomeCol, subject = NULL) {
    checkColumn(cells, outcomeCol, "outcomeCol")
    checkColumn(cells, "imageID", "imageID")
    if (!is.null(subject)) checkSubjects(cells, subject)
    outcome = cells[[outcomeCol]]

    # one value per image: each image is one independent observation, so
    # repeated measures / mixed models are not supported
    perImage = split(outcome, cells$imageID)
    valuesPerImage = lapply(perImage, function(v) unique(v[!is.na(v)]))
    mixed = names(valuesPerImage)[lengths(valuesPerImage) > 1]
    if (length(mixed) > 0) {
        shown = utils::head(mixed, 5)
        examples = vapply(shown, function(id) {
            sprintf("%s (%s)", id, paste(valuesPerImage[[id]], collapse = ", "))
        }, character(1))
        stop("`outcomeCol = \"", outcomeCol, "\"` must have a single value per ",
            "image, but ", length(mixed), " image", if (length(mixed) > 1) "s",
            " ", if (length(mixed) > 1) "have" else "has", " more than one: ",
            paste(examples, collapse = "; "), if (length(mixed) > 5) "; ...",
            ". The models treat each image as one independent observation, so ",
            "repeated measures and mixed models are not supported. Please make ",
            "sure there is exactly one outcome value per image, and one image ",
            "per sample (for example one image per patient).", call. = FALSE)
    }

    # exactly two conditions
    levels = sort(unique(outcome[!is.na(outcome)]))
    shownLevels = paste(utils::head(levels, 6), collapse = ", ")
    if (length(levels) > 2) {
        stop("`outcomeCol = \"", outcomeCol, "\"` has ", length(levels),
            " levels (", shownLevels, if (length(levels) > 6) ", ...", "). The ",
            "models currently compare two conditions only. Please subset the ",
            "data to two conditions before running.", call. = FALSE)
    }
    if (length(levels) < 2) {
        stop("`outcomeCol = \"", outcomeCol, "\"` has fewer than two levels",
            if (length(levels) == 1) paste0(" (only ", shownLevels, ")"),
            ". Two conditions are needed to compare.", call. = FALSE)
    }
    invisible(levels)
}

# images whose smaller side is below 2 * r, i.e. too small to hold a
# neighbourhood of radius r. Returns a data.frame (imageID, xRange, yRange), one
# row per such image.
smallImages = function(cells, r) {
    xRange = tapply(cells$x, cells$imageID, function(v) diff(range(v)))
    yRange = tapply(cells$y, cells$imageID, function(v) diff(range(v)))
    small = which(pmin(xRange, yRange) < 2 * r)
    data.frame(imageID = names(xRange)[small], xRange = as.numeric(xRange[small]),
        yRange = as.numeric(yRange[small]), stringsAsFactors = FALSE)
}

# tells the user which images are being dropped as too small for radius r, and
# stops if that is every image
reportSmallImages = function(small, nImages, r) {
    if (nrow(small) == 0) {
        return(invisible())
    }
    if (nrow(small) == nImages) {
        stop("All ", nImages, " images are smaller than 2 x r = ", 2 * r,
            " units in at least one direction. Use a smaller `r`.", call. = FALSE)
    }
    shown = utils::head(small, 10)
    message("Dropping ", nrow(small), " of ", nImages, " images smaller than ",
        "2 x r = ", 2 * r, " units in at least one direction: ",
        paste0(shown$imageID, " (", round(shown$xRange), " x ", round(shown$yRange), ")",
            collapse = ", "), if (nrow(small) > 10) ", ...")
    invisible()
}

# removes images that are too small for radius r from the per-cell data, with a
# message listing them
dropSmallImages = function(cells, r) {
    small = smallImages(cells, r)
    reportSmallImages(small, length(unique(cells$imageID)), r)
    cells[!cells$imageID %in% small$imageID, , drop = FALSE]
}

# checks that images and subjects map one-to-one: every image belongs to a
# single subject, and every subject contributes a single image. Several images
# from one subject are not independent, and the models treat every image as an
# independent observation (they are not mixed models).
checkSubjects = function(cells, subject) {
    checkColumn(cells, subject, "subject")
    checkColumn(cells, "imageID", "imageID")
    ids = as.character(cells[[subject]])
    if (anyNA(ids)) {
        noSubject = unique(as.character(cells$imageID)[is.na(ids)])
        stop("The `subject` column '", subject, "' has missing values (in ",
            length(noSubject), " image", if (length(noSubject) > 1) "s", ": ",
            paste(utils::head(noSubject, 5), collapse = ", "),
            if (length(noSubject) > 5) ", ...",
            "). Every cell needs a subject so that images from the same ",
            "subject can be detected.", call. = FALSE)
    }

    map = unique(data.frame(imageID = as.character(cells$imageID), subject = ids,
        stringsAsFactors = FALSE))
    describe = function(by, of, label) {
        groups = split(map[[of]], map[[by]])
        groups = groups[lengths(groups) > 1]
        shown = utils::head(names(groups), 5)
        list(n = length(groups), text = paste0(vapply(shown, function(g) {
            sprintf("%s (%s)", g, paste(utils::head(groups[[g]], 4), collapse = ", "))
        }, character(1)), collapse = "; "), more = length(groups) > 5)
    }

    # every image must belong to one subject
    perImage = describe("imageID", "subject")
    if (perImage$n > 0) {
        stop("Every image must belong to a single subject, but ", perImage$n,
            " image", if (perImage$n > 1) "s", " ", if (perImage$n > 1) "have" else "has",
            " more than one value in the `subject` column '", subject, "': ",
            perImage$text, if (perImage$more) "; ...", ".", call. = FALSE)
    }

    # every subject must contribute one image
    perSubject = describe("subject", "imageID")
    if (perSubject$n > 0) {
        stop("`subject = \"", subject, "\"` has more than one image for ",
            perSubject$n, " subject", if (perSubject$n > 1) "s", ": ",
            perSubject$text, if (perSubject$more) "; ...",
            ". Images from the same subject are not independent, and the models ",
            "treat every image as an independent observation (mixed models are ",
            "not supported). Please keep one image per subject, for example by ",
            "subsetting to one image per patient.", call. = FALSE)
    }
    invisible(TRUE)
}

# A copy of the per-cell data with the standard columns the analysis code reads
# (cellID, cellType, imageID, x, y) filled from the columns the user named.
# `cols` is a named list: names are the standard columns, values the user's
# column names (each checked to be a single string). Checks that every column
# exists, that the coordinates are numeric and complete, and that cell IDs are
# unique.
standardiseCells = function(cells, cols) {
    for (std in names(cols)) checkColumn(cells, cols[[std]], std)
    for (std in names(cols)) cells[[std]] = cells[[cols[[std]]]]

    for (coord in intersect(c("x", "y"), names(cols))) {
        if (!is.numeric(cells[[coord]]) || anyNA(cells[[coord]])) {
            stop("`", coord, " = \"", cols[[coord]], "\"` must be a numeric ",
                "column without missing values", call. = FALSE)
        }
    }
    if ("cellID" %in% names(cols) &&
        (anyNA(cells$cellID) || anyDuplicated(cells$cellID))) {
        stop("The `cellID` column '", cols[["cellID"]], "' must contain unique, ",
            "non-missing cell IDs", if (anyDuplicated(cells$cellID))
                paste0(" (e.g. '", cells$cellID[anyDuplicated(cells$cellID)],
                    "' is repeated)"), call. = FALSE)
    }
    cells
}

# checks that `k`, the number of clusters, is a single whole number of at least 2
checkNumberOfClusters = function(k, arg = "k") {
    if (!is.numeric(k) || length(k) != 1 || is.na(k) || k != round(k) || k < 2) {
        stop("`", arg, "` (the number of niches) must be a single whole number ",
            "of at least 2", call. = FALSE)
    }
    invisible(k)
}

# checks that `r` is a non-empty vector of positive, finite radii
checkRadii = function(r, arg = "r") {
    if (!is.numeric(r) || length(r) == 0 || anyNA(r) || any(!is.finite(r)) || any(r <= 0)) {
        stop("`", arg, "` must be a vector of one or more positive numbers (radii)",
            call. = FALSE)
    }
    invisible(r)
}

# checks that `seed` is NULL (no fixed seed) or a single number
checkSeed = function(seed, arg = "seed") {
    if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || is.na(seed))) {
        stop("`", arg, "` must be NULL or a single number", call. = FALSE)
    }
    invisible(seed)
}

# checks that every label is a cell type or a group of the tree
checkTreeLabels = function(labels, cellTypeTree, arg = "toLabels") {
    treeInfo = buildCellTypeMap(cellTypeTree)
    valid = union(treeInfo$cellMap$cellType, treeInfo$groupsWithChildren)
    if (!is.character(labels) || length(labels) == 0 || anyNA(labels)) {
        stop("`", arg, "` must be NULL or a character vector of cell types and ",
            "groups in cellTypeTree", call. = FALSE)
    }
    unknown = setdiff(labels, valid)
    if (length(unknown) > 0) {
        stop("`", arg, "` has ", length(unknown), " label",
            if (length(unknown) > 1) "s", " that ", if (length(unknown) > 1) "are" else "is",
            " not a cell type or group in cellTypeTree: ",
            paste0("'", utils::head(unknown, 5), "'", collapse = ", "),
            if (length(unknown) > 5) ", ...", ". Available labels: ",
            paste(utils::head(valid, 10), collapse = ", "),
            if (length(valid) > 10) ", ...", call. = FALSE)
    }
    invisible(labels)
}

# per-cell data without any group / lineage columns from an earlier annotation,
# so that annotating it again does not report a replaced column
withoutAnnotation = function(cells) {
    cells[, !names(cells) %in% c("group", "lineage"), drop = FALSE]
}

# one-line reminder that images from one subject cannot be detected when no
# `subject` column was given
messageIfNoSubject = function(subject) {
    if (is.null(subject)) {
        message("`subject` was not given, so it is not checked that every image ",
            "comes from a different subject. Several images from one subject ",
            "would be treated as independent samples.")
    }
    invisible()
}

# checks the arguments the proportion pipelines share: `spe` is a
# SummarizedExperiment, `confList` has the structure that
# bootstrapLineageConfidence() returns for groups of the tree, and `skipGroups`
# names groups of the tree (the tree itself, not confList, so the same
# skipGroups can be passed here and to bootstrapLineageConfidence()).
checkProportionArgs = function(spe, cellTypeTree, confList, skipGroups) {
    checkSpe(spe)
    # also validates the tree
    groups = buildCellTypeMap(cellTypeTree)$groupsWithChildren
    showGroups = function() paste(utils::head(groups, 10), collapse = ", ")

    if (!is.list(confList) || is.data.frame(confList) || length(confList) == 0 ||
        is.null(names(confList)) || any(!nzchar(names(confList)))) {
        stop("`confList` must be a named list with one data.frame per group, as ",
            "returned by bootstrapLineageConfidence()", call. = FALSE)
    }
    for (g in names(confList)) {
        conf = confList[[g]]
        if (!is.data.frame(conf) || !"cellID" %in% names(conf) || ncol(conf) < 2) {
            stop("`confList[[\"", g, "\"]]` must be a data.frame with a `cellID` ",
                "column and one column of confidence scores per subtype, as ",
                "returned by bootstrapLineageConfidence()", call. = FALSE)
        }
    }
    notGroups = setdiff(names(confList), groups)
    if (length(notGroups) > 0) {
        stop("`confList` has ", length(notGroups), " group", if (length(notGroups) > 1) "s",
            " that ", if (length(notGroups) > 1) "are" else "is", " not a group in ",
            "cellTypeTree: ", paste0("'", utils::head(notGroups, 5), "'", collapse = ", "),
            ". Groups in the tree: ", showGroups(), call. = FALSE)
    }

    if (!is.character(skipGroups) || anyNA(skipGroups)) {
        stop("`skipGroups` must be a character vector of group names (use ",
            "character(0) to skip none)", call. = FALSE)
    }
    unknownSkip = setdiff(skipGroups, groups)
    if (length(unknownSkip) > 0) {
        stop("`skipGroups` has ", length(unknownSkip), " name", if (length(unknownSkip) > 1) "s",
            " that ", if (length(unknownSkip) > 1) "are" else "is", " not a group in ",
            "cellTypeTree: ", paste0("'", utils::head(unknownSkip, 5), "'", collapse = ", "),
            ". Groups in the tree: ", showGroups(), call. = FALSE)
    }
    if (length(setdiff(names(confList), skipGroups)) == 0) {
        stop("Every group in `confList` is in `skipGroups`, so there is nothing ",
            "left to analyse", call. = FALSE)
    }
    invisible(TRUE)
}

# checks that `spe` is a SummarizedExperiment (which includes SpatialExperiment
# and SingleCellExperiment)
checkSpe = function(spe) {
    if (!methods::is(spe, "SummarizedExperiment")) {
        stop("`spe` must be a SummarizedExperiment, for example a ",
            "SpatialExperiment", call. = FALSE)
    }
    invisible(spe)
}

# checks that `cellSize`, the size of the highlighted cells in a spatial plot,
# is a single positive number
checkCellSize = function(cellSize) {
    if (!is.numeric(cellSize) || length(cellSize) != 1 || is.na(cellSize) ||
        !is.finite(cellSize) || cellSize <= 0) {
        stop("`cellSize` must be a single positive number", call. = FALSE)
    }
    invisible(cellSize)
}

# checks that `gene`, the argument called `arg`, is a single string naming a row
# of the count matrix
checkGene = function(counts, gene, arg, assayName) {
    if (!is.character(gene) || length(gene) != 1 || is.na(gene)) {
        stop("`", arg, "` must be a single string naming a gene", call. = FALSE)
    }
    if (!gene %in% rownames(counts)) {
        stop("`", arg, " = \"", gene, "\"` is not a gene (a row) of the '",
            assayName, "' assay", call. = FALSE)
    }
    invisible(gene)
}

# ---- argument checks for the plot functions ----

# the groups of the cell type tree that have cells in the data. Without a tree
# the `group` column of the data is used.
groupsInData = function(cells, cellTypeTree = NULL) {
    if (is.null(cellTypeTree)) {
        return(unique(as.character(cells$group[!is.na(cells$group)])))
    }
    cellMap = buildCellTypeMap(cellTypeTree)$cellMap
    unique(cellMap$group[cellMap$cellType %in% cells$cellType & !is.na(cellMap$group)])
}

# checks that `imageID`, the image to plot, is a single image of the data
checkImage = function(cells, imageID) {
    images = unique(as.character(cells$imageID))
    if (!is.character(imageID) || length(imageID) != 1 || is.na(imageID)) {
        stop("`imageID` must be a single string", call. = FALSE)
    }
    if (!imageID %in% images) {
        stop("`imageID = \"", imageID, "\"` is not an image in the data. Images: ",
            paste(utils::head(images, 10), collapse = ", "),
            if (length(images) > 10) ", ...", call. = FALSE)
    }
    invisible(imageID)
}

# checks that every group the user passed exists in the data. `groups` is a
# named list: the names are the arguments, the values what the user passed.
checkGroupsInData = function(cells, cellTypeTree, groups) {
    available = groupsInData(cells, cellTypeTree)
    for (arg in names(groups)) {
        g = groups[[arg]]
        if (!is.character(g) || length(g) != 1 || is.na(g)) {
            stop("`", arg, "` must be a single group name", call. = FALSE)
        }
        if (!g %in% available) {
            stop("`", arg, " = \"", g, "\"` is not a group in the data. Groups: ",
                paste(utils::head(available, 10), collapse = ", "),
                if (length(available) > 10) ", ...", call. = FALSE)
        }
    }
    invisible(TRUE)
}

# checks that every label the user passed is a cell type or a group of the data
checkLabelsInData = function(cells, cellTypeTree, labels) {
    available = union(unique(as.character(cells$cellType)),
        groupsInData(cells, cellTypeTree))
    for (arg in names(labels)) {
        l = labels[[arg]]
        if (!is.character(l) || length(l) != 1 || is.na(l)) {
            stop("`", arg, "` must be a single cell type or group name", call. = FALSE)
        }
        if (!l %in% available) {
            stop("`", arg, " = \"", l, "\"` is neither a cell type nor a group in ",
                "the data", call. = FALSE)
        }
    }
    invisible(TRUE)
}

# checks that `subtype` (the argument called `arg`) has confidence scores for
# `group` in confList; the group itself stands for the whole group
checkSubtype = function(confList, group, subtype, arg) {
    if (!is.character(subtype) || length(subtype) != 1 || is.na(subtype)) {
        stop("`", arg, "` must be a single cell type name", call. = FALSE)
    }
    if (identical(group, subtype)) {
        return(invisible(TRUE))
    }
    conf = confList[[group]]
    if (!is.data.frame(conf)) {
        stop("`confList` has no confidence scores for the group '", group, "'",
            call. = FALSE)
    }
    available = setdiff(names(conf), "cellID")
    if (!subtype %in% available) {
        stop("`", arg, " = \"", subtype, "\"` is not a subtype of '", group,
            "' in `confList`. Subtypes: ", paste(available, collapse = ", "),
            call. = FALSE)
    }
    invisible(TRUE)
}

# checks that `niche` is one of the niches in the column `nicheCol`
checkNiche = function(cells, niche, nicheCol) {
    niches = sort(unique(as.character(cells[[nicheCol]][!is.na(cells[[nicheCol]])])))
    if (!is.character(niche) || length(niche) != 1 || is.na(niche)) {
        stop("`niche` must be a single niche name", call. = FALSE)
    }
    if (!niche %in% niches) {
        stop("`niche = \"", niche, "\"` is not a niche in `", nicheCol, "`. Niches: ",
            paste(niches, collapse = ", "), call. = FALSE)
    }
    invisible(niche)
}
