# Unweighted co-localisation of one pair of labels (the cross-type L-function,
# centred so that 0 means no attraction or repulsion). `counts` has one row per
# cell that has at least one neighbour within r, and one column per neighbouring cell
# type holding its edge-corrected neighbour count; fromLeaves / toLeaves are the
# cell types the from / to labels cover (one type for a leaf, several for a
# group). nFrom / nTo are the numbers of ALL cells of those types in the image,
# including cells with no neighbour within r. Such cells add nothing to the
# numerator, but they must be in the denominators: counting only cells that
# have a neighbour (as Statial does) biases the result towards attraction in
# sparse images.
calcLfunc = function(counts, fromLeaves, toLeaves, r, area, nFrom, nTo) {
    relevantCols = intersect(toLeaves, names(counts))

    # if no cell of a "to" type has a neighbour within r there is nothing to sum
    numerator = if (length(relevantCols) == 0) {
        0
    } else {
        sum(as.matrix(counts[cellTypeI %in% fromLeaves, ..relevantCols]), na.rm = TRUE)
    }

    Lfunction(numerator, r, area, nFrom, nTo)
}

# Per-image co-localisation (the cross-type L-function) of every pair in
# `pairs`, unweighted (original labels) and weighted (bootstrap confidences).
# Used by runColocalisationModels(), which documents the measures; not
# exported. `cells` is a per-cell data frame with imageID, cellID, cellType, x,
# y and group (ungrouped cell types have group equal to their cellType).
# Images whose smaller side is below 2 * r are dropped, with a message. Returns
# list(unweighted, weighted) of data.frames (imageID, from, to, colocalisation).
calcLfuncPerImage = function(cells, confList, pairs, r = 50, BPPARAM = BiocParallel::SerialParam(),
                              weightBy = c("both", "i", "j")) {
    weightBy = match.arg(weightBy)
    imageList = split(cells, cells$imageID)

    # skip images that are too small for the radius, and say which
    small = smallImages(cells, r)
    reportSmallImages(small, length(imageList), r)
    imageList = imageList[!names(imageList) %in% small$imageID]

    # maps a label (leaf or group) to the leaf types it covers - shared by both
    # the unweighted and weighted calculations below
    leafToGroup = setNames(cells$group[match(unique(cells$cellType), cells$cellType)], unique(cells$cellType)) # map leaf labels to group (named vector)
    leavesOf = function(label) { # return the leaves a group covers
        if (label %in% names(leafToGroup)) {
            label
        } # if it's a leaf itself, return the leaf label
        else {
            names(leafToGroup)[!is.na(leafToGroup) & leafToGroup == label]
        } # if it doesn't, pull in group labels
    }

    hardWeight = function(label) as.numeric(cells$cellType %in% leavesOf(label)) # assign hard weight (0/1)

    softWeight = function(label) {
        ownGroup = leafToGroup[label] # obtain group
        # use the bootstrap confidences only for a "genuine" child of a group:
        # the label has a group other than itself and weights were computed for
        # it
        if (!is.na(ownGroup) && ownGroup != label && label %in% colnames(confList[[ownGroup]])) {
            w = rep(0, nrow(cells)) # assign 0 weight to start
            idx = which(cells$group == ownGroup)
            conf = confList[[ownGroup]]
            w[idx] = conf[[label]][match(cells$cellID[idx], conf$cellID)] # assign weights from pre-computed
            w
        } else {
            hardWeight(label) # if it's not a genuine group, assign weight of 1
        }
    }
    buildWeight = function(label, useSoft) if (useSoft) softWeight(label) else hardWeight(label)

    # cache weight vectors by unique (label, useSoft) key - a group's own weight
    # only needs building once even though many pairs reference it
    weightCache = new.env(parent = emptyenv())
    getWeight = function(label, useSoft) {
        key = paste(label, useSoft)
        if (is.null(weightCache[[key]])) weightCache[[key]] = buildWeight(label, useSoft)
        weightCache[[key]]
    }

    # build every pair's weight vectors once, up front
    pairWeights = lapply(seq_len(nrow(pairs)), function(k) {
        list(
            weightI = getWeight(pairs$from[k], weightBy %in% c("both", "i")),
            weightJ = getWeight(pairs$to[k], weightBy %in% c("both", "j"))
        )
    })

    # one weight vector per unique `to` label
    toFirstOccurrence = match(unique(pairs$to), pairs$to)
    toWeightByLabel = setNames(
        lapply(toFirstOccurrence, function(k) pairWeights[[k]]$weightJ),
        pairs$to[toFirstOccurrence]
    )

    rowsByImage = split(seq_len(nrow(cells)), cells$imageID)


    perImage = BiocParallel::bplapply(names(imageList), function(imgID) {
        image = imageList[[imgID]]

        # build this image's closePairs geometry
        cellID = factor(seq_len(nrow(image)))
        ow = Statial::makeWindow(image, window = "square")
        imagePPP = spatstat.geom::ppp(x = image$x, y = image$y, window = ow)

        cp = spatstat.geom::closepairs(imagePPP, r, what = "indices", distinct = TRUE) |>
            data.frame()
        cp$i = factor(cp$i, levels = cellID)

        edgeDat = data.frame(
            i = factor(cellID, levels = cellID),
            edge = 1 / borderEdge(imagePPP, r)
        )
        cp = as.data.table(dplyr::left_join(as.data.table(cp), edgeDat, by = "i"))

        area = spatstat.geom::area(imagePPP)

        # unweighted counts table
        cp$cellTypeI = image$cellType[as.integer(as.character(cp$i))]
        cp$cellTypeJ = image$cellType[cp$j]

        counts = cp[, .(n = sum(edge, na.rm = TRUE)), by = .(i, cellTypeI, cellTypeJ)]
        counts = dcast(counts, i + cellTypeI ~ cellTypeJ, value.var = "n", fill = 0)

        typeCounts = table(image$cellType)
        imgRows = rowsByImage[[imgID]]

        neighborWeightSum = lapply(toWeightByLabel, function(wFull) {
            cp$wj = wFull[imgRows][cp$j] # store "to" cell weight for this image
            # aggregate the "to" weights over the neighbours of every i (from)
            # cell; a cell with no such neighbour within radius r gets 0
            agg = cp[, .(s = sum(edge * wj, na.rm = TRUE)), by = i]
            out = numeric(nrow(image))
            out[as.integer(as.character(agg$i))] = agg$s
            out
        })

        pairResults = lapply(seq_len(nrow(pairs)), function(k) {
            from = pairs$from[k]
            to = pairs$to[k]
            fromLeaves = leavesOf(from)
            toLeaves = leavesOf(to)

            # -- unweighted --
            nFrom = sum(typeCounts[names(typeCounts) %in% fromLeaves])
            nTo = sum(typeCounts[names(typeCounts) %in% toLeaves])
            unweightedL = if (nFrom == 0 || nTo == 0) {
                NA
            } else {
                calcLfunc(counts, fromLeaves, toLeaves, r, area, nFrom, nTo)
            }

            # -- weighted --
            weightI = pairWeights[[k]]$weightI[imgRows]
            weightJ = pairWeights[[k]]$weightJ[imgRows]
            numerator = sum(weightI * neighborWeightSum[[to]], na.rm = TRUE)
            sumFrom = sum(weightI, na.rm = TRUE)
            sumTo = sum(weightJ, na.rm = TRUE)
            weightedL = if (sumFrom == 0 || sumTo == 0) {
                NA
            } else {
                Lfunction(numerator, r, area, sumFrom, sumTo)
            }

            list(
                unweighted = data.frame(imageID = imgID, from = from, to = to, colocalisation = unweightedL),
                weighted = data.frame(imageID = imgID, from = from, to = to, colocalisation = weightedL)
            )
        })

        list(
            unweighted = dplyr::bind_rows(lapply(pairResults, `[[`, "unweighted")),
            weighted = dplyr::bind_rows(lapply(pairResults, `[[`, "weighted"))
        )
    }, BPPARAM = BPPARAM)

    list(
        unweighted = dplyr::bind_rows(lapply(perImage, `[[`, "unweighted")),
        weighted = dplyr::bind_rows(lapply(perImage, `[[`, "weighted"))
    )
}

# Relative co-localisation for every valid pair in one image, with hard or soft
# weights
# (fromWeightCache / toWeightCache). nParentCache and lambdaParentCache hold,
# for every unique parent, the number of parent cells in the image and the
# local parent density around each cell. They are computed once per image and
# shared by the unweighted and weighted calls, because thousands of pairs share
# only a dozen or so parents and recomputing the neighbour aggregation per pair
# is very slow. neighborSum is likewise built once per image and passed in.
kontextualPairResults = function(neighborSum, imgRows, nCellsTotal, validPairs,
                                  nParentCache, lambdaParentCache, fromWeightCache, toWeightCache, r) {
    pairResults = lapply(seq_len(nrow(validPairs)), function(k) {
        # for every from/to/parent
        # "parent" is only "to" parents for now
        # relative co-localisation is only defined for parents of the j (to)
        # cells
        from = validPairs$from[k]
        to = validPairs$to[k]
        parent = validPairs$parent[k]

        wFrom = fromWeightCache[[from]] # pull "from" membership
        wTo = toWeightCache[[to]] # pull "to" membership

        nParent = nParentCache[[parent]] # number of parent cells in image (precomputed once per parent)
        nChild2 = sum(wTo[imgRows], na.rm = TRUE) # number of j cells in image

        if (nParent == 0 || nChild2 == 0) {
            return(data.frame(imageID = NA, from = from, to = to, parent = parent, colocalisation = NA))
        }

        # local density of the parent population around every cell (precomputed
        # once per parent)
        lambdaParent = lambdaParentCache[[parent]]
        # expected density of child 2
        lambdaChild2 = lambdaParent * (nChild2 / nParent)

        wAdj = numeric(nCellsTotal)
        # is cell k a to cell? 0 or 1 / expected to-density around cell k
        wAdj[imgRows] = ifelse(lambdaChild2 > 0, wTo[imgRows] / lambdaChild2, 0)

        # aggregated densities of j cells for each focal cell i
        neighborAdjSum = neighborSum(wAdj)

        # "from" weights * i's own parent density * i's neighbours' adjusted
        # to-weight
        numerator = sum(wFrom[imgRows] * lambdaParent * neighborAdjSum, na.rm = TRUE)
        denominator = sum(wFrom[imgRows] * lambdaParent, na.rm = TRUE)

        if (denominator == 0) {
            return(data.frame(imageID = NA, from = from, to = to, parent = parent, colocalisation = NA))
        }

        data.frame(
            imageID = NA, from = from, to = to, parent = parent,
            colocalisation = sqrt(numerator / denominator / pi) - r
        )
    })
    dplyr::bind_rows(pairResults)
}

# Relative co-localisation (Kontextual): as calcLfuncPerImage(), but compared
# with a parent population. Only pairs whose `to` has a parent other than itself
# are returned (see parentOverrides). Returns list(unweighted, weighted) of
# data.frames (imageID, from, to, parent, colocalisation).
calcKontextualPerImage = function(cells, confList, pairs, cellTypeTree, r = 50, BPPARAM = BiocParallel::SerialParam(),
                                   weightBy = c("both", "i", "j"), parentOverrides = character(0)) {
    weightBy = match.arg(weightBy)
    imageList = split(cells, cells$imageID)

    # skip images that are too small for the radius, and say which
    small = smallImages(cells, r)
    reportSmallImages(small, length(imageList), r)
    imageList = imageList[!names(imageList) %in% small$imageID]

    leafToGroup = setNames(cells$group[match(unique(cells$cellType), cells$cellType)], unique(cells$cellType))
    leavesOf = function(label) {
        if (label %in% names(leafToGroup)) {
            label
        } else {
            names(leafToGroup)[!is.na(leafToGroup) & leafToGroup == label]
        }
    }
    hardWeight = function(label) as.numeric(cells$cellType %in% leavesOf(label))
    softWeight = function(label) {
        ownGroup = leafToGroup[label]
        if (!is.na(ownGroup) && ownGroup != label && label %in% colnames(confList[[ownGroup]])) {
            w = rep(0, nrow(cells))
            idx = which(cells$group == ownGroup)
            conf = confList[[ownGroup]]
            w[idx] = conf[[label]][match(cells$cellID[idx], conf$cellID)]
            w
        } else {
            hardWeight(label)
        }
    }

    # leavesOf() only reaches one level up (a leaf's immediate group). An
    # override parent like "Mesenchymal" sits several levels above any leaf, so
    # fall back to resolveTreeLeaves() (a full tree walk) whenever leavesOf()
    # comes up empty. This only happens for override parents; every ordinary
    # leaf or group still resolves through the fast leafToGroup path above.
    hardWeightAny = function(label) {
        direct = leavesOf(label)
        covers = if (length(direct) > 0) direct else resolveTreeLeaves(cellTypeTree, label)
        as.numeric(cells$cellType %in% covers)
    }

    validPairs = pairs |>
        dplyr::mutate(parent = dplyr::if_else(to %in% names(parentOverrides),
            unname(parentOverrides[to]),
            unname(leafToGroup[to])
        )) |>
        dplyr::filter(!is.na(parent), parent != to)

    emptyResult = data.frame(
        imageID = character(), from = character(), to = character(),
        parent = character(), colocalisation = numeric()
    )
    if (nrow(validPairs) == 0) {
        return(list(unweighted = emptyResult, weighted = emptyResult))
    }

    weightForFrom = if (weightBy %in% c("both", "i")) softWeight else hardWeight
    weightForTo = if (weightBy %in% c("both", "j")) softWeight else hardWeight

    # parent membership is always hard, and is shared between the unweighted and
    # weighted variants: which group a cell belongs to is not what is being
    # bootstrapped, only which leaf within the group it is
    uniqueParents = unique(validPairs$parent)
    parentWeightCache = setNames(lapply(uniqueParents, hardWeightAny), uniqueParents)
    fromWeightCacheHard = setNames(lapply(unique(validPairs$from), hardWeight), unique(validPairs$from))
    toWeightCacheHard = setNames(lapply(unique(validPairs$to), hardWeight), unique(validPairs$to))
    fromWeightCacheSoft = setNames(lapply(unique(validPairs$from), weightForFrom), unique(validPairs$from))
    toWeightCacheSoft = setNames(lapply(unique(validPairs$to), weightForTo), unique(validPairs$to))

    rowsByImage = split(seq_len(nrow(cells)), cells$imageID)
    nCellsTotal = nrow(cells)


    perImage = BiocParallel::bplapply(names(imageList), function(imgID) {
        image = imageList[[imgID]]

        # build this image's closePairs geometry 
        cellID = factor(seq_len(nrow(image)))
        ow = Statial::makeWindow(image, window = "square")
        imagePPP = spatstat.geom::ppp(x = image$x, y = image$y, window = ow)

        cp = spatstat.geom::closepairs(imagePPP, r, what = "indices", distinct = TRUE) |>
            data.frame()
        cp$i = factor(cp$i, levels = cellID)

        edgeDat = data.frame(
            i = factor(cellID, levels = cellID),
            edge = 1 / borderEdge(imagePPP, r)
        )
        cp = as.data.table(dplyr::left_join(as.data.table(cp), edgeDat, by = "i"))

        imgRows = rowsByImage[[imgID]]
        n = length(imgRows)

        # edge-corrected sum of `weightVec` over each cell's neighbours within
        # r, shared by every parent-density lookup below and by each pair's
        # wAdj aggregation
        neighborSum = function(weightVec) {
            cp$wj = weightVec[imgRows][cp$j]
            agg = cp[, .(s = sum(edge * wj, na.rm = TRUE)), by = i]
            out = numeric(n)
            out[as.integer(as.character(agg$i))] = agg$s
            out
        }

        # local parent density, computed ONCE per unique parent for this image
        # and shared by every pair that uses it and by the unweighted and
        # weighted calls below (recomputing it per pair is what made relative
        # co-localisation scale so badly on trees with many groups and cell
        # types)
        nParentCache = setNames(lapply(uniqueParents, function(p) sum(parentWeightCache[[p]][imgRows], na.rm = TRUE)), uniqueParents)
        lambdaParentCache = setNames(lapply(uniqueParents, function(p) neighborSum(parentWeightCache[[p]]) / (pi * r^2)), uniqueParents)

        unweighted = kontextualPairResults(
            neighborSum, imgRows, nCellsTotal, validPairs,
            nParentCache, lambdaParentCache, fromWeightCacheHard, toWeightCacheHard, r
        )
        weighted = kontextualPairResults(
            neighborSum, imgRows, nCellsTotal, validPairs,
            nParentCache, lambdaParentCache, fromWeightCacheSoft, toWeightCacheSoft, r
        )
        unweighted$imageID = imgID
        weighted$imageID = imgID

        list(unweighted = unweighted, weighted = weighted)
    }, BPPARAM = BPPARAM)

    list(
        unweighted = dplyr::bind_rows(lapply(perImage, `[[`, "unweighted")),
        weighted = dplyr::bind_rows(lapply(perImage, `[[`, "weighted"))
    )
}

# fits colocalisation ~ outcome separately for every distinct (from, to) pair;
# works the same whether the pairs are leaf x leaf, group x group or leaf x group
runLModels = function(perImage, cells, outcomeCol, subject = NULL) {
    checkOutcome(cells, outcomeCol, subject)
    outcomeDat = dplyr::distinct(cells, imageID, .data[[outcomeCol]])
    dat = merge(perImage, outcomeDat, by = "imageID")

    pairList = dplyr::distinct(dat, from, to)

    lapply(seq_len(nrow(pairList)), function(k) {
        from = pairList$from[k]
        to = pairList$to[k]
        pairDat = dat |>
            dplyr::filter(from == !!from, to == !!to, !is.na(.data[[outcomeCol]]))

        # skip pairs with too little data to fit a meaningful model
        if (length(unique(pairDat[[outcomeCol]])) < 2 || nrow(pairDat) < 5) {
            return(data.frame(from = from, to = to, estimate = NA, statistic = NA, p.value = NA))
        }

        fit = tryCatch(lm(as.formula(paste0("colocalisation ~ `", outcomeCol, "`")), data = pairDat),
            error = function(e) NULL
        )
        if (is.null(fit)) {
            return(data.frame(from = from, to = to, estimate = NA, statistic = NA, p.value = NA))
        }

        coef = summary(fit)$coefficients
        data.frame(
            from = from, to = to, estimate = coef[2, "Estimate"],
            statistic = coef[2, "t value"], p.value = coef[2, "Pr(>|t|)"]
        )
    }) |>
        dplyr::bind_rows() |>
        dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))
}

#' Co-localisation analysis with soft cell type weights
#'
#' Measures, in each image and for every pair of leaf cell types and groups,
#' how strongly the cells of one type are attracted to or avoid the cells of the
#' other (co-localisation), both with hard (original-label) membership and with
#' bootstrap-confidence (soft) membership, then fits a linear model of each
#' measure against an outcome. A second measure, relative co-localisation,
#' compares the co-localisation with a parent population.
#'
#' @details
#' Co-localisation measures, for a pair of labels in an image, how strongly the
#' `from` cells are attracted to (positive values) or avoid (negative values)
#' the `to` cells within radius `r`; 0 means no attraction or repulsion. It is
#' the cross-type L-function. The unweighted version uses the original cell
#' type labels; the weighted version uses the bootstrap confidences as soft
#' membership weights. Cells near the image border are edge corrected, as in
#' Statial. Cell counts include cells with no neighbour within `r`, which
#' Statial leaves out; this matters in sparse images, where leaving them out
#' biases the result towards attraction. As in Statial, the normalising
#' constant is the product of the two cell counts, so for pairs whose from and
#' to cells overlap (for example a cell type with itself) the result is biased
#' slightly downwards, by about `-r / (2n)` for `n` cells.
#'
#' Relative co-localisation measures co-localisation relative to a parent
#' population, so that attraction that only arises because both cell types
#' prefer the same region is not counted. It follows the Kontextual method of
#' Statial, and on images without isolated cells the unweighted results equal
#' those of Statial's `Kontextual()`. By default the parent of a `to` cell type
#' is its own group.
#'
#' Every pair of labels is compared in both directions: group with group, group
#' with cell type, and cell type with cell type, where the groups are those
#' bootstrapped in `confList`. Relative co-localisation only has results for
#' pairs whose `to` has a parent other than itself (by default its own group).
#' Pairs whose `to` is a group or an ungrouped cell type are therefore absent
#' from the `relative*` results and models, not `NA`, unless
#' `autoOverrideParents` supplies a broader parent for them (see
#' [autoParentOverrides()]). Parents given in `parentOverrides` take precedence
#' over automatically found ones for the same label.
#'
#' @param spe A `SummarizedExperiment` (e.g. `SpatialExperiment`) with a
#'   cell ID, cell type, image and coordinates for every cell (see
#'   [getCellData()]). Images whose smaller side is below `2 * r` are dropped,
#'   with a message listing them.
#' @param cellTypeTree A cell type tree, see [cellTypeTree].
#' @param confList Confidence list from [bootstrapLineageConfidence()].
#' @param outcomeCol Name of the `colData(spe)` column with the per-image
#'   outcome. It must have exactly two levels (subset your data to two
#'   conditions first) and a single value per image: each image is modelled as
#'   one independent observation, so repeated measures and mixed models are
#'   not supported. An error is raised otherwise, before any computation.
#' @param subject Name of the `colData(spe)` column identifying the subject
#'   (for example the patient) each image comes from. Several images from the
#'   same subject are not independent, so when `subject` is given the images
#'   and subjects must map one-to-one: an error is raised if an image has more
#'   than one subject or a subject has more than one image, and you should keep
#'   one image per subject. With `NULL` (the default) this cannot be checked,
#'   and a message says so.
#' @param r Radius within which cells count as neighbours, in the units of the
#'   coordinates; a single positive number.
#' @param BPPARAM A [BiocParallel::BiocParallelParam-class] object saying how
#'   to run the per-image computations in parallel, for example
#'   [BiocParallel::MulticoreParam()]. The default runs serially.
#' @param weightBy Which cell in each pair carries a soft weight: `"both"`,
#'   `"i"` (the from cell) or `"j"` (the to cell).
#' @param parentOverrides Named character vector mapping a `to` label to the
#'   ancestor node (a group, wrapper or lineage in `cellTypeTree`) that
#'   relative co-localisation should compare against instead of the label's
#'   own group, for
#'   example `c("Fibroblast lineage" = "Mesenchymal")`.
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
#' @return A list with `weighted`, `unweighted`, `relativeWeighted` and
#'   `relativeUnweighted` (per image and pair), and `modelWeighted`,
#'   `modelUnweighted`, `modelRelativeWeighted` and
#'   `modelRelativeUnweighted` (one linear model per pair). In the models,
#'   `estimate` is the difference in co-localisation for the second level of the
#'   outcome compared with the first (alphabetical order for text, factor order for
#'   factors). The per-image results have columns `imageID`, `from`, `to` and
#'   `colocalisation` (and `parent` for the relative version), whatever the
#'   columns of `spe` are called.
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
runColocalisationModels = function(spe, cellTypeTree, confList, outcomeCol, r = 50, 
                           BPPARAM = BiocParallel::SerialParam(), # nolint
                           weightBy = c("both", "i", "j"), parentOverrides = character(0),
                           autoOverrideParents = TRUE, subject = NULL,
                           cellType = "cellType", cellID = "cellID",
                           imageID = "imageID", x = "x", y = "y") {
    # check the arguments before doing any work
    checkSpe(spe)
    weightBy = match.arg(weightBy)
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
    cells = annotateCellTypeTree(cells, cellTypeTree)
    cells$group[is.na(cells$group)] = cells$cellType[is.na(cells$group)]

    # drop images too small for `r` (once, with a message), then check the
    # outcome before any of the slow computation starts
    cells = dropSmallImages(cells, r)
    checkOutcome(cells, outcomeCol, subject)
    messageIfNoSubject(subject)

    message("Generating comparisons...")
    comparisons = generateComparisons(cellTypeTree, confList)
    pairs = dplyr::bind_rows(comparisons$groupGroup, comparisons$groupLeaf, comparisons$leafLeaf)

    effectiveOverrides = if (autoOverrideParents) {
        auto = autoParentOverrides(cellTypeTree)
        c(auto[!names(auto) %in% names(parentOverrides)], parentOverrides)
    } else {
        parentOverrides
    }

    message("Computing unweighted + weighted co-localisation...")
    coloc = calcLfuncPerImage(cells, confList, pairs, r = r, BPPARAM = BPPARAM, weightBy = weightBy)

    message("Computing unweighted + weighted relative co-localisation...")
    relative = calcKontextualPerImage(cells, confList, pairs, cellTypeTree,
        r = r, BPPARAM = BPPARAM,
        weightBy = weightBy, parentOverrides = effectiveOverrides
    )

    message("Fitting unweighted models...")
    modelUnweighted = runLModels(coloc$unweighted, cells, outcomeCol, subject)
    modelRelativeUnweighted = runLModels(relative$unweighted, cells, outcomeCol, subject)

    message("Fitting weighted models...")
    modelWeighted = runLModels(coloc$weighted, cells, outcomeCol, subject)
    modelRelativeWeighted = runLModels(relative$weighted, cells, outcomeCol, subject)

    list(
        weighted = coloc$weighted, unweighted = coloc$unweighted,
        modelWeighted = modelWeighted, modelUnweighted = modelUnweighted,
        relativeWeighted = relative$weighted, relativeUnweighted = relative$unweighted,
        modelRelativeWeighted = modelRelativeWeighted,
        modelRelativeUnweighted = modelRelativeUnweighted
    )
}
