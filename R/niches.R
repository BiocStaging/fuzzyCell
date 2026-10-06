# Close pairs, their distances and the edge-correction weights for every image,
# up to the largest radius in r. Images too small for that radius are dropped.
makeClosePairsGeometryNiches = function(cells, r, BPPARAM = BiocParallel::SerialParam()) {
    r = sort(unique(r))
    imageList = split(cells, cells$imageID)

    # skip images that are too small for the largest radius, and say which
    small = smallImages(cells, max(r))
    reportSmallImages(small, length(imageList), max(r))
    imageList = imageList[!names(imageList) %in% small$imageID]

    BiocParallel::bplapply(imageList, function(image) {
        cellID = factor(seq_len(nrow(image)))
        ow = Statial::makeWindow(image, window = "square")
        imagePPP = spatstat.geom::ppp(x = image$x, y = image$y, window = ow)

        closePairs = spatstat.geom::closepairs(imagePPP, max(r), what = "ijd", distinct = TRUE) |>
            data.frame()
        closePairs$i = factor(closePairs$i, levels = cellID)
        closePairs = as.data.table(closePairs)

        # edge correction: needs to be rerun for every radius
        edge = do.call(cbind, lapply(r, function(rr) 1 / borderEdge(imagePPP, rr)))
        colnames(edge) = as.character(r)

        list(closePairs = closePairs, area = spatstat.geom::area(imagePPP), edge = edge)
    }, BPPARAM = BPPARAM)
}

# unweighted per-cell LISA curves: for every cell, the local co-localisation of its
# neighbours of each `toLabels` label at each radius, using hard cell types
calcUnweightedLisaCurves = function(cells, closePairsGeom, toLabels, r, BPPARAM = BiocParallel::SerialParam()) {
    cellsByImage = split(cells, cells$imageID)

    # maps a leaf to the group it's covered by
    leafToGroup = setNames(cells$group[match(unique(cells$cellType), cells$cellType)], unique(cells$cellType))

    # maps a label (leaf or group) to the leaf types it covers
    leavesOf = function(label) {
        if (label %in% names(leafToGroup)) {
            # a cell type that is a leaf covers just itself
            label
        } else {
            # a group label covers all of its leaves
            names(leafToGroup)[!is.na(leafToGroup) & leafToGroup == label]
        }
    }

    r = sort(unique(r))
    rChr = as.character(r)

    perImage = BiocParallel::bplapply(names(closePairsGeom), function(imgID) {
        # extract image info
        cp = closePairsGeom[[imgID]]$closePairs
        area = closePairsGeom[[imgID]]$area
        edge = closePairsGeom[[imgID]]$edge
        imgCells = cellsByImage[[imgID]]
        nCells = nrow(imgCells)

        cp$cellTypeJ = imgCells$cellType[cp$j]
        typeCounts = table(imgCells$cellType)

        lisaCols = lapply(rChr, function(rl) {
            # close pairs within radius rl, and how many cells of each cell
            # type are within that radius of every cell i
            cpR = cp[d <= as.numeric(rl)]
            neighborCounts = cpR[, .N, by = .(i, cellTypeJ)]

            vapply(toLabels, function(to) {
                # the cell types to count (one for a leaf, several for a group)
                toLeaves = leavesOf(to)
                # total number of those cells in the whole image
                nTo = sum(typeCounts[names(typeCounts) %in% toLeaves])

                # to-type entirely absent from this image - no meaningful local co-localisation
                if (nTo == 0) {
                    return(rep(NA_real_, nCells))
                }

                # neighbour counts of the cell types of interest, summed by i
                matching = neighborCounts[cellTypeJ %in% toLeaves, .(n = sum(N)), by = i]
                rawCount = rep(0, nCells)
                rawCount[as.integer(as.character(matching$i))] = matching$n

                numerator = rawCount * edge[, rl] # edge correction
                Lfunction(numerator, as.numeric(rl), area, sumWeightsChild1 = 1, sumWeightsChild2 = nTo)
            }, numeric(nCells))
        })

        lisaMat = do.call(cbind, lisaCols)
        colnames(lisaMat) = paste(rep(rChr, each = length(toLabels)), rep(toLabels, times = length(rChr)), sep = "_")
        rownames(lisaMat) = imgCells$cellID
        lisaMat
    }, BPPARAM = BPPARAM)

    do.call(rbind, perImage)
}

# weighted per-cell LISA curves: as above, but neighbours are counted with their
# bootstrap confidence of being the `to` label
calcWeightedLisaCurves = function(cells, confList, closePairsGeom, toLabels, r, BPPARAM = BiocParallel::SerialParam()) {
    # maps a leaf to the group it's covered by
    leafToGroup = setNames(cells$group[match(unique(cells$cellType), cells$cellType)], unique(cells$cellType))

    # maps a label (leaf or group) to the leaf types it covers
    leavesOf = function(label) {
        if (label %in% names(leafToGroup)) {
            label
        } else {
            names(leafToGroup)[!is.na(leafToGroup) & leafToGroup == label]
        }
    }

    # hard weight: 1 for every cell belonging to this label, leaf or group alike
    hardWeight = function(label) as.numeric(cells$cellType %in% leavesOf(label))

    # soft weight: the real bootstrap confidence, but only when the label is a
    # genuine child of a group with a matching confList column - otherwise falls
    # back to hard weights
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

    # store weights
    toWeight = setNames(lapply(toLabels, softWeight), toLabels)

    r = sort(unique(r))
    rChr = as.character(r)

    rowsByImage = split(seq_len(nrow(cells)), cells$imageID)
    cellsByImage = split(cells, cells$imageID)

    perImage = BiocParallel::bplapply(names(closePairsGeom), function(imgID) {
        # extract image info
        cp = closePairsGeom[[imgID]]$closePairs
        area = closePairsGeom[[imgID]]$area
        edge = closePairsGeom[[imgID]]$edge
        imgCells = cellsByImage[[imgID]]
        nCells = nrow(imgCells)
        imgRows = rowsByImage[[imgID]]

        lisaCols = lapply(rChr, function(rl) {
            cpR = cp[d <= as.numeric(rl)]

            vapply(toLabels, function(to) {
                # total "to" cell weight in image
                wFull = toWeight[[to]][imgRows]
                sumTo = sum(wFull, na.rm = TRUE)

                # to-type has zero total soft weight in this image - no
                # meaningful local co-localisation
                if (sumTo == 0) {
                    return(rep(NA_real_, nCells))
                }

                # weight of the j cell of every close pair
                cpR$wj = wFull[cpR$j]

                # sum the j weights for each i cell
                agg = cpR[, .(s = sum(wj, na.rm = TRUE)), by = i]

                weightedCount = rep(0, nCells)
                weightedCount[as.integer(as.character(agg$i))] = agg$s

                numerator = weightedCount * edge[, rl]
                Lfunction(numerator, as.numeric(rl), area, sumWeightsChild1 = 1, sumWeightsChild2 = sumTo)
            }, numeric(nCells))
        })

        lisaMat = do.call(cbind, lisaCols)
        colnames(lisaMat) = paste(rep(rChr, each = length(toLabels)), rep(toLabels, times = length(rChr)), sep = "_")
        rownames(lisaMat) = imgCells$cellID
        lisaMat
    }, BPPARAM = BPPARAM)

    do.call(rbind, perImage)
}

# k-means clusters on LISA curves; features that are NA (a `to` type absent from
# the image) are set to 0, i.e. no attraction or repulsion
clusterNiches = function(lisaMat, k, ...) {
    lisaMat[is.na(lisaMat)] = 0
    km = kmeans(lisaMat, centers = k, ...)
    setNames(paste0("niche", km$cluster), rownames(lisaMat))
}

#' Niche identification
#'
#' Defines niches by k-means clustering of per-cell feature vectors, once with
#' hard cell type membership and once with bootstrap-confidence weights. The
#' features of each cell are its local co-localisation curves at several radii.
#'
#' @details
#' For every cell, every label in `toLabels` and every radius in `r`, the
#' feature is the local co-localisation of the cell's neighbours with that label,
#' centred so that 0 means no attraction or repulsion. Where a label is absent
#' from an image the feature is set to 0. The weighted features count each
#' neighbour with its bootstrap confidence of having the label, instead of 0 or
#' 1.
#'
#' Images whose smaller side is below `2 * max(r)` are too small to hold a
#' neighbourhood and are dropped, with a message; their cells get no niche.
#'
#' Niche names (`niche1`, `niche2`, ...) are arbitrary. K-means starts from
#' random cluster centres, so the niches depend on the random seed. The
#' unweighted and weighted clusterings are both started from the same `seed`,
#' so that if the weights made no difference the two sets of niches would be
#' identical.
#'
#' @inheritParams runColocalisationModels
#' @param toLabels Labels (leaves and/or groups) to build features from;
#'   defaults to every leaf and group in the tree. An error is raised for a
#'   label that is not in the tree.
#' @param r Radii for the local co-localisation curves: one or more positive numbers.
#' @param k Number of niches, a whole number of at least 2 and below the number
#'   of cells.
#' @param seed Random seed for the k-means clustering, as in
#'   [bootstrapLineageConfidence()], or `NULL` for no fixed seed. The same seed
#'   is used for the unweighted and weighted clusterings.
#' @param cellType,cellID,imageID Names of the `colData(spe)` columns holding
#'   each cell's cell type label, a unique cell ID, and the image it is in.
#' @param x,y Names of the columns holding the cell coordinates. For a
#'   `SpatialExperiment` whose coordinates are in `spatialCoords()`, the
#'   defaults `"x"` and `"y"` refer to its first two coordinates.
#'
#' @return A list with the feature matrices (`lisaMatUnweighted`,
#'   `lisaMatWeighted`) and named cell-to-niche vectors `nicheUnweighted` and
#'   `nicheWeighted`.
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 10,
#'     seed = 1)
#' niches = runNicheIdentification(spe, heartCellTypeTree, confList,
#'     r = c(20, 50), k = 5)
#' table(niches$nicheWeighted)
#' }
#' @name nicheIdentification
#' @export
runNicheIdentification = function(spe, cellTypeTree, confList, toLabels = NULL, r = c(20, 50, 100),
                                   k, BPPARAM = BiocParallel::SerialParam(), seed = 51773,
                                   cellType = "cellType", cellID = "cellID",
                                   imageID = "imageID", x = "x", y = "y") {
    # check the arguments before doing any work
    checkSpe(spe)
    if (missing(k)) stop("`k`, the number of niches, must be given", call. = FALSE)
    checkNumberOfClusters(k)
    checkRadii(r)
    checkSeed(seed)
    if (is.null(toLabels)) {
        treeInfo = buildCellTypeMap(cellTypeTree)
        toLabels = union(treeInfo$cellMap$cellType, treeInfo$groupsWithChildren)
    } else {
        checkTreeLabels(toLabels, cellTypeTree)
    }

    cells = standardiseCells(getCellData(spe), list(cellID = cellID, cellType = cellType,
        imageID = imageID, x = x, y = y))
    cells = annotateCellTypeTree(cells, cellTypeTree)
    cells$group[is.na(cells$group)] = cells$cellType[is.na(cells$group)]

    # drop images too small for the largest radius (once, with a message)
    cells = dropSmallImages(cells, max(r))
    if (k >= nrow(cells)) {
        stop("`k` = ", k, " must be smaller than the number of cells (",
            nrow(cells), ")", call. = FALSE)
    }

    message("Building close-pairs geometry...")
    closePairsGeom = makeClosePairsGeometryNiches(cells, r = r, BPPARAM = BPPARAM)

    message("Computing unweighted LISA curves...")
    lisaMatUnweighted = calcUnweightedLisaCurves(cells, closePairsGeom, toLabels, r = r, BPPARAM = BPPARAM)

    message("Computing weighted LISA curves...")
    lisaMatWeighted = calcWeightedLisaCurves(cells, confList, closePairsGeom, toLabels, r = r, BPPARAM = BPPARAM)

    # the same seed for both clusterings, so that they differ only through the
    # features and not through different random starts
    message("Clustering niches...")
    if (!is.null(seed)) set.seed(seed)
    nicheUnweighted = clusterNiches(lisaMatUnweighted, k = k)
    if (!is.null(seed)) set.seed(seed)
    nicheWeighted = clusterNiches(lisaMatWeighted, k = k)

    list(
        lisaMatUnweighted = lisaMatUnweighted, lisaMatWeighted = lisaMatWeighted,
        nicheUnweighted = nicheUnweighted, nicheWeighted = nicheWeighted
    )
}

#' Test niche proportions against an outcome
#'
#' Computes the proportion of each image's cells in each niche and tests every
#' niche against a two-level outcome with a Wilcoxon or t-test.
#'
#' @details
#' The `estimate` is the difference in the niche proportion for the second
#' level of the outcome compared with the first (alphabetical order for text,
#' factor order for factors), for both tests, as in [runColocalisationModels()]. For
#' `test = "wilcox"` this is the Hodges-Lehmann estimate of the shift in
#' location. Each image is one observation, so every image must come from a
#' different subject if `subject` is given.
#'
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`) with a cell
#'   ID, an image and the outcome for every cell in `colData`.
#' @param cellID,imageID Names of the `colData(spe)` columns holding a unique
#'   cell ID (matching the names of the niche vectors) and the image each cell
#'   is in.
#' @param nicheUnweighted,nicheWeighted Named cell-to-niche vectors, as
#'   returned by [runNicheIdentification()].
#' @param outcomeCol Name of the per-image outcome column in `colData(spe)`. It
#'   must have exactly two levels (subset your data to two conditions first) and
#'   a single value per image: repeated measures and mixed models are not
#'   supported. An error is raised otherwise, before any computation.
#' @param test `"wilcox"` or `"ttest"`.
#' @param subject Name of the `colData(spe)` column identifying the subject (for
#'   example the patient) each image comes from, or `NULL`. If given, every
#'   image must belong to one subject and every subject must contribute only
#'   one image, as in [runColocalisationModels()]. With `NULL` this cannot be checked,
#'   and a message says so.
#'
#' @return A list with `testUnweighted`, `testWeighted` (one row per niche)
#'   and `propUnweighted`, `propWeighted` (per image niche proportions, to
#'   pass to [nicheProportionBoxplot()]).
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 10,
#'     seed = 1)
#' niches = runNicheIdentification(spe, heartCellTypeTree, confList,
#'     r = c(20, 50), k = 5)
#' res = testNicheOutcome(spe, niches$nicheUnweighted,
#'     niches$nicheWeighted, outcomeCol = "resolution_broad",
#'     subject = "subject")
#' res$testWeighted
#' }
#' @export
testNicheOutcome = function(spe, nicheUnweighted, nicheWeighted, outcomeCol,
                             test = c("wilcox", "ttest"), subject = NULL,
                             cellID = "cellID", imageID = "imageID") {
    test = match.arg(test)

    # check the arguments before doing any work
    checkSpe(spe)
    cells = standardiseCells(getCellData(spe), list(cellID = cellID, imageID = imageID))
    for (niche in c("nicheUnweighted", "nicheWeighted")) {
        vec = get(niche)
        if (is.null(names(vec)) || !any(cells$cellID %in% names(vec))) {
            stop("`", niche, "` must be a vector of niches named by cell ID, with ",
                "at least some of the cell IDs in `cells`", call. = FALSE)
        }
    }
    checkOutcome(cells, outcomeCol, subject)
    messageIfNoSubject(subject)

    cells$nicheUnweighted = nicheUnweighted[match(cells$cellID, names(nicheUnweighted))]
    cells$nicheWeighted = nicheWeighted[match(cells$cellID, names(nicheWeighted))]

    # per-image proportion of each niche, long format (imageID, niche, proportion)
    nicheProportions = function(cells, niche) {
        props = table(sample = cells$imageID, cluster = cells[[niche]])
        props = props / rowSums(props)
        data.frame(
            imageID = rep(rownames(props), times = ncol(props)),
            niche = rep(colnames(props), each = nrow(props)),
            proportion = as.vector(props), stringsAsFactors = FALSE
        )
    }

    runTest = function(propDat, outcomeDat, outcomeCol) {
        dat = dplyr::inner_join(propDat, outcomeDat, by = "imageID")
        formula = as.formula(paste0("proportion ~ `", outcomeCol, "`"))

        lapply(unique(dat$niche), function(nicheLabel) {
            sub = dat |> dplyr::filter(niche == nicheLabel, !is.na(proportion), !is.na(.data[[outcomeCol]]))

            if (length(unique(sub[[outcomeCol]])) != 2 || nrow(sub) < 5) {
                return(data.frame(niche = nicheLabel, estimate = NA, statistic = NA, p.value = NA))
            }

            fit = tryCatch(
                {
                    if (test == "wilcox") {
                        stats::wilcox.test(formula, data = sub, conf.int = TRUE)
                    } else {
                        stats::t.test(formula, data = sub)
                    }
                },
                error = function(e) NULL
            )
            if (is.null(fit)) {
                return(data.frame(niche = nicheLabel, estimate = NA, statistic = NA, p.value = NA))
            }

            # both tests report the second level minus the first. t.test() gives
            # the two group means, so the difference is taken here;
            # wilcox.test() estimates first minus second, so its sign is flipped
            estimate = if (test == "wilcox") -unname(fit$estimate) else unname(diff(fit$estimate))

            data.frame(
                niche = nicheLabel, estimate = estimate,
                statistic = unname(fit$statistic), p.value = fit$p.value
            )
        }) |>
            dplyr::bind_rows() |>
            dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) |>
            dplyr::arrange(p.value)
    }

    outcomeDat = dplyr::distinct(cells, imageID, .data[[outcomeCol]])

    propUnweighted = nicheProportions(cells, "nicheUnweighted")
    propWeighted = nicheProportions(cells, "nicheWeighted")
    # recorded so that nicheProportionBoxplot() can say which table it was given
    attr(propUnweighted, "nicheType") = "unweighted"
    attr(propWeighted, "nicheType") = "weighted"

    testUnweighted = runTest(propUnweighted, outcomeDat, outcomeCol)
    testWeighted = runTest(propWeighted, outcomeDat, outcomeCol)

    list(
        testUnweighted = testUnweighted, testWeighted = testWeighted,
        propUnweighted = propUnweighted, propWeighted = propWeighted
    )
}
