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
#' features of each cell are its local indicators of spatial association from
#' lisaClust (`lisaClust::lisa()`) at several radii.
#'
#' @details
#' For every cell, every label in `toLabels` and every radius in `r`, the
#' feature is the local co-localisation of the cell's neighbours with that label
#' (lisaClust's local K or L function, centred so that 0 means no attraction or
#' repulsion). The weighted features count each neighbour with its bootstrap
#' confidence of having the label, instead of 0 or 1 (see [confidenceWeights()]),
#' and a label's density is its summed weight per unit area.
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
#' @param r Radii of the local curves: one or more positive numbers.
#' @param k Number of niches, a whole number of at least 2 and below the number
#'   of cells.
#' @param window,lisaFunc The image window and the local function (`"K"` or
#'   `"L"`), as in `lisaClust::lisa()`.
#' @param seed Random seed for the k-means clustering, as in
#'   [bootstrapLineageConfidence()], or `NULL` for no fixed seed. The same seed
#'   is used for the unweighted and weighted clusterings.
#'
#' @return A list with the feature matrices (`lisaMatUnweighted`,
#'   `lisaMatWeighted`) and named cell-to-niche vectors `nicheUnweighted` and
#'   `nicheWeighted`.
#'
#' @examples
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' niches = runNicheIdentification(lungSubset, lungCellTypeTree, lungConfList,
#'     r = 50, k = 6)
#' table(weighted = niches$nicheWeighted, unweighted = niches$nicheUnweighted)
#' @name nicheIdentification
#' @export
runNicheIdentification = function(spe, cellTypeTree, confList, toLabels = NULL, r = c(20, 50, 100),
                                   k, window = "convex", lisaFunc = "K", cores = 1, seed = 51773,
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
    if (k >= nrow(cells)) {
        stop("`k` = ", k, " must be smaller than the number of cells (",
            nrow(cells), ")", call. = FALSE)
    }

    lisaMat = list()
    for (v in c("unweighted", "weighted")) {
        message("Computing ", v, " LISA curves...")
        w = weightMatrix(cells, cellTypeTree, confList, v == "weighted")
        lisaMat[[v]] = lisaClust::lisa(cells[, c("cellID", "imageID", "cellType", "x", "y")], r = r,
            window = window, lisaFunc = lisaFunc, cores = cores,
            cellWeights = w[, intersect(colnames(w), toLabels), drop = FALSE])
    }

    # the same seed for both clusterings, so that they differ only through the
    # features and not through different random starts
    message("Clustering niches...")
    cluster = function(mat) {
        localSeed(seed)
        clusterNiches(mat, k = k)
    }
    nicheUnweighted = cluster(lisaMat$unweighted)
    nicheWeighted = cluster(lisaMat$weighted)

    list(
        lisaMatUnweighted = lisaMat$unweighted, lisaMatWeighted = lisaMat$weighted,
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
#' data(lungSubset)
#' data(lungCellTypeTree)
#' data(lungConfList)
#' niches = runNicheIdentification(lungSubset, lungCellTypeTree, lungConfList,
#'     r = 50, k = 6)
#' res = testNicheOutcome(lungSubset, niches$nicheUnweighted,
#'     niches$nicheWeighted, outcomeCol = "disease_status", subject = "subject")
#' res$testWeighted
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
