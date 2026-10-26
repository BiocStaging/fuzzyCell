# Train the chosen classifier on a (already scaled/subset) feature matrix. 
# ntree only applies to classifier = "rf".
# A random forest is grown as 25 sub-forests on the workers of BPPARAM and merged
# with randomForest::combine(); the number of sub-forests does not depend on the
# number of workers, so the forest does not either.
trainClassifier = function(x, y, classifier = c("pcaLogistic", "logistic", "rf", "lda", "svm"), ntree = 500,
                           BPPARAM = BiocParallel::SerialParam(), nPC = 50) {
    classifier = match.arg(classifier)
    y = factor(y)
    if (classifier %in% c("logistic", "pcaLogistic")) {
        # pcaLogistic: the cells' first nPC principal components (of the training cells) in place of the genes
        proj = NULL
        k = min(nPC, ncol(x) - 1, nrow(x) - 1)
        # glmnet needs at least two predictors, so with fewer than three genes the genes are used as they are
        if (classifier == "pcaLogistic" && k >= 2) {
            pc = irlba::prcomp_irlba(x, n = k, center = TRUE, scale. = FALSE)
            proj = list(genes = colnames(x), center = pc$center, rotation = pc$rotation)
            x = pc$x
            colnames(x) = paste0("PC", seq_len(ncol(x)))
        }
        # ridge multinomial logistic regression, the penalty chosen by 3-fold
        # cross-validated deviance; kept as its coefficients
        cv = glmnet::cv.glmnet(x, y, family = "multinomial", alpha = 0, nfolds = 3, type.measure = "deviance")
        co = stats::coef(cv, s = "lambda.min")
        beta = do.call(cbind, lapply(co, function(b) as.numeric(b[-1])))
        colnames(beta) = names(co)
        return(list(a0 = vapply(co, function(b) b[1], 0), beta = beta, genes = colnames(x), proj = proj))
    }
    if (classifier == "rf") {
        sizes = diff(round(seq(0, ntree, length.out = min(25, ntree) + 1)))
        forests = BiocParallel::bplapply(sizes, function(nt) randomForest::randomForest(x, y, ntree = nt),
            BPPARAM = BPPARAM)
        return(do.call(randomForest::combine, forests))
    }
    switch(classifier,
        lda = MASS::lda(x, y),
        svm = e1071::svm(x, y, kernel = "radial")
    )
}

# Multinomial resampling of every cell's counts at its own total, in one pass:
# each of a cell's transcripts is drawn with replacement from its observed
# transcripts. `layout` is from resampleLayout(). The cells are resampled in
# blocks of about `blockSize` transcripts, which bounds the memory a resample
# needs; the draws, and so the counts, do not depend on the blocks.
resampleLayout = function(countMat, blockSize = 2^22) {
    # whole numbers stored with floating-point error would shift the transcript boundaries
    countMat = round(countMat)
    nz = which(countMat > 0, arr.ind = TRUE)
    nz = nz[order(nz[, 1], nz[, 2]), , drop = FALSE]
    n = countMat[nz]
    totals = as.numeric(rowSums(countMat))
    list(nz = nz, ends = cumsum(n), start = c(0, cumsum(totals))[seq_len(nrow(countMat))],
        totals = totals, blocks = split(seq_along(totals), cumsum(totals) %/% blockSize),
        dim = dim(countMat), dimnames = dimnames(countMat))
}
resampleCounts = function(layout) {
    hits = integer(nrow(layout$nz))
    for (cells in layout$blocks) {
        cell = rep.int(cells, layout$totals[cells])
        if (length(cell) == 0) next
        # a uniform position among the cell's transcripts, then the gene it falls in
        pos = layout$start[cell] + ceiling(stats::runif(length(cell)) * layout$totals[cell])
        entry = findInterval(pos - 1, layout$ends) + 1
        # a block's cells own a contiguous run of the nonzero entries
        lo = min(entry)
        hits[lo:max(entry)] = tabulate(entry - lo + 1L, max(entry) - lo + 1L)
    }
    out = matrix(0L, layout$dim[1], layout$dim[2], dimnames = layout$dimnames)
    out[layout$nz] = hits
    out
}

# Predict class labels, uniformly as a character vector regardless of classifier
# type. MASS::lda returns a list ($class among other things), e1071::svm and
# randomForest::randomForest return the predicted factor directly.
# logPrior: optional named log class priors added to the classifier's log posterior (its training classes are
# balanced, so its own prior is uniform).
# The linear predictors of a logistic fit (projecting onto its principal components first for pcaLogistic).
logisticEta = function(fit, newx) {
    if (!is.null(fit$proj)) {
        newx = sweep(newx[, fit$proj$genes, drop = FALSE], 2, fit$proj$center) %*% fit$proj$rotation
        colnames(newx) = fit$genes
    }
    sweep(newx[, fit$genes, drop = FALSE] %*% fit$beta, 2, fit$a0, "+")
}

predictClassifierClass = function(fit, newx, classifier = c("pcaLogistic", "logistic", "rf", "lda", "svm"), logPrior = NULL) {
    classifier = match.arg(classifier)
    if (classifier %in% c("logistic", "pcaLogistic")) {
        eta = logisticEta(fit, newx)
        if (!is.null(logPrior)) eta = sweep(eta, 2, logPrior[colnames(fit$beta)], "+")
        return(colnames(fit$beta)[max.col(eta, ties.method = "first")])
    }
    if (!is.null(logPrior)) {
        if (classifier == "svm") stop("class priors are not available for classifier = \"svm\"", call. = FALSE)
        p = if (classifier == "lda") predict(fit, newx)$posterior else predict(fit, newx, type = "prob")
        p = sweep(p, 2, exp(logPrior[colnames(p)]), "*")
        return(colnames(p)[max.col(p, ties.method = "first")])
    }
    switch(classifier,
        lda = as.character(predict(fit, newx)$class),
        svm = as.character(predict(fit, newx)),
        rf  = as.character(predict(fit, newx, type = "response"))
    )
}


#' Bootstrap classification confidence for one group
#'
#' The single-group worker behind bootstrapLineageConfidence().
#'  Trains a classifier on one group's cells, then repeatedly 
#' resamples each cell's counts (multinomial resampling) and
#' re-predicts its class. The fraction of resamples assigning a cell to each
#' class is its confidence for that class.
#'
#' @param countMat A cells x genes matrix of raw (non-negative integer) counts
#'   for the cells of one group. Normalised or log-transformed values are
#'   rejected, since the resampling needs each cell's total transcript count.
#' @param labels Character vector of the original cell type label for each
#'   row of `countMat`.
#' @param nBoot Number of bootstrap resamples.
#' @param nPerClass Minimum number of training cells per class; classes
#'   smaller than this are up-sampled with replacement.
#' @param BPPARAM A [BiocParallel::BiocParallelParam-class] object saying
#'   how to run the bootstrap resamples, for example
#'   [BiocParallel::MulticoreParam()]. The default runs serially. If `seed` is
#'   given it is set on `BPPARAM`, and the results are identical whichever
#'   backend and number of workers is used.
#' @param seed Random seed for reproducibility, or `NULL`. Recommended for reproducibility.
#' @param classifier One of `"pcaLogistic"` (ridge multinomial logistic
#'   regression, via glmnet, with the penalty chosen by cross-validation, on the
#'   first 50 principal components of the training cells, via irlba; the
#'   default), `"logistic"` (the same on all genes), `"rf"` (random forest, via randomForest), `"lda"` (linear
#'   discriminant analysis, via MASS) or `"svm"` (support vector machine, via
#'   e1071).
#' @param ntree Number of trees; only used when `classifier = "rf"`.
#' @param priors The classifier is trained with balanced classes, which gives
#'   every class the same prior. `"balanced"` keeps those equal priors;
#'   `"observed"` weights each class by its share of the cells, so a resample
#'   is assigned to its most probable class given how common the classes are.
#'
#' @return A numeric matrix with one row per cell (in the order of
#'   `countMat`) and one column per class in `labels`. Entry `[i, k]` is the
#'   fraction of bootstrap resamples in which cell `i` was classified as class
#'   `k`, so each row sums to 1.
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' cells = annotateCellTypeTree(getCellData(spe), heartCellTypeTree)
#' idx = which(cells$group %in% "T cells")
#' countMat = t(SummarizedExperiment::assay(spe, "counts")[, idx])
#' freq = bootstrapClassify(countMat, cells$cellType[idx], nBoot = 10, seed = 1)
#' head(freq)
#' }
#' @export
bootstrapClassify = function(countMat, labels, nBoot = 100, nPerClass = 50, BPPARAM = BiocParallel::SerialParam(),
                            seed = NULL, classifier = c("pcaLogistic", "logistic", "rf", "lda", "svm"), ntree = 500,
                            priors = c("balanced", "observed")) {
    classifier = match.arg(classifier)
    priors = match.arg(priors)
    labels = as.character(labels)
    # the logistic fits find their genes by name, so unnamed genes are named by position
    if (is.null(colnames(countMat))) colnames(countMat) = paste0("gene", seq_len(ncol(countMat)))

    if (nrow(countMat) != length(labels)) {
        stop("countMat has ", nrow(countMat), " rows but labels has ",
            length(labels), " elements; they must match (one row per cell)",
            call. = FALSE)
    }
    # the resampling treats each cell's total as a number of transcripts, so
    # normalised or log-transformed values would give meaningless (and
    # silently truncated) totals. a small tolerance allows whole numbers that
    # are stored as doubles with floating-point error.
    if (anyNA(countMat) || any(countMat < 0) ||
        any(abs(countMat - round(countMat)) > 1e-6)) {
        stop("countMat must contain raw non-negative integer counts, with no ",
            "missing values. Did you pass normalised or log-transformed data?",
            call. = FALSE)
    }

    # a class needs at least 2 cells: with one, the classifier has no
    # within-class variance to learn from
    classSizes = table(labels)
    tooSmall = names(classSizes)[classSizes < 2]
    if (length(tooSmall) > 0) {
        stop("Each class needs at least 2 cells to train the classifier, but ",
            paste0("'", tooSmall, "' has ", classSizes[tooSmall], collapse = ", "),
            ". Drop or merge these cell types.", call. = FALSE)
    }

    scaledMat = scaleCountMat(countMat)
    if (!is.null(seed)) set.seed(seed)

    # balance classes for training: use every real cell in the rarest class,
    # falling back to nPerClass (upsampled with replacement) only when even
    # the smallest class is smaller than that.
    targetSize = max(min(classSizes), nPerClass)
    balancedIdx = unlist(lapply(unique(labels), function(lab) {
        idx = which(labels == lab)
        # sample positions, not values: sample(k) alone would draw from 1:k
        idx[sample.int(length(idx), size = targetSize,
            replace = length(idx) < targetSize)]
    }))

    message("Minimum class size: ", targetSize)

    keep = apply(scaledMat, 2, var) > 0 # drop zero-variance genes
    if (!any(keep, na.rm = TRUE)) {
        stop("Every gene has zero variance across these cells, so there is ",
            "nothing to classify on. Check that countMat is a cells x genes ",
            "matrix of counts.", call. = FALSE)
    }
    keep = keep & !is.na(keep)
    # the forest's own random streams, apart from those of the resamples below
    if (!is.null(seed)) BiocParallel::bpRNGseed(BPPARAM) = seed + 7919L
    fit = trainClassifier(scaledMat[balancedIdx, keep, drop = FALSE], labels[balancedIdx], classifier, ntree = ntree,
        BPPARAM = BPPARAM)
    message("Finished training classifier.")


    # plain set.seed() before bplapply() does NOT make forked MulticoreParam
    # workers reproducible - the seed has to be set on the BPPARAM itself. Each
    # resample gets its own random stream, so results do not depend on the
    # number of workers.
    if (!is.null(seed)) BiocParallel::bpRNGseed(BPPARAM) = seed

    # observed priors: the classes' shares among these cells
    logPrior = if (priors == "observed") log(c(table(labels)) / length(labels))
    layout = resampleLayout(countMat)
    preds = BiocParallel::bplapply(seq_len(nBoot), function(b) {
        bootMat = resampleCounts(layout)
        predictClassifierClass(fit, scaleCountMat(bootMat)[, keep, drop = FALSE], classifier, logPrior)
    }, BPPARAM = BPPARAM)

    predMatrix = do.call(cbind, preds)

    allLabels = unique(labels)
    freqMatrix = vapply(
        allLabels, function(lab) rowMeans(predMatrix == lab),
        numeric(nrow(predMatrix))
    )
    colnames(freqMatrix) = allLabels
    freqMatrix
}

# Runs bootstrapClassify() for one lineage group in a fresh callr subprocess. A
# clean subprocess only carries the countMat/labels explicitly passed in as args, so
# the fork-based nBoot-level parallelism inside it forks from a small heap instead
# of whatever has accumulated in this session (the spe object, prior lineages'
# results) - that accumulated heap is what gets duplicated into every worker otherwise,
# ballooning RSS on large lineages.
bootstrapLineageInSubprocess = function(countMat, labels, nBoot, nPerClass, BPPARAM, seed, wd,
                                         classifier = "pcaLogistic", ntree = 500, priors = "balanced") {
    callr::r(
        func = function(countMat, labels, nBoot, nPerClass, BPPARAM, seed, classifier, ntree, priors) {
            fuzzyCell::bootstrapClassify(countMat, labels,
                nBoot = nBoot, nPerClass = nPerClass, BPPARAM = BPPARAM,
                seed = seed, classifier = classifier, ntree = ntree, priors = priors
            )
        },
        args = list(
            countMat = countMat, labels = labels, nBoot = nBoot,
            nPerClass = nPerClass, BPPARAM = BPPARAM, seed = seed,
            classifier = classifier, ntree = ntree, priors = priors
        ),
        wd = wd,
        stdout = "", stderr = "" # inherit this session's streams so progress bars still show live
    )
}

# The subprocess loads fuzzyCell from the library, not this session's copy, so
# checks that it is installed and is the same version as the one running here
# (e.g. under devtools::load_all() the installed copy may be older).
checkSubprocessInstall = function() {
    installed = tryCatch(utils::packageDescription("fuzzyCell", lib.loc = .libPaths())$Version,
        warning = function(w) NULL, error = function(e) NULL)
    if (is.null(installed)) {
        stop("`subprocess = TRUE` runs each group in a fresh R session, which needs fuzzyCell to be installed. ",
            "Install it, or use `subprocess = FALSE`.", call. = FALSE)
    }
    running = as.character(getNamespaceVersion("fuzzyCell"))
    if (!identical(installed, running)) {
        stop("`subprocess = TRUE` would run the installed fuzzyCell ", installed, ", but this session is running ",
            running, ". Install this version, or use `subprocess = FALSE`.", call. = FALSE)
    }
    invisible(TRUE)
}

#' Bootstrap confidence scores for every group in a cell type tree
#'
#' For each group in `cellTypeTree` that has children (see [cellTypeTree]),
#' runs [bootstrapClassify()] on that group's cells to get, for every cell, a
#' confidence score for each of the group's child cell types. The result is
#' the `confList` used by the downstream analyses, such as
#' [runProportionsPipeline()] and [runColocalisationModels()]. Cell types without a
#' group (ungrouped leaves) are not refined and keep hard labels.
#'
#' @inheritParams bootstrapClassify
#' @param spe A `SummarizedExperiment` (e.g. `SpatialExperiment`) with a counts
#'   assay and per-cell cell type labels and unique cell IDs in `colData`.
#' @param cellType Name of the `colData(spe)` column holding each cell's cell
#'   type label.
#' @param cellID Name of the `colData(spe)` column holding unique cell IDs.
#'   They are returned in the `cellID` column of each result.
#' @param assayName Name of the assay with the raw counts.
#' @param cellTypeTree A cell type tree, see [cellTypeTree].
#' @param skipGroups Names of groups to leave out of the result.
#' @param subprocess If `TRUE`, each group is bootstrapped in a fresh R session
#'   (via callr), which only holds that group's counts. With a forking backend
#'   such as [BiocParallel::MulticoreParam()], every worker starts as a copy of
#'   the session it is forked from, so this keeps the workers small on large
#'   data. It requires the same version of fuzzyCell to be installed.
#' @param seed Random seed, or `NULL` for no fixed seed. Each group is
#'   bootstrapped with `seed` plus its position in the list of groups, so every
#'   group is reproducible independently of the others.
#'
#' @return A named list with one element per group (named by group). Groups
#'   with no cells, or with only one cell type present, cannot be bootstrapped:
#'   they are skipped with a message and left out of the list. Each
#'   element is a `data.frame` with one row per cell in that group, a `cellID`
#'   column, and one column per child cell type of the group holding the
#'   cell's confidence score for that cell type (between 0 and 1, rows summing
#'   to 1 across the child columns).
#'   Each data frame also carries a `"settings"` attribute recording how it
#'   was computed (`nBoot`, `nPerClass`, `classifier`, `priors`, `ntree`, the
#'   group's `seed`, its `children` and the `fuzzyCellVersion`).
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 10,
#'     seed = 1)
#' names(confList)
#' head(confList[["T cells"]])
#' }
#' @export
bootstrapLineageConfidence = function(spe, cellTypeTree, nBoot = 100, nPerClass = 50,
                                       BPPARAM = BiocParallel::SerialParam(), seed = 51773,
                                       skipGroups = character(0), subprocess = FALSE,
                                       classifier = c("pcaLogistic", "logistic", "rf", "lda", "svm"), ntree = 500,
                                       priors = c("balanced", "observed"),
                                       cellType = "cellType", cellID = "cellID",
                                       assayName = "counts") {
    classifier = match.arg(classifier)
    priors = match.arg(priors)
    if (!is.logical(subprocess) || length(subprocess) != 1 || is.na(subprocess)) {
        stop("`subprocess` must be TRUE or FALSE", call. = FALSE)
    }
    if (subprocess) checkSubprocessInstall()

    checkSpe(spe)
    checkColumn(spe, cellType, "cellType")
    checkColumn(spe, cellID, "cellID")
    checkAssay(spe, assayName)
    ids = spe[[cellID]]
    if (anyNA(ids) || anyDuplicated(ids)) {
        stop("The `cellID` column '", cellID, "' must contain unique, non-missing ",
            "cell IDs", if (anyDuplicated(ids)) paste0(" (e.g. '", ids[anyDuplicated(ids)], "' is repeated)"),
            call. = FALSE)
    }

    spe = annotateCellTypeTree(spe, cellTypeTree, cellType = cellType)
    groups = setdiff(buildCellTypeMap(cellTypeTree)$groupsWithChildren, skipGroups)

    confList = lapply(seq_along(groups), function(gi) {
        g = groups[gi]

        idx = which(spe$group == g)
        labels = spe[[cellType]][idx]
        cellIds = spe[[cellID]][idx]

        # a group can only be bootstrapped if it has cells of at least 2 cell
        # types; otherwise leave it out (its cells keep their hard labels
        # downstream) and carry on with the other groups
        nClasses = length(unique(labels))
        if (length(idx) == 0 || nClasses < 2) {
            reason = if (length(idx) == 0) "no cells in the data" else
                sprintf("only 1 cell type present ('%s'), needs at least 2", labels[1])
            message(sprintf("[%d/%d] Lineage: %s - skipping: %s", gi, length(groups), g, reason))
            return(NULL)
        }

        # distinct-but-deterministic seed per lineage: makes each lineage
        # independently reproducible regardless of the other groups.
        # stays NULL (no fixed seed) when seed is NULL.
        groupSeed = if (is.null(seed)) NULL else seed + gi

        # how this group is bootstrapped, recorded with the result
        settings = list(
            nBoot = nBoot, nPerClass = nPerClass, classifier = classifier, priors = priors,
            ntree = ntree, seed = groupSeed, children = unique(labels),
            fuzzyCellVersion = tryCatch(as.character(utils::packageVersion("fuzzyCell")),
                error = function(e) NA_character_)
        )

        message(sprintf(
            "[%d/%d] Lineage: %s (%d cells, %d classes)",
            gi, length(groups), g, length(idx), length(unique(labels))
        ))
        startTime = Sys.time()

        countMat = denseCounts(assay(spe[, idx], assayName))

        freqMatrix = if (subprocess) {
            bootstrapLineageInSubprocess(countMat, labels, nBoot, nPerClass, BPPARAM, groupSeed,
                wd = getwd(), classifier = classifier, ntree = ntree, priors = priors
            )
        } else {
            bootstrapClassify(countMat, labels,
                nBoot = nBoot, nPerClass = nPerClass,
                BPPARAM = BPPARAM, seed = groupSeed, classifier = classifier, ntree = ntree, priors = priors
            )
        }

        message(sprintf(
            "[%d/%d] Lineage: %s done in %.1fs",
            gi, length(groups), g, as.numeric(Sys.time() - startTime, units = "secs")
        ))

        result = data.frame(cellID = cellIds, freqMatrix, check.names = FALSE)

        attr(result, "settings") = settings
        result
    })
    names(confList) = groups

    # drop groups that were skipped above
    confList = confList[!vapply(confList, is.null, logical(1))]
    if (length(confList) == 0) {
        stop("None of the groups could be bootstrapped: each had no cells or ",
            "only one cell type present", call. = FALSE)
    }
    confList
}

# A cells x genes dense count matrix from a genes x cells assay, dense or sparse
# (a sparse assay is made dense one group at a time).
denseCounts = function(m) {
    if (methods::is(m, "sparseMatrix")) m = Matrix::as.matrix(m)
    t(as.matrix(m))
}
