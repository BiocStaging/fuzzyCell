# Train the chosen classifier on a (already scaled/subset) feature matrix. 
# ntree only applies to classifier = "rf".
trainClassifier = function(x, y, classifier = c("lda", "svm", "rf"), ntree = 500) {
    classifier = match.arg(classifier)
    y = factor(y)
    switch(classifier,
        lda = MASS::lda(x, y),
        svm = e1071::svm(x, y, kernel = "radial"),
        rf  = randomForest::randomForest(x, y, ntree = ntree)
    )
}

# Predict class labels, uniformly as a character vector regardless of classifier
# type. MASS::lda returns a list ($class among other things), e1071::svm and
# randomForest::randomForest return the predicted factor directly.
predictClassifierClass = function(fit, newx, classifier = c("lda", "svm", "rf")) {
    classifier = match.arg(classifier)
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
#' @param classifier One of `"lda"` (linear discriminant analysis, via
#'   MASS), `"svm"` (support vector machine, via e1071) or `"rf"` (random
#'   forest, via randomForest).
#' @param ntree Number of trees; only used when `classifier = "rf"`.
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
                            seed = NULL, classifier = c("lda", "svm", "rf"), ntree = 500) {
    classifier = match.arg(classifier)
    labels = as.character(labels)

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
    fit = trainClassifier(scaledMat[balancedIdx, keep, drop = FALSE], labels[balancedIdx], classifier, ntree = ntree)
    message("Finished training classifier.")

    totals = rowSums(countMat)

    # plain set.seed() before bplapply() does NOT make forked MulticoreParam
    # workers reproducible - the seed has to be set on the BPPARAM itself. Each
    # resample gets its own random stream, so results do not depend on the
    # number of workers.
    if (!is.null(seed)) BiocParallel::bpRNGseed(BPPARAM) = seed

    preds = BiocParallel::bplapply(seq_len(nBoot), function(b) {
        bootMat = t(vapply(seq_len(nrow(countMat)), function(i) {
            if (totals[i] == 0) {
                return(integer(ncol(countMat)))
            }
            as.integer(stats::rmultinom(1, totals[i], countMat[i, ]))
        }, integer(ncol(countMat))))
        colnames(bootMat) = colnames(countMat)

        predictClassifierClass(fit, scaleCountMat(bootMat)[, keep, drop = FALSE], classifier)
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
# of  whatever has accumulated in this session (the spe object, prior lineages'
# results) - that accumulated heap is what gets duplicated into every worker otherwise,
# ballooning RSS on large lineages.
bootstrapLineageInSubprocess = function(countMat, labels, nBoot, nPerClass, BPPARAM, seed, wd,
                                         classifier = "lda", ntree = 500) {
    callr::r(
        func = function(countMat, labels, nBoot, nPerClass, BPPARAM, seed, classifier, ntree) {
            set.seed(seed)
            scBoot::bootstrapClassify(countMat, labels,
                nBoot = nBoot, nPerClass = nPerClass, BPPARAM = BPPARAM,
                seed = seed, classifier = classifier, ntree = ntree
            )
        },
        args = list(
            countMat = countMat, labels = labels, nBoot = nBoot,
            nPerClass = nPerClass, BPPARAM = BPPARAM, seed = seed,
            classifier = classifier, ntree = ntree
        ),
        wd = wd,
        stdout = "", stderr = "" # inherit this session's streams so progress bars still show live
    )
}

# Reassembles a confList purely from already-written checkpoint files, without
# recomputing anything or needing spe/counts data at all. skipGroups must match
# whatever was passed to the original bootstrapLineageConfidence() call, so the two lists
# of "groups that should have a checkpoint" line up.
#' Reload bootstrap confidence scores from checkpoint files
#'
#' Rebuilds the `confList` returned by [bootstrapLineageConfidence()] purely
#' from checkpoint files written by an earlier call with `checkpointDir` set,
#' without recomputing anything or needing the expression data.
#'
#' @param cellTypeTree The cell type tree used in the original call, see
#'   [cellTypeTree].
#' @param checkpointDir The directory the checkpoint files (one `.rds` file
#'   per group) were written to.
#' @param skipGroups Names of groups to leave out. Must match the
#'   `skipGroups` of the original call, so that the groups expected to have a
#'   checkpoint line up.
#'
#' @return A named list with one `data.frame` per group, in the same format as
#'   the result of [bootstrapLineageConfidence()]. An error is raised if a
#'   group's checkpoint file is missing.
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' dir = tempfile()
#' confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 10,
#'     seed = 1, checkpointDir = dir)
#' confList2 = loadCheckpointedConfList(heartCellTypeTree, dir)
#' identical(confList, confList2)
#' }
#' @export
loadCheckpointedConfList = function(cellTypeTree, checkpointDir, skipGroups = character(0)) {
    groups = setdiff(buildCellTypeMap(cellTypeTree)$groupsWithChildren, skipGroups)

    confList = lapply(groups, function(g) {
        checkpointFile = file.path(checkpointDir, paste0(make.names(g), ".rds"))
        if (!file.exists(checkpointFile)) stop("No checkpoint file for lineage: ", g, " (", checkpointFile, ")")
        readRDS(checkpointFile)
    })
    names(confList) = groups
    confList
}


# Describes how the settings stored in a checkpoint differ from those of the
# current call: returns a character vector with one entry per difference
# (empty if they match). `ntree` only matters for random forests.
checkpointSettingsDiff = function(stored, current) {
    fmt = function(x) if (is.null(x)) "NULL" else paste(x, collapse = ", ")
    same = function(a, b) isTRUE(all.equal(a, b))
    diffs = character(0)
    for (nm in c("nBoot", "nPerClass", "classifier", "ntree", "seed")) {
        if (nm == "ntree" && !"rf" %in% c(stored$classifier, current$classifier)) next
        if (!same(stored[[nm]], current[[nm]])) {
            diffs = c(diffs, sprintf("%s (checkpoint: %s, this call: %s)",
                nm, fmt(stored[[nm]]), fmt(current[[nm]])))
        }
    }
    if (!setequal(stored$children, current$children)) {
        diffs = c(diffs, sprintf("child cell types (checkpoint: %s; this call: %s)",
            fmt(stored$children), fmt(current$children)))
    }
    diffs
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
#' Use [loadCheckpointedConfList()] to reload the result from checkpoint files.
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
#' @param checkpointDir If given, each group's result is saved there as
#'   soon as it finishes, and re-running resumes from existing files. Each
#'   group then runs in a fresh `callr` subprocess, which requires scBoot to
#'   be installed.
#' @param skipGroups Names of groups to leave out of the result.
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
#'   was computed (`nBoot`, `nPerClass`, `classifier`, `ntree`, the group's
#'   `seed`, its `children` and the `scBootVersion`), which is also saved in
#'   checkpoint files. When resuming from a checkpoint, a message reports
#'   the stored classifier and any settings that differ from the current call.
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
                                       BPPARAM = BiocParallel::SerialParam(), seed = 51773, checkpointDir = NULL,
                                       skipGroups = character(0),
                                       classifier = c("lda", "svm", "rf"), ntree = 500,
                                       cellType = "cellType", cellID = "cellID",
                                       assayName = "counts") {
    classifier = match.arg(classifier)

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

    if (!is.null(checkpointDir)) dir.create(checkpointDir, recursive = TRUE, showWarnings = FALSE)

    confList = lapply(seq_along(groups), function(gi) {
        g = groups[gi]
        checkpointFile = if (!is.null(checkpointDir)) file.path(checkpointDir, paste0(make.names(g), ".rds"))

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
        # independently reproducible regardless of run/resume order.
        # stays NULL (no fixed seed) when seed is NULL.
        groupSeed = if (is.null(seed)) NULL else seed + gi

        # how this group is bootstrapped, recorded with the result so a
        # checkpoint can later be compared with the settings of a re-run
        settings = list(
            nBoot = nBoot, nPerClass = nPerClass, classifier = classifier,
            ntree = ntree, seed = groupSeed, children = unique(labels),
            scBootVersion = tryCatch(as.character(utils::packageVersion("scBoot")),
                error = function(e) NA_character_)
        )

        if (!is.null(checkpointFile) && file.exists(checkpointFile)) {
            saved = readRDS(checkpointFile)
            stored = attr(saved, "settings")
            prefix = sprintf("[%d/%d] Lineage: %s - already checkpointed, loading from disk", gi, length(groups), g)
            if (is.null(stored)) {
                message(prefix, " (no settings recorded, so the classifier and parameters cannot be verified)")
            } else {
                model = if (identical(stored$classifier, "rf")) {
                    sprintf("rf, ntree = %s", stored$ntree)
                } else {
                    stored$classifier
                }
                message(prefix, sprintf(" (classifier = %s, nBoot = %s)", model, stored$nBoot))
                diffs = checkpointSettingsDiff(stored, settings)
                if (length(diffs) > 0) {
                    message("  Note: the checkpoint for '", g, "' differs from this call: ",
                        paste(diffs, collapse = "; "), ". Delete ", basename(checkpointFile),
                        " to recompute it.")
                }
            }
            return(saved)
        }

        message(sprintf(
            "[%d/%d] Lineage: %s (%d cells, %d classes)",
            gi, length(groups), g, length(idx), length(unique(labels))
        ))
        startTime = Sys.time()

        countMat = assay(spe[, idx], assayName) |>
            t() |>
            as.matrix()

        freqMatrix = if (!is.null(checkpointDir)) {
            bootstrapLineageInSubprocess(countMat, labels, nBoot, nPerClass, BPPARAM, groupSeed,
                wd = getwd(), classifier = classifier, ntree = ntree
            )
        } else {
            bootstrapClassify(countMat, labels,
                nBoot = nBoot, nPerClass = nPerClass,
                BPPARAM = BPPARAM, seed = groupSeed, classifier = classifier, ntree = ntree
            )
        }

        message(sprintf(
            "[%d/%d] Lineage: %s done in %.1fs",
            gi, length(groups), g, as.numeric(Sys.time() - startTime, units = "secs")
        ))

        result = data.frame(cellID = cellIds, freqMatrix, check.names = FALSE)

        attr(result, "settings") = settings
        if (!is.null(checkpointFile)) saveRDS(result, checkpointFile)
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
