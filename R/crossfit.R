# Class probabilities from a trained classifier, one column per class.
predictClassifierProb = function(fit, newx, classifier = c("pcaLogistic", "logistic", "rf", "lda")) {
    classifier = match.arg(classifier)
    if (classifier %in% c("logistic", "pcaLogistic")) {
        eta = logisticEta(fit, newx)
        eta = eta - apply(eta, 1, max)
        p = exp(eta)
        return(p / rowSums(p))
    }
    switch(classifier,
        rf = predict(fit, newx, type = "prob"),
        lda = predict(fit, newx)$posterior
    )
}

# Cross-fitted class probabilities for the cells of one group: the cells are
# split into folds, and each fold is predicted by a classifier trained on the
# other folds with balanced classes (as in bootstrapClassify()), so no cell is
# scored by a classifier that saw it.
crossfitClassify = function(countMat, labels, nFolds = 5, nPerClass = 50, seed = NULL,
                            classifier = c("pcaLogistic", "logistic", "rf", "lda"), ntree = 500,
                            priors = c("balanced", "observed"), BPPARAM = BiocParallel::SerialParam()) {
    classifier = match.arg(classifier)
    priors = match.arg(priors)
    labels = as.character(labels)
    if (!is.numeric(nFolds) || length(nFolds) != 1 || is.na(nFolds) || nFolds != round(nFolds) ||
        nFolds < 2 || nFolds > length(labels)) {
        stop("`nFolds` must be a whole number between 2 and the number of cells (", length(labels), ")",
            call. = FALSE)
    }
    # every training split needs every class, so a class needs cells in at least 2 folds
    classSizes = table(labels)
    tooSmall = names(classSizes)[classSizes < 2]
    if (length(tooSmall) > 0) {
        stop("Each class needs at least 2 cells, so that every fold's training cells include it, but ",
            paste0("'", tooSmall, "' has ", classSizes[tooSmall], collapse = ", "),
            ". Drop or merge these cell types.", call. = FALSE)
    }
    if (is.null(colnames(countMat))) colnames(countMat) = paste0("gene", seq_len(ncol(countMat)))
    blasThreads = setBlasThreads(1L)
    on.exit(setBlasThreads(blasThreads), add = TRUE)
    if (!is.null(seed)) set.seed(seed)
    scaledMat = scaleCountMat(countMat)
    keep = apply(scaledMat, 2, var) > 0
    keep = keep & !is.na(keep)
    scaledMat = scaledMat[, keep, drop = FALSE]
    classes = sort(unique(labels))
    # stratified folds: the cells of each class, shuffled, are dealt round the folds in turn, so a class of
    # at least 2 cells has cells in at least 2 folds and the fold sizes differ by at most 1
    dealt = unlist(lapply(classes, function(lab) {
        idx = which(labels == lab)
        idx[sample.int(length(idx))]
    }))
    fold = integer(length(labels))
    fold[dealt] = rep(seq_len(nFolds), length.out = length(labels))
    out = matrix(NA_real_, length(labels), length(classes), dimnames = list(NULL, classes))
    # one fold per task; each draws its balanced training set from its own random stream
    seeds = sample.int(.Machine$integer.max, nFolds)
    fits = parLapply(seq_len(nFolds), function(k) {
        set.seed(seeds[k])
        train = which(fold != k)
        sizes = table(labels[train])
        targetSize = max(min(sizes), nPerClass)
        balancedIdx = unlist(lapply(names(sizes), function(lab) {
            idx = train[labels[train] == lab]
            idx[sample.int(length(idx), size = targetSize, replace = length(idx) < targetSize)]
        }))
        fit = trainClassifier(scaledMat[balancedIdx, , drop = FALSE], labels[balancedIdx], classifier, ntree = ntree)
        test = which(fold == k)
        list(test = test, p = predictClassifierProb(fit, scaledMat[test, , drop = FALSE], classifier))
    }, BPPARAM)
    for (f in fits) out[f$test, colnames(f$p)] = f$p
    out[is.na(out)] = 0
    if (priors == "observed") {
        out = sweep(out, 2, c(table(factor(labels, levels = classes))) / length(labels), "*")
        out = out / rowSums(out)
    }
    out
}

#' Cross-fitted classification probabilities for every group of a tree
#'
#' An alternative to [bootstrapLineageConfidence()]: for each group with
#' children, a classifier is trained on the group's cells with balanced classes
#' and gives each cell the probability of each child cell type. The cells are
#' split into `nFolds` folds and each fold is scored by a classifier trained on
#' the other folds, so no cell is scored by a classifier that saw it. Because
#' the classes are balanced in training, the probabilities give every child cell
#' type the same prior weight.
#'
#' @inheritParams bootstrapLineageConfidence
#' @param nFolds Number of folds. The folds are stratified by cell type, so
#'   every fold's classifier is trained on every child cell type; each child cell
#'   type therefore needs at least 2 cells.
#' @param classifier `"pcaLogistic"` (the default), `"logistic"`, `"rf"` or `"lda"`, as in
#'   [bootstrapClassify()].
#' @param priors `"balanced"` (equal priors, as the classifier is trained) or
#'   `"observed"` (each class weighted by its share of the group's cells).
#' @param cores Number of cores to fit the folds on in parallel, as in
#'   [bootstrapClassify()].
#'
#' @return A confidence list in the format of [bootstrapLineageConfidence()],
#'   with probabilities in place of bootstrap confidences.
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' data(heartCellTypeTree)
#' probList = crossfitLineageProbabilities(spe, heartCellTypeTree, seed = 1)
#' head(probList[["Vascular"]])
#' }
#' @export
crossfitLineageProbabilities = function(spe, cellTypeTree, nFolds = 5, nPerClass = 50, seed = 51773,
                                        classifier = c("pcaLogistic", "logistic", "rf", "lda"), ntree = 500,
                                        priors = c("balanced", "observed"), cores = 1,
                                        skipGroups = character(0), cellType = "cellType",
                                        cellID = "cellID", assayName = "counts") {
    classifier = match.arg(classifier)
    priors = match.arg(priors)
    checkCores(cores)
    BPPARAM = parallelParam(cores)
    checkSpe(spe)
    checkColumn(spe, cellType, "cellType")
    checkColumn(spe, cellID, "cellID")
    checkAssay(spe, assayName)
    spe = annotateCellTypeTree(spe, cellTypeTree, cellType = cellType)
    groups = setdiff(buildCellTypeMap(cellTypeTree)$groupsWithChildren, skipGroups)
    probList = lapply(seq_along(groups), function(gi) {
        idx = which(spe$group == groups[gi])
        labels = spe[[cellType]][idx]
        if (length(idx) == 0 || length(unique(labels)) < 2) return(NULL)
        countMat = denseCounts(assay(spe[, idx], assayName))
        p = crossfitClassify(countMat, labels, nFolds = nFolds, nPerClass = nPerClass,
            seed = if (is.null(seed)) NULL else seed + gi, classifier = classifier, ntree = ntree,
            priors = priors, BPPARAM = BPPARAM)
        message(sprintf("[%d/%d] Lineage: %s done", gi, length(groups), groups[gi]))
        data.frame(cellID = spe[[cellID]][idx], p, check.names = FALSE)
    })
    names(probList) = groups
    probList[!vapply(probList, is.null, logical(1))]
}
