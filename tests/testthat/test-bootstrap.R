## small hand-made count matrix: two classes with different marker genes
makeCounts = function(nPerClass = 40, nGenes = 8) {
    set.seed(1)
    labels = rep(c("A", "B"), each = nPerClass)
    mu = matrix(1, 2 * nPerClass, nGenes)
    mu[labels == "A", 1:3] = 6
    mu[labels == "B", 4:6] = 6
    counts = matrix(stats::rpois(length(mu), mu), nrow = length(labels))
    colnames(counts) = paste0("g", seq_len(nGenes))
    list(counts = counts, labels = labels)
}

test_that("bootstrapClassify returns row-stochastic confidences", {
    d = makeCounts()
    freq = suppressMessages(suppressWarnings(
        bootstrapClassify(d$counts, d$labels, nBoot = 5, seed = 1)))
    expect_equal(dim(freq), c(80, 2))
    expect_equal(unname(rowSums(freq)), rep(1, 80))
    expect_setequal(colnames(freq), c("A", "B"))
})

test_that("bootstrapClassify is reproducible given a seed", {
    d = makeCounts()
    f1 = suppressMessages(suppressWarnings(
        bootstrapClassify(d$counts, d$labels, nBoot = 5, seed = 1)))
    f2 = suppressMessages(suppressWarnings(
        bootstrapClassify(d$counts, d$labels, nBoot = 5, seed = 1)))
    expect_equal(f1, f2)
})

test_that("bootstrapClassify recovers well-separated classes", {
    d = makeCounts()
    freq = suppressMessages(suppressWarnings(
        bootstrapClassify(d$counts, d$labels, nBoot = 10, seed = 1)))
    assigned = colnames(freq)[max.col(freq)]
    expect_gt(mean(assigned == d$labels), 0.9)
})

test_that("bootstrapClassify errors clearly on classes with fewer than 2 cells", {
    d = makeCounts()
    labels = d$labels
    labels[1] = "singleton"
    expect_error(
        suppressMessages(bootstrapClassify(d$counts, labels, nBoot = 2)),
        "'singleton' has 1")
})

test_that("bootstrapClassify rejects counts that are not raw integers", {
    d = makeCounts()
    expect_error(
        bootstrapClassify(log1p(d$counts), d$labels, nBoot = 2),
        "raw non-negative integer counts")
    withNA = d$counts
    withNA[1, 1] = NA
    expect_error(bootstrapClassify(withNA, d$labels, nBoot = 2),
        "raw non-negative integer counts")
    negative = d$counts
    negative[1, 1] = -1
    expect_error(bootstrapClassify(negative, d$labels, nBoot = 2),
        "raw non-negative integer counts")
    expect_error(bootstrapClassify(d$counts[-1, ], d$labels, nBoot = 2),
        "must match")
})

test_that("whole numbers stored as doubles are accepted", {
    d = makeCounts()
    shaky = d$counts * 0.1 / 0.1  # floating-point error, still whole numbers
    freq = suppressMessages(suppressWarnings(
        bootstrapClassify(shaky, d$labels, nBoot = 2, seed = 1)))
    expect_equal(dim(freq), c(80, 2))
})

test_that("bootstrapClassify errors when no gene varies", {
    d = makeCounts()
    expect_error(
        suppressMessages(bootstrapClassify(d$counts * 0, d$labels, nBoot = 2)),
        "zero variance")
})

test_that("bootstrapLineageConfidence records its settings", {
    d = makeCounts()
    se = SummarizedExperiment::SummarizedExperiment(
        list(counts = t(d$counts)),
        colData = S4Vectors::DataFrame(cellID = paste0("c", 1:80),
            cellType = d$labels))
    tree = list(Grp = list(children = c("A", "B")))
    conf = suppressMessages(suppressWarnings(
        bootstrapLineageConfidence(se, tree, nBoot = 3, nPerClass = 20,
            seed = 10)))
    settings = attr(conf$Grp, "settings")
    expect_equal(settings[c("nBoot", "nPerClass", "classifier", "priors",
        "ntree", "seed")], list(nBoot = 3, nPerClass = 20,
        classifier = "pcaLogistic", priors = "balanced", ntree = 500, seed = 11))
    expect_setequal(settings$children, c("A", "B"))
    expect_true("fuzzyCellVersion" %in% names(settings))
})

test_that("bootstrapLineageConfidence requires a SummarizedExperiment", {
    tree = list(Grp = list(children = c("A", "B")))
    for (bad in list(data.frame(cellID = "a", cellType = "A"), NULL, "spe", 1)) {
        expect_error(bootstrapLineageConfidence(bad, tree),
            "`spe` must be a SummarizedExperiment, for example a SpatialExperiment")
    }
})

test_that("bootstrapLineageConfidence validates its input", {
    d = makeCounts()
    mk = function(colData, assays = list(counts = t(d$counts))) {
        SummarizedExperiment::SummarizedExperiment(assays,
            colData = S4Vectors::DataFrame(colData))
    }
    tree = list(Grp = list(children = c("A", "B")))
    cd = data.frame(cellID = paste0("c", 1:80), cellType = d$labels)

    # columns / assay that do not exist
    expect_error(bootstrapLineageConfidence(mk(cd["cellType"]), tree),
        "`cellID = \"cellID\"` is not a column of colData\\(spe\\)")
    expect_error(bootstrapLineageConfidence(mk(cd["cellID"]), tree),
        "`cellType = \"cellType\"` is not a column of colData\\(spe\\)")
    expect_error(bootstrapLineageConfidence(
        mk(cd, assays = list(logcounts = t(d$counts))), tree),
        "`assayName = \"counts\"` is not an assay of spe. Available assays: logcounts")

    # the error lists the columns that are available
    expect_error(bootstrapLineageConfidence(mk(cd), tree, cellType = "ct"),
        "Available columns: cellID, cellType")

    # arguments of the wrong type
    expect_error(bootstrapLineageConfidence(mk(cd), tree, cellType = c("a", "b")),
        "`cellType` must be a single string")
    expect_error(bootstrapLineageConfidence(mk(cd), tree, cellID = NA_character_),
        "`cellID` must be a single string")
    expect_error(bootstrapLineageConfidence(mk(cd), tree, assayName = 1),
        "`assayName` must be a single string")

    # cell IDs must be unique
    dupIds = cd
    dupIds$cellID[2] = dupIds$cellID[1]
    expect_error(bootstrapLineageConfidence(mk(dupIds), tree),
        "unique, non-missing cell IDs \\(e.g. 'c1' is repeated\\)")
})

test_that("custom column and assay names give the same result as the defaults", {
    d = makeCounts()
    tree = list(Grp = list(children = c("A", "B")))
    se = SummarizedExperiment::SummarizedExperiment(list(counts = t(d$counts)),
        colData = S4Vectors::DataFrame(cellID = paste0("c", 1:80),
            cellType = d$labels))
    custom = SummarizedExperiment::SummarizedExperiment(
        list(raw = t(d$counts)),
        colData = S4Vectors::DataFrame(barcode = paste0("c", 1:80),
            ct = d$labels))
    run = function(...) suppressMessages(suppressWarnings(
        bootstrapLineageConfidence(..., nBoot = 3, seed = 1)))
    expect_equal(
        run(custom, tree, cellType = "ct", cellID = "barcode",
            assayName = "raw"),
        run(se, tree))
})

test_that("groups that cannot be bootstrapped are skipped with a message", {
    d = makeCounts()
    se = SummarizedExperiment::SummarizedExperiment(list(counts = t(d$counts)),
        colData = S4Vectors::DataFrame(cellID = paste0("c", 1:80),
            cellType = d$labels))
    # "Ghost" has no cells in the data
    tree = list(Grp = list(children = c("A", "B")),
        Ghost = list(children = c("X", "Y")))
    expect_message(
        conf <- suppressWarnings(bootstrapLineageConfidence(se, tree,
            nBoot = 2, seed = 1)),
        "skipping: no cells in the data")
    expect_named(conf, "Grp")

    # "Solo" only has cells of one of its cell types (A; there are no Q)
    solo = list(Solo = list(children = c("A", "Q")), Other = "B")
    expect_message(
        expect_error(bootstrapLineageConfidence(se, solo),
            "None of the groups could be bootstrapped"),
        "only 1 cell type present \\('A'\\)")
})

test_that("all classifiers run", {
    d = makeCounts()
    for (cl in c("lda", "svm", "rf")) {
        freq = suppressMessages(suppressWarnings(
            bootstrapClassify(d$counts, d$labels, nBoot = 3, seed = 1,
                classifier = cl, ntree = 20)))
        expect_equal(unname(rowSums(freq)), rep(1, 80), info = cl)
    }
})

test_that("bootstrapLineageConfidence returns one table per group (full data)", {
    skip_unless_full()
    d = fullTestData()
    expect_true(all(c("T cells", "Macrophage lineage") %in%
        names(d$confList)))
    expect_true(all(c("cellID", "CD8+ T cells") %in%
        names(d$confList[["T cells"]])))
})

test_that("crossfitLineageProbabilities gives probabilities for every cell of each group", {
    d = miniPipelineData()
    p = suppressMessages(crossfitLineageProbabilities(d$spe, d$tree, nFolds = 3, seed = 1))
    expect_named(p, "G")
    expect_setequal(names(p$G), c("cellID", "A", "B"))
    expect_equal(nrow(p$G), sum(d$spe$cellType %in% c("A", "B")))
    expect_equal(unname(rowSums(p$G[, c("A", "B")])), rep(1, nrow(p$G)))
    # A and B differ in genes 1-3 vs 4-6: most cells get their own label
    own = ifelse(d$spe$cellType[match(p$G$cellID, d$spe$cellID)] == "A", p$G$A, p$G$B)
    expect_gt(mean(own > 0.5), 0.7)
})

test_that("observed priors move fewer cells into a rare class", {
    set.seed(5)
    n = c(A = 400, B = 40)
    lab = rep(names(n), n)
    mu = rbind(A = c(3, 3, 1, 1), B = c(2, 2, 2, 2))
    counts = t(vapply(lab, function(l) rpois(4, mu[l, ] * 3), numeric(4)))
    colnames(counts) = paste0("g", 1:4)
    bal = suppressMessages(bootstrapClassify(counts, lab, nBoot = 10, seed = 1))
    obs = suppressMessages(bootstrapClassify(counts, lab, nBoot = 10, seed = 1, priors = "observed"))
    expect_lt(mean(obs[lab == "A", "B"]), mean(bal[lab == "A", "B"]))
    p1 = crossfitClassify(counts, lab, nFolds = 3, seed = 1)
    p2 = crossfitClassify(counts, lab, nFolds = 3, seed = 1, priors = "observed")
    expect_lt(sum(p2[, "B"]), sum(p1[, "B"]))
})

test_that("pcaLogistic classifies on principal components", {
    d = miniPipelineData()
    p = suppressMessages(crossfitLineageProbabilities(d$spe, d$tree, nFolds = 3, seed = 1, classifier = "pcaLogistic"))
    own = ifelse(d$spe$cellType[match(p$G$cellID, d$spe$cellID)] == "A", p$G$A, p$G$B)
    expect_gt(mean(own > 0.5), 0.7)
    b = suppressMessages(bootstrapLineageConfidence(d$spe, d$tree, nBoot = 3, seed = 1, classifier = "pcaLogistic"))
    expect_equal(unname(rowSums(b$G[, c("A", "B")])), rep(1, nrow(b$G)))
})

test_that("pcaLogistic fits the genes themselves when only two genes vary", {
    d = makeCounts()
    counts = d$counts[, c(1, 4)]
    freq = suppressMessages(bootstrapClassify(counts, d$labels, nBoot = 5, seed = 1))
    expect_equal(unname(rowSums(freq)), rep(1, nrow(freq)))
})

test_that("an unnamed count matrix gives the same confidences as a named one", {
    d = makeCounts()
    unnamed = unname(d$counts)
    for (classifier in c("pcaLogistic", "logistic")) {
        f1 = suppressMessages(bootstrapClassify(d$counts, d$labels, nBoot = 5, seed = 1, classifier = classifier))
        f2 = suppressMessages(bootstrapClassify(unnamed, d$labels, nBoot = 5, seed = 1, classifier = classifier))
        expect_equal(f2, f1)
    }
    p1 = crossfitClassify(d$counts, d$labels, nFolds = 3, seed = 1)
    p2 = crossfitClassify(unnamed, d$labels, nFolds = 3, seed = 1)
    expect_equal(p2, p1)
})

test_that("whole-number counts with floating-point error are resampled as their rounded values", {
    x = rbind(c(6 * 0.1 / 0.1, 1))
    set.seed(1)
    draws = replicate(500, resampleCounts(resampleLayout(x))[1, 2])
    expect_gt(mean(draws), 0)
    expect_equal(sum(resampleCounts(resampleLayout(x))), 7)
})

test_that("resampling in blocks of cells gives the same counts", {
    d = makeCounts()
    set.seed(1)
    whole = resampleCounts(resampleLayout(d$counts))
    set.seed(1)
    blocked = resampleCounts(resampleLayout(d$counts, blockSize = 30))
    expect_gt(length(resampleLayout(d$counts, blockSize = 30)$blocks), 1)
    expect_identical(blocked, whole)
})

test_that("crossfit folds keep every class in every training split", {
    d = makeCounts()
    lab = d$labels
    lab[lab == "B"][-(1:2)] = "A" # B has only 2 cells
    for (seed in 1:5) {
        p = crossfitClassify(d$counts, lab, nFolds = 2, seed = seed)
        expect_equal(unname(rowSums(p)), rep(1, nrow(p)))
        expect_true(all(p[, "B"] > 0))
    }
})

test_that("crossfitClassify rejects singleton classes and impossible fold counts", {
    d = makeCounts()
    lab = d$labels
    lab[lab == "B"][-1] = "A"
    expect_error(crossfitClassify(d$counts, lab, nFolds = 3), "'B' has 1")
    for (k in list(1, 2.5, NA, c(2, 3), nrow(d$counts) + 1)) {
        expect_error(crossfitClassify(d$counts, d$labels, nFolds = k), "`nFolds` must be")
    }
})
