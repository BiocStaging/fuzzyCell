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
    expect_equal(settings[c("nBoot", "nPerClass", "classifier", "ntree",
        "seed")], list(nBoot = 3, nPerClass = 20, classifier = "lda",
        ntree = 500, seed = 11))
    expect_setequal(settings$children, c("A", "B"))
    expect_true("scBootVersion" %in% names(settings))
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

test_that("checkpointSettingsDiff reports only real differences", {
    stored = list(nBoot = 20, nPerClass = 50, classifier = "rf", ntree = 100,
        seed = 5, children = c("A", "B"))
    expect_length(scBoot:::checkpointSettingsDiff(stored, stored), 0)
    # integer vs numeric is not a difference; child order does not matter
    same = modifyList(stored, list(nBoot = 20L, children = c("B", "A")))
    expect_length(scBoot:::checkpointSettingsDiff(stored, same), 0)
    diffs = scBoot:::checkpointSettingsDiff(stored,
        modifyList(stored, list(nBoot = 100, classifier = "lda")))
    expect_length(diffs, 2)
    expect_match(diffs[1], "nBoot \\(checkpoint: 20, this call: 100\\)")
    # ntree is ignored unless a random forest is involved
    lda1 = modifyList(stored, list(classifier = "lda", ntree = 100))
    lda2 = modifyList(lda1, list(ntree = 500))
    expect_length(scBoot:::checkpointSettingsDiff(lda1, lda2), 0)
    # NULL seeds
    expect_length(scBoot:::checkpointSettingsDiff(
        modifyList(stored, list(seed = NULL)),
        modifyList(stored, list(seed = NULL))), 0)
    expect_length(scBoot:::checkpointSettingsDiff(stored,
        modifyList(stored, list(seed = NULL))), 1)
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

test_that("checkpointing round-trips (full data)", {
    skip_unless_full()
    skip_if(nzchar(Sys.getenv("DEVTOOLS_LOAD")), "package not installed")
    d = fullTestData()
    dir = tempfile()
    conf = suppressMessages(suppressWarnings(
        bootstrapLineageConfidence(d$spe, d$heartCellTypeTree, nBoot = 5,
            seed = 1, checkpointDir = dir)))
    expect_equal(loadCheckpointedConfList(d$heartCellTypeTree, dir), conf)
})
