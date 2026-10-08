d = miniPipelineData()
quiet = function(x) suppressWarnings(suppressMessages(x))

test_that("clusterNiches labels every row with one of k niches", {
    set.seed(1)
    m = rbind(matrix(rnorm(60), 30, 2), matrix(rnorm(60, 5), 30, 2))
    rownames(m) = paste0("c", seq_len(nrow(m)))
    niches = stBoot:::clusterNiches(m, k = 2)
    expect_length(niches, 60)
    expect_setequal(unique(niches), c("niche1", "niche2"))
    expect_equal(names(niches), rownames(m))
})

## ---- runNicheIdentification -----------------------------------------------
# confidences that are exactly the hard labels, so the weighted and unweighted
# features are identical
hardConf = function() {
    cells = getCellData(d$spe)
    g = cells[cells$cellType %in% c("A", "B"), ]
    list(G = data.frame(cellID = g$cellID, A = as.numeric(g$cellType == "A"),
        B = as.numeric(g$cellType == "B")))
}

test_that("niche identification assigns every cell to one of k niches", {
    res = quiet(runNicheIdentification(d$spe, d$tree, d$conf, r = c(20, 50),
        k = 3))
    expect_length(res$nicheWeighted, ncol(d$spe))
    expect_lte(length(unique(res$nicheWeighted)), 3)
    expect_equal(nrow(res$lisaMatWeighted), ncol(d$spe))
})

test_that("results are reproducible given a seed", {
    run = function(seed) quiet(runNicheIdentification(d$spe, d$tree, d$conf,
        r = c(20, 50), k = 3, seed = seed))
    expect_identical(run(7), run(7))
    # the default seed is fixed too
    expect_identical(
        quiet(runNicheIdentification(d$spe, d$tree, d$conf, r = c(20, 50), k = 3)),
        quiet(runNicheIdentification(d$spe, d$tree, d$conf, r = c(20, 50), k = 3)))
})

test_that("identical features give identical unweighted and weighted niches", {
    res = quiet(runNicheIdentification(d$spe, d$tree, hardConf(),
        r = c(20, 50), k = 4, seed = 1))
    expect_equal(res$lisaMatUnweighted, res$lisaMatWeighted)
    # both clusterings start from the same seed, so the niches match exactly
    expect_identical(res$nicheUnweighted, res$nicheWeighted)
})

test_that("images too small for the largest radius are dropped with a message", {
    tiny = d$spe[, d$spe$imageID == "img1"]
    SpatialExperiment::spatialCoords(tiny) =
        SpatialExperiment::spatialCoords(tiny) / 3
    tiny$imageID = "tiny"
    colnames(tiny) = paste0("t_", colnames(tiny))
    tiny$cellID = colnames(tiny)
    spe = BiocGenerics::cbind(d$spe, tiny)
    # the tiny image is about 100 units wide, so 2 * 60 = 120 drops it
    msgs = testthat::capture_messages(
        res <- suppressWarnings(runNicheIdentification(spe, d$tree, d$conf,
            r = c(20, 60), k = 3)))
    expect_equal(sum(grepl("Dropping 1 of 7 images smaller than 2 x r = 120", msgs)), 1)
    expect_false(any(grepl("^t_", names(res$nicheWeighted))))
    expect_length(res$nicheWeighted, ncol(d$spe))
})

## ---- testNicheOutcome -----------------------------------------------------
# 12 images of 40 cells; niche1 is more common in the "case" images
nicheData = function() {
    set.seed(4)
    imgs = paste0("img", 1:12)
    outcome = rep(c("ctrl", "case"), each = 6)
    cells = do.call(rbind, lapply(seq_along(imgs), function(i) data.frame(
        cellID = paste0(imgs[i], "_", 1:40), imageID = imgs[i],
        outcome = outcome[i], subject = paste0("p", i))))
    pNiche1 = ifelse(outcome == "case", 0.7, 0.3)
    niche = setNames(unlist(lapply(seq_along(imgs), function(i) {
        ifelse(runif(40) < pNiche1[i], "niche1", "niche2")
    })), cells$cellID)
    list(cells = cells, niche = niche)
}

test_that("both tests report the second level minus the first", {
    nd = nicheData()
    run = function(test) quiet(testNicheOutcome(asSpe(nd$cells), nd$niche, nd$niche,
        "outcome", test = test))
    wilcox = run("wilcox")
    ttest = run("ttest")
    est = function(x) x$estimate[x$niche == "niche1"]

    # niche1 is more common in "case", the FIRST level alphabetically, so the
    # second level (ctrl) minus the first (case) is negative for both tests
    expect_lt(est(wilcox$testWeighted), 0)
    expect_lt(est(ttest$testWeighted), 0)

    # and the t-test estimate is exactly that difference in mean proportions
    props = dplyr::filter(ttest$propWeighted, niche == "niche1")
    byGroup = tapply(props$proportion,
        nd$cells$outcome[match(props$imageID, nd$cells$imageID)], mean)
    expect_equal(est(ttest$testWeighted), unname(byGroup["ctrl"] - byGroup["case"]))
})

test_that("testNicheOutcome requires a SummarizedExperiment", {
    nd = nicheData()
    for (bad in list(nd$cells, NULL, "spe", 1)) {
        expect_error(testNicheOutcome(bad, nd$niche, nd$niche, "outcome"),
            "`spe` must be a SummarizedExperiment, for example a SpatialExperiment")
    }
})

test_that("testNicheOutcome checks the outcome before computing", {
    nd = nicheData()
    three = nd$cells
    three$outcome = rep(c("a", "b", "c"), each = 4 * 40)
    expect_error(testNicheOutcome(asSpe(three), nd$niche, nd$niche, "outcome"),
        "has 3 levels \\(a, b, c\\)\\. The models currently compare two conditions only\\. Please subset the data to two conditions")

    mixed = nd$cells
    mixed$outcome = rep(c("a", "b"), 6 * 40)
    expect_error(testNicheOutcome(asSpe(mixed), nd$niche, nd$niche, "outcome"),
        "must have a single value per image")
    expect_error(testNicheOutcome(asSpe(mixed), nd$niche, nd$niche, "outcome"),
        "mixed models are not supported")

    expect_error(testNicheOutcome(asSpe(nd$cells), nd$niche, nd$niche, "nope"),
        "`outcomeCol = \"nope\"` is not a column of the data")
})

test_that("testNicheOutcome checks subjects", {
    nd = nicheData()
    shared = nd$cells
    shared$subject[shared$imageID == "img2"] = "p1"
    expect_error(testNicheOutcome(asSpe(shared), nd$niche, nd$niche, "outcome",
        subject = "subject"),
        "has more than one image for 1 subject: p1 \\(img1, img2\\)")
    expect_error(testNicheOutcome(asSpe(nd$cells), nd$niche, nd$niche, "outcome",
        subject = "patient"), "`subject = \"patient\"` is not a column")

    # one image per subject runs, with no 'not given' message
    msgs = testthat::capture_messages(res <- testNicheOutcome(asSpe(nd$cells),
        nd$niche, nd$niche, "outcome", subject = "subject"))
    expect_false(any(grepl("`subject` was not given", msgs)))
    expect_named(res, c("testUnweighted", "testWeighted", "propUnweighted",
        "propWeighted"))

    # without subject it says so, once
    msgs = testthat::capture_messages(testNicheOutcome(asSpe(nd$cells), nd$niche,
        nd$niche, "outcome"))
    expect_equal(sum(grepl("`subject` was not given", msgs)), 1)
})

test_that("niche identification and tests run on the real data (full data)", {
    skip_unless_full()
    f = fullTestData()
    res = suppressMessages(runNicheIdentification(f$spe, f$heartCellTypeTree,
        f$confList, r = c(20, 50), k = 4))
    expect_length(res$nicheWeighted, ncol(f$spe))
    expect_lte(length(unique(res$nicheWeighted)), 4)

    tst = suppressMessages(testNicheOutcome(f$spe,
        res$nicheUnweighted, res$nicheWeighted, "resolution_broad",
        subject = "subject"))
    expect_named(tst, c("testUnweighted", "testWeighted", "propUnweighted",
        "propWeighted"))
})

## ---- argument validation ----------------------------------------------------
run = function(...) {
    args = list(...)
    base = list(d$spe, d$tree, d$conf, r = c(20, 50), k = 3)
    base[names(args)] = args
    do.call(runNicheIdentification, base)
}

test_that("runNicheIdentification requires a SummarizedExperiment", {
    for (bad in list(getCellData(d$spe), NULL, "spe", 1)) {
        expect_error(runNicheIdentification(bad, d$tree, d$conf, r = c(20, 50), k = 3),
            "`spe` must be a SummarizedExperiment, for example a SpatialExperiment")
    }
    # the check comes before the others
    expect_error(runNicheIdentification(getCellData(d$spe), d$tree, d$conf),
        "`spe` must be a SummarizedExperiment")
})

test_that("k must be a whole number of at least 2 and below the cell count", {
    expect_error(runNicheIdentification(d$spe, d$tree, d$conf, r = c(20, 50)),
        "`k`, the number of niches, must be given")
    msg = "`k` \\(the number of niches\\) must be a single whole number of at least 2"
    for (bad in list(1, 0, -3, 2.5, NA, "3", c(2, 3), NULL)) {
        expect_error(run(k = bad), msg)
    }
    expect_error(run(k = ncol(d$spe)),
        "`k` = 420 must be smaller than the number of cells \\(420\\)")
})

test_that("r must be positive numbers", {
    msg = "`r` must be a vector of one or more positive numbers \\(radii\\)"
    for (bad in list(0, -5, c(20, NA), "50", numeric(0), Inf, NULL)) {
        expect_error(run(r = bad), msg)
    }
})

test_that("seed must be NULL or a single number", {
    msg = "`seed` must be NULL or a single number"
    expect_error(run(seed = "a"), msg)
    expect_error(run(seed = c(1, 2)), msg)
    expect_error(run(seed = NA_real_), msg)
    expect_no_error(quiet(run(seed = NULL)))
})

test_that("toLabels must be cell types or groups in the tree", {
    expect_error(run(toLabels = c("A", "Typo", "alsoWrong")),
        "`toLabels` has 2 labels that are not a cell type or group in cellTypeTree: 'Typo', 'alsoWrong'. Available labels: A, B, C, G")
    expect_error(run(toLabels = 5), "`toLabels` must be NULL or a character vector")
    expect_error(run(toLabels = character(0)), "`toLabels` must be NULL or")
    # cell types and groups are both fine, and give one column per label and radius
    res = quiet(run(toLabels = c("A", "G")))
    expect_equal(ncol(res$lisaMatWeighted), 2 * 2)
})

test_that("arguments are checked before any computation", {
    msgs = testthat::capture_messages(expect_error(run(k = 1)))
    expect_false(any(grepl("Building close-pairs", msgs)))
})

## ---- column arguments -------------------------------------------------------
test_that("custom column names give the same result as the defaults", {
    default = quiet(run(seed = 3))
    custom = quiet(runNicheIdentification(renamedSpe(d$spe), d$tree, d$conf,
        r = c(20, 50), k = 3, seed = 3, cellType = "ct", cellID = "barcode",
        imageID = "sample", x = "px", y = "py"))
    expect_equal(custom, default)
})

test_that("the column arguments are checked", {
    spe = renamedSpe(d$spe)
    # with the defaults the renamed columns are not found
    expect_error(quiet(runNicheIdentification(spe, d$tree, d$conf, r = c(20, 50),
        k = 3)), "`cellID = \"cellID\"` is not a column of the data")
    expect_error(run(cellType = "nope"),
        "`cellType = \"nope\"` is not a column of the data")
    expect_error(run(imageID = c("a", "b")),
        "`imageID` must be a single string naming a column of the data")

    # coordinates must be numeric and complete
    bad = d$spe
    bad$px = as.character(SpatialExperiment::spatialCoords(bad)[, 1])
    bad$py = SpatialExperiment::spatialCoords(bad)[, 2]
    expect_error(runNicheIdentification(bad, d$tree, d$conf, r = c(20, 50), k = 3,
        x = "px", y = "py"), "`x = \"px\"` must be a numeric column without missing values")
    bad$px = SpatialExperiment::spatialCoords(bad)[, 1]
    bad$px[3] = NA
    expect_error(runNicheIdentification(bad, d$tree, d$conf, r = c(20, 50), k = 3,
        x = "px", y = "py"), "must be a numeric column without missing values")

    # cell IDs must be unique
    dup = d$spe
    dup$cellID[2] = dup$cellID[1]
    expect_error(runNicheIdentification(dup, d$tree, d$conf, r = c(20, 50), k = 3),
        "must contain unique, non-missing cell IDs \\(e.g. 'i1_1' is repeated\\)")
})

test_that("testNicheOutcome takes custom column names and checks the niches", {
    nd = nicheData()
    names(nd$cells)[names(nd$cells) == "cellID"] = "barcode"
    names(nd$cells)[names(nd$cells) == "imageID"] = "sample"
    custom = quiet(testNicheOutcome(asSpe(nd$cells), nd$niche, nd$niche, "outcome",
        cellID = "barcode", imageID = "sample"))
    nd2 = nicheData()
    default = quiet(testNicheOutcome(asSpe(nd2$cells), nd2$niche, nd2$niche, "outcome"))
    expect_equal(custom$testWeighted, default$testWeighted)

    expect_error(testNicheOutcome(asSpe(nd$cells), nd$niche, nd$niche, "outcome"),
        "`cellID = \"cellID\"` is not a column of the data")

    nd2 = nicheData()
    expect_error(testNicheOutcome(asSpe(nd2$cells), unname(nd2$niche), nd2$niche, "outcome"),
        "`nicheUnweighted` must be a vector of niches named by cell ID")
    wrongIds = setNames(nd2$niche, paste0("x", seq_along(nd2$niche)))
    expect_error(testNicheOutcome(asSpe(nd2$cells), nd2$niche, wrongIds, "outcome"),
        "`nicheWeighted` must be a vector of niches named by cell ID")
})
