d = miniPipelineData()
quiet = function(x) suppressWarnings(suppressMessages(x))

# a confidence list that gives every cell confidence 1 in its own label
hardConf = function() {
    g = getCellData(d$spe)
    g = g[g$cellType %in% c("A", "B"), ]
    list(G = data.frame(cellID = g$cellID, A = as.numeric(g$cellType == "A"),
        B = as.numeric(g$cellType == "B")))
}

## ---- weights -------------------------------------------------------------
test_that("confidenceWeights gives subtypes their confidences and groups 1", {
    w = confidenceWeights(d$spe, d$tree, d$conf)
    expect_setequal(colnames(w), c("A", "B", "C", "G"))
    expect_equal(nrow(w), ncol(d$spe))
    inG = d$spe$cellType %in% c("A", "B")
    expect_equal(unname(rowSums(w[inG, c("A", "B")])), rep(1, sum(inG)))
    expect_true(all(w[inG, "G"] == 1) && all(w[!inG, "G"] == 0))
    expect_true(all(w[!inG, "C"] == 1))
    conf = d$conf$G
    i = match(conf$cellID, d$spe$cellID)
    expect_equal(unname(w[i, "A"]), conf$A)

    h = confidenceWeights(d$spe, d$tree, d$conf, weighted = FALSE)
    expect_true(all(h %in% c(0, 1)))
    expect_equal(unname(h[cbind(seq_len(nrow(h)), match(d$spe$cellType, colnames(h)))]), rep(1, nrow(h)))
    expect_equal(confidenceWeights(d$spe, d$tree, hardConf()), h)
})

test_that("confidenceWeights needs a confidence for every cell of a group", {
    conf = d$conf
    conf$G = conf$G[-1, ]
    expect_error(confidenceWeights(d$spe, d$tree, conf), "has no confidence for 1 of the group's cells")
})

## ---- outcome and subject checks --------------------------------------------
test_that("runColocalisationModels checks the outcome before computing anything", {
    spe = d$spe
    spe$three = c(a = "a", a = "a", b = "b", b = "b", c = "c", c = "c")[
        match(spe$imageID, paste0("img", 1:6))]
    spe$three = rep(c("a", "a", "b", "b", "c", "c"), each = 70)
    msgs = testthat::capture_messages(
        expect_error(runColocalisationModels(spe, d$tree, d$conf, "three"),
            "has 3 levels"))
    expect_false(any(grepl("Generating comparisons|Computing", msgs)))

    expect_error(runColocalisationModels(spe, d$tree, d$conf, "nope"),
        "is not a column of colData\\(spe\\)")

    spe$mixed = rep(c("a", "b"), 3 * 70)
    expect_error(runColocalisationModels(spe, d$tree, d$conf, "mixed"),
        "must have a single value per image")
})

test_that("images of one subject are combined; an image needs a single subject", {
    spe = d$spe
    spe$patient = paste0("p", rep(c(1, 1, 2, 3, 4, 4), each = 70))
    res = quiet(runColocalisationModels(spe, d$tree, d$conf, "outcome", subject = "patient"))
    expect_true(nrow(res$modelWeighted) > 0)

    spe$patient[2] = "other"
    expect_error(runColocalisationModels(spe, d$tree, d$conf, "outcome", subject = "patient"),
        "Every image must belong to a single subject")
    spe$patient[2] = NA
    expect_error(runColocalisationModels(spe, d$tree, d$conf, "outcome", subject = "patient"),
        "has missing values")
    expect_error(runColocalisationModels(spe, d$tree, d$conf, "outcome", subject = "nope"),
        "is not a column of colData\\(spe\\)")
})

test_that("runColocalisationModels says when subject was not given", {
    msgs = testthat::capture_messages(
        suppressWarnings(runColocalisationModels(d$spe, d$tree, d$conf, "outcome")))
    expect_equal(sum(grepl("`subject` was not given", msgs)), 1)
})

## ---- results -----------------------------------------------------------------
test_that("runColocalisationModels returns tests, per-image effects and the spicyR results", {
    res = quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome"))
    expect_named(res, c("weighted", "unweighted", "modelWeighted",
        "modelUnweighted", "relativeWeighted", "relativeUnweighted",
        "modelRelativeWeighted", "modelRelativeUnweighted", "spicy"))
    expect_true(all(c("from", "to", "estimate", "se", "p.value", "p.adj", "side") %in% names(res$modelWeighted)))
    expect_true(all(c("from", "to", "parent", "estimate", "p.value") %in% names(res$modelRelativeWeighted)))
    # no label is paired with itself or with a group that contains it
    pairs = paste(res$modelWeighted$from, res$modelWeighted$to)
    expect_false(any(c("A A", "A G", "G A", "B G") %in% pairs))
    # (a group with the only other cell type leaves no randomness under the null: not tested)
    expect_true(all(c("A B", "A C", "B C", "C A") %in% pairs))
    expect_true(all(c("imageID", "from", "to", "colocalisation") %in% names(res$weighted)))
})

test_that("confidence 1 in the original labels gives the unweighted results", {
    res = quiet(runColocalisationModels(d$spe, d$tree, hardConf(), "outcome"))
    expect_equal(res$modelWeighted, res$modelUnweighted)
    expect_equal(res$modelRelativeWeighted, res$modelRelativeUnweighted)
})

test_that("the unweighted test is spicyR's test of the cell types", {
    res = quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome"))
    cells = getCellData(d$spe)
    direct = quiet(spicyR::spicy(cells, condition = "outcome", r = 50, from = "A", to = c("B", "C")))
    tab = direct$cellResults
    m = res$modelUnweighted[res$modelUnweighted$from == "A" & res$modelUnweighted$to %in% c("B", "C"), ]
    m = m[match(tab$to, m$to), ]
    expect_equal(m$estimate, tab$excess_difference)
    expect_equal(m$p.value, tab$p_value)
})

## ---- other arguments and column names -----------------------------------------
test_that("runColocalisationModels requires a SummarizedExperiment", {
    for (bad in list(getCellData(d$spe), NULL, "spe", 1)) {
        expect_error(runColocalisationModels(bad, d$tree, d$conf, "outcome"),
            "`spe` must be a SummarizedExperiment, for example a SpatialExperiment")
    }
})

test_that("r and autoOverrideParents are checked before computing", {
    msgs = testthat::capture_messages(
        expect_error(runColocalisationModels(d$spe, d$tree, d$conf, "outcome", r = c(20, 50)),
            "`r` must be a single positive number \\(the radius\\)"))
    expect_false(any(grepl("Generating comparisons", msgs)))
    for (bad in list(0, -5, NA, "50", numeric(0), Inf)) {
        expect_error(runColocalisationModels(d$spe, d$tree, d$conf, "outcome", r = bad),
            "`r` must be")
    }
    for (bad in list(NA, "yes", c(TRUE, FALSE), 1)) {
        expect_error(runColocalisationModels(d$spe, d$tree, d$conf, "outcome",
            autoOverrideParents = bad), "`autoOverrideParents` must be TRUE or FALSE")
    }
})

test_that("runColocalisationModels takes custom column names and gives the same result", {
    default = quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome"))
    custom = quiet(runColocalisationModels(renamedSpe(d$spe), d$tree, d$conf, "outcome",
        cellType = "ct", cellID = "barcode", imageID = "sample", x = "px",
        y = "py"))
    expect_equal(custom, default)
})

test_that("the column arguments of runColocalisationModels are checked", {
    spe = renamedSpe(d$spe)
    expect_error(quiet(runColocalisationModels(spe, d$tree, d$conf, "outcome")),
        "`cellID = \"cellID\"` is not a column of the data. Available columns:")
    expect_error(runColocalisationModels(d$spe, d$tree, d$conf, "outcome", imageID = "nope"),
        "`imageID = \"nope\"` is not a column of the data")
    expect_error(runColocalisationModels(d$spe, d$tree, d$conf, "outcome", cellType = c("a", "b")),
        "`cellType` must be a single string naming a column")

    bad = d$spe
    bad$px = as.character(SpatialExperiment::spatialCoords(bad)[, 1])
    expect_error(runColocalisationModels(bad, d$tree, d$conf, "outcome", x = "px"),
        "`x = \"px\"` must be a numeric column without missing values")

    dup = d$spe
    dup$cellID[2] = dup$cellID[1]
    expect_error(runColocalisationModels(dup, d$tree, d$conf, "outcome"),
        "must contain unique, non-missing cell IDs")
})

test_that("relative results keep their columns when no pair has an eligible parent", {
    res = quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome",
        parentOverrides = c(A = "A", B = "B")))
    for (v in c("relativeWeighted", "relativeUnweighted")) {
        expect_equal(nrow(res[[v]]), 0)
        expect_named(res[[v]], c("imageID", "from", "to", "colocalisation", "parent"))
    }
    for (v in c("modelRelativeWeighted", "modelRelativeUnweighted")) {
        expect_equal(nrow(res[[v]]), 0)
        expect_named(res[[v]], c("from", "to", "parent", "estimate", "se", "statistic", "p.value", "p.adj"))
    }
    expect_error(plotMetricBoxplot(res, d$spe, "A", "B", "outcome", metric = "kontextual"),
        "is not in the results")
})
