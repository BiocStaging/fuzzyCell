d = miniPipelineData()
quiet = function(x) suppressWarnings(suppressMessages(x))

test_that("unweighted proportions sum to one per image", {
    cells = data.frame(imageID = rep(c("a", "b"), c(5, 3)),
        cellType = c("T", "T", "B", "B", "B", "T", "B", "B"))
    props = stBoot:::calcUnweightedProportions(cells)
    expect_equal(unname(rowSums(props)), c(1, 1))
    expect_equal(props["a", "T"], 2 / 5, ignore_attr = TRUE)
})

## ---- runProportionModels ------------------------------------------------------
imgs = paste0("img", 1:9)
props = matrix(seq(0.1, 0.9, length.out = 18), 9, 2,
    dimnames = list(imgs, c("T", "B")))
perImage = function(outcome, subjects = imgs) {
    data.frame(imageID = rep(imgs, each = 2), outcome = rep(outcome, each = 2),
        subject = rep(subjects, each = 2))
}

test_that("runProportionModels skips cell types with too little data", {
    p = matrix(c(0.2, 0.3, 0.4, 0.1), 2, 2, dimnames = list(c("a", "b"), c("T", "B")))
    cells = data.frame(imageID = c("a", "b"), outcome = c("x", "y"))
    res = stBoot:::runProportionModels(p, cells, "outcome")
    expect_equal(nrow(res), 2)
    expect_true(all(is.na(res$p.value)))
})

test_that("runProportionModels estimates second level minus first", {
    # images 1-4 are "ctrl", 5-9 are "case"; alphabetically "case" is the first
    # level, so the estimate is mean(ctrl) - mean(case)
    cells = perImage(rep(c("ctrl", "case"), c(4, 5)))
    res = stBoot:::runProportionModels(props, cells, "outcome")
    expect_equal(res$estimate[res$cellType == "T"],
        mean(props[1:4, "T"]) - mean(props[5:9, "T"]))
})

test_that("runProportionModels checks the outcome and subjects", {
    expect_error(stBoot:::runProportionModels(props, perImage(rep(c("a", "b", "c"), each = 3)),
        "outcome"), "has 3 levels \\(a, b, c\\)\\. The models currently compare two conditions only")
    mixed = data.frame(imageID = rep(imgs, each = 2), outcome = rep(c("a", "b"), 9))
    expect_error(stBoot:::runProportionModels(props, mixed, "outcome"),
        "must have a single value per image")
    expect_error(stBoot:::runProportionModels(props, perImage(rep(c("a", "b"), c(4, 5))), "nope"),
        "`outcomeCol = \"nope\"` is not a column of the data")
    shared = perImage(rep(c("a", "b"), c(4, 5)), c("p1", "p1", paste0("p", 3:9)))
    expect_error(stBoot:::runProportionModels(props, shared, "outcome", subject = "subject"),
        "has more than one image for 1 subject: p1 \\(img1, img2\\)")
    expect_no_error(stBoot:::runProportionModels(props, perImage(rep(c("a", "b"), c(4, 5)),
        paste0("p", 1:9)), "outcome", subject = "subject"))
})

test_that("runProportionModels accepts column names that are not syntactic", {
    cells = perImage(rep(c("ctrl", "case"), c(4, 5)))
    names(cells)[2] = "my outcome"
    expect_equal(nrow(stBoot:::runProportionModels(props, cells, "my outcome")), 2)
})

## ---- the pipelines ------------------------------------------------------------
# confidences that are exactly the hard labels: the weighted analysis must then
# give exactly the unweighted results
hardConf = function() {
    cells = getCellData(d$spe)
    g = cells[cells$cellType %in% c("A", "B"), ]
    list(G = data.frame(cellID = g$cellID, A = as.numeric(g$cellType == "A"),
        B = as.numeric(g$cellType == "B")))
}

test_that("weighted results cover the same cell types as the unweighted ones", {
    res = quiet(runProportionsPipeline(d$spe, d$tree, d$conf, "outcome"))
    expect_named(res, c("proportions", "weighted", "unweighted",
        "modelWeighted", "modelUnweighted"))
    # unweighted: every cell type plus the bootstrapped group G
    expect_setequal(res$modelUnweighted$cellType, c("A", "B", "C", "G"))
    expect_setequal(res$modelWeighted$cellType, c("A", "B", "C", "G"))
    expect_setequal(res$weighted$cellType, c("A", "B", "C", "G"))
    expect_true(all(res$weighted$proportionWeighted >= 0))
    # cell types without confidences take the unweighted proportion
    byType = function(x, type, col) x[[col]][x$cellType == type]
    expect_equal(byType(res$weighted, "C", "proportionWeighted"),
        byType(res$unweighted, "C", "proportionUnweighted"))
    expect_equal(byType(res$weighted, "G", "proportionWeighted"),
        byType(res$unweighted, "G", "proportionUnweighted"))
    # and lineage is the group for subtypes, the type itself otherwise
    lin = unique(res$modelWeighted[, c("cellType", "lineage")])
    expect_equal(setNames(lin$lineage, lin$cellType)[c("A", "B", "C", "G")],
        c(A = "G", B = "G", C = "C", G = "G"))
})

test_that("with hard confidences the weighted and unweighted models are identical", {
    res = quiet(runProportionsPipeline(d$spe, d$tree, hardConf(), "outcome"))
    cols = c("cellType", "estimate", "statistic", "p.value", "p.adj")
    u = res$modelUnweighted[order(res$modelUnweighted$cellType), cols]
    w = res$modelWeighted[order(res$modelWeighted$cellType), cols]
    expect_equal(w, u, ignore_attr = TRUE)
})

test_that("pipelines do not report replaced group columns", {
    msgs = testthat::capture_messages(
        runProportionsPipeline(d$spe, d$tree, d$conf, "outcome"))
    expect_false(any(grepl("Replacing", msgs)))
    msgs = testthat::capture_messages(
        runProportionsRelativeToParentPipeline(d$spe, d$tree, d$conf, "outcome"))
    expect_false(any(grepl("Replacing", msgs)))
})

test_that("proportions relative to parent are within [0, 1]", {
    res = quiet(runProportionsRelativeToParentPipeline(d$spe, d$tree, d$conf,
        "outcome"))
    expect_true(all(res$proportions$proportionUnweighted <= 1, na.rm = TRUE))
    expect_true(all(res$proportions$proportionWeighted <= 1 + 1e-8, na.rm = TRUE))
    expect_named(res, c("proportions", "modelWeighted", "modelUnweighted"))
})

test_that("calcProportionsRelativeToParent works on annotated or plain cells", {
    cells = getCellData(d$spe)
    plain = stBoot:::calcProportionsRelativeToParent(cells, d$tree, d$conf)
    annotated = stBoot:::calcProportionsRelativeToParent(
        annotateCellTypeTree(cells, d$tree), d$tree, d$conf)
    expect_equal(plain, annotated)
})

test_that("the pipelines check their arguments before computing", {
    for (f in list(runProportionsPipeline, runProportionsRelativeToParentPipeline)) {
        spe = d$spe
        spe$three = rep(c("a", "a", "b", "b", "c", "c"), each = 70)
        msgs = testthat::capture_messages(
            expect_error(f(spe, d$tree, d$conf, "three"), "has 3 levels"))
        expect_false(any(grepl("Computing", msgs)))
        expect_error(f(spe, d$tree, d$conf, "nope"),
            "is not a column of colData\\(spe\\)")
        spe$patient = paste0("p", rep(c(1, 1, 2, 3, 4, 5), each = 70))
        expect_error(f(spe, d$tree, d$conf, "outcome", subject = "patient"),
            "has more than one image for 1 subject: p1 \\(img1, img2\\)")
        msgs = testthat::capture_messages(suppressWarnings(f(d$spe, d$tree,
            d$conf, "outcome")))
        expect_equal(sum(grepl("`subject` was not given", msgs)), 1)
    }
})

test_that("the pipeline arguments are validated", {
    pipelines = list(runProportionsPipeline, runProportionsRelativeToParentPipeline)
    for (f in pipelines) {
        # spe
        expect_error(f(getCellData(d$spe), d$tree, d$conf, "outcome"),
            "`spe` must be a SummarizedExperiment")
        expect_error(f("spe", d$tree, d$conf, "outcome"),
            "`spe` must be a SummarizedExperiment")

        # confList: structure
        msg = "`confList` must be a named list with one data.frame per group"
        for (bad in list(NULL, list(), d$conf[[1]], unname(d$conf), "G", 5)) {
            expect_error(f(d$spe, d$tree, bad, "outcome"), msg)
        }
        expect_error(f(d$spe, d$tree, list(G = 1), "outcome"),
            "`confList\\[\\[\"G\"\\]\\]` must be a data.frame with a `cellID` column")
        noId = list(G = d$conf$G[, c("A", "B")])
        expect_error(f(d$spe, d$tree, noId, "outcome"), "with a `cellID` column")
        onlyId = list(G = d$conf$G["cellID"])
        expect_error(f(d$spe, d$tree, onlyId, "outcome"), "with a `cellID` column")

        # confList: groups must be groups of the tree
        wrong = setNames(d$conf, "Typo")
        expect_error(f(d$spe, d$tree, wrong, "outcome"),
            "`confList` has 1 group that is not a group in cellTypeTree: 'Typo'. Groups in the tree: G")

        # skipGroups
        expect_error(f(d$spe, d$tree, d$conf, "outcome", skipGroups = 1),
            "`skipGroups` must be a character vector of group names")
        expect_error(f(d$spe, d$tree, d$conf, "outcome", skipGroups = NA_character_),
            "`skipGroups` must be a character vector")
        expect_error(f(d$spe, d$tree, d$conf, "outcome", skipGroups = "Typo"),
            "`skipGroups` has 1 name that is not a group in cellTypeTree: 'Typo'. Groups in the tree: G")
        expect_error(f(d$spe, d$tree, d$conf, "outcome", skipGroups = "G"),
            "Every group in `confList` is in `skipGroups`")

        # an invalid tree is reported before anything else
        expect_error(f(d$spe, list(), d$conf, "outcome"),
            "cellTypeTree must be a non-empty named list")
    }
})

test_that("skipGroups naming a group that confList lacks is allowed", {
    # the tree has a second group, as when it was skipped in bootstrapping
    tree2 = list(G = list(children = c("A", "B")), H = list(children = c("C", "D")))
    for (f in list(runProportionsPipeline, runProportionsRelativeToParentPipeline)) {
        expect_no_error(quiet(f(d$spe, tree2, d$conf, "outcome", skipGroups = "H")))
    }
})

test_that("the arguments are checked before any computation", {
    msgs = testthat::capture_messages(
        expect_error(runProportionsPipeline(d$spe, d$tree, d$conf, "outcome",
            skipGroups = "Typo")))
    expect_false(any(grepl("Computing", msgs)))
})

test_that("the pipelines take custom column names", {
    for (f in list(runProportionsPipeline, runProportionsRelativeToParentPipeline)) {
        default = quiet(f(d$spe, d$tree, d$conf, "outcome"))
        custom = quiet(f(renamedSpe(d$spe), d$tree, d$conf, "outcome",
            cellType = "ct", cellID = "barcode", imageID = "sample"))
        expect_equal(custom, default)
        expect_error(quiet(f(renamedSpe(d$spe), d$tree, d$conf, "outcome")),
            "`cellID = \"cellID\"` is not a column of the data")
    }
})

test_that("the dead weighted-proportion helper is gone", {
    expect_false("calcWeightedProportions" %in% getNamespaceExports("stBoot"))
})

test_that("proportion pipelines return models (full data)", {
    skip_unless_full()
    f = fullTestData()
    res = suppressMessages(runProportionsPipeline(f$spe, f$heartCellTypeTree,
        f$confList, outcomeCol = "resolution_broad", subject = "subject"))
    expect_named(res, c("proportions", "weighted", "unweighted",
        "modelWeighted", "modelUnweighted"))
    expect_setequal(res$modelWeighted$cellType, res$modelUnweighted$cellType)
    expect_true(all(res$weighted$proportionWeighted >= 0))

    rel = suppressMessages(runProportionsRelativeToParentPipeline(f$spe,
        f$heartCellTypeTree, f$confList, outcomeCol = "resolution_broad",
        subject = "subject"))
    expect_true(all(rel$proportions$proportionUnweighted <= 1, na.rm = TRUE))
})
