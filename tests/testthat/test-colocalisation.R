d = miniPipelineData()
quiet = function(x) suppressWarnings(suppressMessages(x))

## ---- outcome checks ---------------------------------------------------
Lf = data.frame(imageID = paste0("img", 1:9), from = "A", to = "B",
    colocalisation = c(1, 2, 3, 4, 5, 6, 7, 8, 9))
perImage = function(outcome) {
    data.frame(imageID = rep(paste0("img", 1:9), each = 2),
        outcome = rep(outcome, each = 2))
}

test_that("runLModels fits a model for two conditions", {
    cells = perImage(rep(c("ctrl", "case"), c(4, 5)))
    res = fuzzyCell:::runLModels(Lf, cells, "outcome")
    expect_equal(nrow(res), 1)
    # estimate is the second level (alphabetically: 'ctrl') minus the first
    expect_equal(res$estimate, mean(Lf$colocalisation[1:4]) - mean(Lf$colocalisation[5:9]))
})

test_that("runLModels accepts column names that are not syntactic", {
    cells = perImage(rep(c("ctrl", "case"), c(4, 5)))
    names(cells)[2] = "my outcome"
    expect_equal(nrow(fuzzyCell:::runLModels(Lf, cells, "my outcome")), 1)
})

test_that("outcomes with more than two levels are rejected", {
    cells = perImage(rep(c("low", "mid", "high"), each = 3))
    expect_error(fuzzyCell:::runLModels(Lf, cells, "outcome"),
        "has 3 levels \\(high, low, mid\\)\\. The models currently compare two conditions only\\. Please subset the data to two conditions")
    # a continuous outcome has many levels
    expect_error(fuzzyCell:::runLModels(Lf, perImage(c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5)),
        "outcome"), "has 9 levels")
})

test_that("outcomes with fewer than two levels are rejected", {
    expect_error(fuzzyCell:::runLModels(Lf, perImage(rep("only", 9)), "outcome"),
        "fewer than two levels \\(only only\\)")
    expect_error(fuzzyCell:::runLModels(Lf, perImage(rep(NA_character_, 9)), "outcome"),
        "fewer than two levels")
})

test_that("each image must have a single outcome value", {
    cells = data.frame(imageID = rep(paste0("img", 1:9), each = 2),
        outcome = rep(c("a", "b"), 9))
    expect_error(fuzzyCell:::runLModels(Lf, cells, "outcome"),
        "must have a single value per image, but 9 images have more than one: img1 \\(a, b\\); img2 \\(a, b\\); img3 \\(a, b\\); img4 \\(a, b\\); img5 \\(a, b\\); \\.\\.\\.")
    expect_error(fuzzyCell:::runLModels(Lf, cells, "outcome"),
        "repeated measures and mixed models are not supported. Please make sure there is exactly one outcome value per image, and one image per sample \\(for example one image per patient\\)\\.")

    # a missing value alongside one real value is fine
    ok = perImage(rep(c("a", "b"), c(4, 5)))
    ok$outcome[1] = NA
    expect_no_error(fuzzyCell:::runLModels(Lf, ok, "outcome"))
})

test_that("outcomeCol must be a column", {
    cells = perImage(rep(c("a", "b"), c(4, 5)))
    expect_error(fuzzyCell:::runLModels(Lf, cells, "nope"),
        "`outcomeCol = \"nope\"` is not a column of the data. Available columns: imageID, outcome")
    expect_error(fuzzyCell:::runLModels(Lf, cells, c("outcome", "imageID")),
        "`outcomeCol` must be a single string")
})

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

## ---- weightBy ---------------------------------------------------------
cells = annotateCellTypeTree(getCellData(d$spe), d$tree)
cells$group[is.na(cells$group)] = cells$cellType[is.na(cells$group)]
pairs = data.frame(from = "A", to = "B")

test_that("weightBy must be one of both, i or j", {
    expect_error(fuzzyCell:::calcLfuncPerImage(cells, d$conf, pairs, weightBy = "bth"),
        "'arg' should be one of")
    expect_error(fuzzyCell:::calcKontextualPerImage(cells, d$conf, pairs, d$tree,
        weightBy = "bth"), "'arg' should be one of")
    expect_no_error(quiet(fuzzyCell:::calcLfuncPerImage(cells, d$conf, pairs,
        weightBy = "j")))
})

## ---- images that are too small ---------------------------------------
corners = function(id, w, h) {
    data.frame(imageID = id, x = c(0, w), y = c(0, h))
}

test_that("smallImages drops images with a side below 2 * r", {
    imgs = rbind(corners("big", 300, 300), corners("narrow", 150, 400),
        corners("short", 400, 150), corners("exact", 200, 200))
    # limit 100: nothing is too small
    expect_equal(nrow(fuzzyCell:::smallImages(imgs, r = 50)), 0)
    # limit 200: an image exactly at the limit stays
    expect_setequal(fuzzyCell:::smallImages(imgs, r = 100)$imageID,
        c("narrow", "short"))
    # limit 302: all of them go
    expect_setequal(fuzzyCell:::smallImages(imgs, r = 151)$imageID,
        c("big", "exact", "narrow", "short"))
})

test_that("dropSmallImages reports what it drops", {
    imgs = rbind(corners("big", 300, 300), corners("narrow", 150, 400),
        corners("fine", 500, 500))
    expect_message(out <- fuzzyCell:::dropSmallImages(imgs, r = 100),
        "Dropping 1 of 3 images smaller than 2 x r = 200 units in at least one direction: narrow \\(150 x 400\\)")
    expect_setequal(unique(out$imageID), c("big", "fine"))
    expect_no_message(fuzzyCell:::dropSmallImages(out, r = 100))
    expect_error(fuzzyCell:::dropSmallImages(corners("tiny", 50, 50), r = 50),
        "All 1 images are smaller than 2 x r = 100 units")
})

test_that("calcLfuncPerImage drops small images with a message", {
    tiny = cells[cells$imageID == "img1", ]
    tiny$imageID = "tiny"
    tiny$cellID = paste0("t_", tiny$cellID)
    tiny$x = tiny$x / 3
    tiny$y = tiny$y / 3
    withTiny = rbind(cells, tiny)
    # the tiny image is about 100 units wide, so 2 * r = 120 drops it
    expect_message(res <- suppressWarnings(fuzzyCell:::calcLfuncPerImage(withTiny, d$conf,
        pairs, r = 60)), "Dropping 1 of 7 images smaller than 2 x r = 120 units")
    expect_false("tiny" %in% res$weighted$imageID)
    expect_setequal(unique(res$weighted$imageID), paste0("img", 1:6))
})

test_that("runColocalisationModels reports dropped images once", {
    tiny = d$spe[, d$spe$imageID == "img1"]
    SpatialExperiment::spatialCoords(tiny) =
        SpatialExperiment::spatialCoords(tiny) / 3
    tiny$imageID = "tiny"
    colnames(tiny) = paste0("t_", colnames(tiny))
    tiny$cellID = colnames(tiny)
    spe = BiocGenerics::cbind(d$spe, tiny)
    msgs = testthat::capture_messages(
        res <- suppressWarnings(runColocalisationModels(spe, d$tree, d$conf, "outcome",
            r = 60)))
    expect_equal(sum(grepl("Dropping", msgs)), 1)
    expect_false("tiny" %in% res$weighted$imageID)
})

## ---- subject: one image per subject ---------------------------------------
subjectCells = function(subjects) {
    cells = perImage(rep(c("ctrl", "case"), c(4, 5)))
    cells$subject = rep(subjects, each = 2)
    cells
}

test_that("one image per subject passes", {
    cells = subjectCells(paste0("p", 1:9))
    expect_equal(nrow(fuzzyCell:::runLModels(Lf, cells, "outcome", subject = "subject")), 1)
    # subject = NULL skips the check
    shared = subjectCells(c("p1", "p1", paste0("p", 3:9)))
    expect_no_error(fuzzyCell:::runLModels(Lf, shared, "outcome"))
})

test_that("several images from one subject are rejected", {
    cells = subjectCells(c("p1", "p1", "p2", "p3", "p3", "p3", paste0("p", 4:6)))
    expect_error(fuzzyCell:::runLModels(Lf, cells, "outcome", subject = "subject"),
        "`subject = \"subject\"` has more than one image for 2 subjects: p1 \\(img1, img2\\); p3 \\(img4, img5, img6\\)")
    expect_error(fuzzyCell:::runLModels(Lf, cells, "outcome", subject = "subject"),
        "mixed models are not supported\\)\\. Please keep one image per subject, for example by subsetting to one image per patient")
})

test_that("an image with more than one subject is rejected", {
    cells = subjectCells(paste0("p", 1:9))
    cells$subject[2] = "other"   # img1 now has two subjects
    expect_error(fuzzyCell:::runLModels(Lf, cells, "outcome", subject = "subject"),
        "Every image must belong to a single subject, but 1 image has more than one value in the `subject` column 'subject': img1 \\(p1, other\\)")
})

test_that("missing or unknown subject columns are rejected", {
    cells = subjectCells(paste0("p", 1:9))
    cells$subject[3] = NA
    expect_error(fuzzyCell:::runLModels(Lf, cells, "outcome", subject = "subject"),
        "has missing values \\(in 1 image: img2\\)")
    expect_error(fuzzyCell:::runLModels(Lf, cells, "outcome", subject = "patient"),
        "`subject = \"patient\"` is not a column of the data")
    expect_error(fuzzyCell:::runLModels(Lf, cells, "outcome", subject = c("a", "b")),
        "`subject` must be a single string")
})

test_that("runColocalisationModels checks subjects before computing anything", {
    spe = d$spe
    # images 1 and 2 are from the same patient
    spe$patient = paste0("p", rep(c(1, 1, 2, 3, 4, 5), each = 70))
    msgs = testthat::capture_messages(
        expect_error(runColocalisationModels(spe, d$tree, d$conf, "outcome",
            subject = "patient"), "has more than one image for 1 subject: p1 \\(img1, img2\\)"))
    expect_false(any(grepl("Generating comparisons|Computing", msgs)))

    expect_error(runColocalisationModels(spe, d$tree, d$conf, "outcome",
        subject = "nope"), "is not a column of colData\\(spe\\)")

    # one image per subject: runs, and without the 'not checked' message
    spe$patient = paste0("p", rep(1:6, each = 70))
    msgs = testthat::capture_messages(
        res <- suppressWarnings(runColocalisationModels(spe, d$tree, d$conf, "outcome",
            subject = "patient")))
    expect_false(any(grepl("`subject` was not given", msgs)))
    expect_named(res, c("weighted", "unweighted", "modelWeighted",
        "modelUnweighted", "relativeWeighted", "relativeUnweighted",
        "modelRelativeWeighted", "modelRelativeUnweighted"))
})

test_that("runColocalisationModels says when subject was not given", {
    msgs = testthat::capture_messages(
        suppressWarnings(runColocalisationModels(d$spe, d$tree, d$conf, "outcome")))
    expect_equal(sum(grepl("`subject` was not given", msgs)), 1)
})

test_that("images dropped as too small do not count against a subject", {
    # img7 is a tiny second image from patient p1
    tiny = d$spe[, d$spe$imageID == "img1"]
    SpatialExperiment::spatialCoords(tiny) =
        SpatialExperiment::spatialCoords(tiny) / 3
    tiny$imageID = "tiny"
    colnames(tiny) = paste0("t_", colnames(tiny))
    tiny$cellID = colnames(tiny)
    spe = BiocGenerics::cbind(d$spe, tiny)
    spe$patient = paste0("p", c(rep(1:6, each = 70), rep(1, 70)))
    expect_no_error(suppressWarnings(suppressMessages(
        runColocalisationModels(spe, d$tree, d$conf, "outcome", subject = "patient",
            r = 60))))
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
