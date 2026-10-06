d = miniPipelineData()
quiet = function(x) suppressWarnings(suppressMessages(x))

# a tree with two groups, so that two groups can be plotted together
tree2 = list(G = list(children = c("A", "B")), H = list(children = c("C", "D")))
cells = getCellData(d$spe)
cc = cells[cells$cellType == "C", ]
conf = list(G = d$conf$G, H = data.frame(cellID = cc$cellID, C = 1, D = 0))
niches = quiet(runNicheIdentification(d$spe, d$tree, d$conf, r = c(20, 50), k = 3))
spe = d$spe
spe$nicheWeighted = niches$nicheWeighted[spe$cellID]

# the plots that draw individual cells, each as a function of cellSize
cellPlots = list(
    plotOriginalLabelSpatial = function(...) plotOriginalLabelSpatial(spe, tree2, "img1", "G", ...),
    plotHardSpatialPair = function(...) plotHardSpatialPair(spe, "img1", "A", "B", tree2, ...),
    plotWeightedSpatial = function(...) plotWeightedSpatial(spe, tree2, conf, "img1", "G", "A", ...),
    plotWeightedSpatialPair_alpha = function(...) plotWeightedSpatialPair(spe, tree2, conf, "img1", "G", "A", "H", "C", ...),
    plotWeightedSpatialPair_colour = function(...) plotWeightedSpatialPair(spe, tree2, conf, "img1", "G", "A", "H", "C", encoding = "colour", ...),
    plotOneMarkerSpatial = function(...) plotOneMarkerSpatial(spe, "img1", "G", "g1", cellTypeTree = tree2, ...),
    plotTwoMarkerSpatial = function(...) plotTwoMarkerSpatial(spe, "img1", "G", "g1", "purple", "H", "g4", "green", cellTypeTree = tree2, ...),
    plotNicheSpatial = function(...) plotNicheSpatial(spe, "img1", "niche1", ...)
)

# point size of every layer of a plot, background layer first
layerSizes = function(p) {
    ggplot2::ggplot_build(p)
    unname(vapply(p$layers, function(l) l$aes_params$size, numeric(1)))
}

test_that("every plot takes a SummarizedExperiment", {
    plots = list(
        function(x) plotOriginalLabelSpatial(x, tree2, "img1", "G"),
        function(x) plotHardSpatialPair(x, "img1", "A", "B", tree2),
        function(x) plotWeightedSpatial(x, tree2, conf, "img1", "G", "A"),
        function(x) plotWeightedSpatialPair(x, tree2, conf, "img1", "G", "A", "H", "C"),
        function(x) plotOneMarkerSpatial(x, "img1", "G", "g1", cellTypeTree = tree2),
        function(x) plotTwoMarkerSpatial(x, "img1", "G", "g1", "purple", "H", "g4", "green", cellTypeTree = tree2),
        function(x) plotNicheSpatial(x, "img1", "niche1"),
        function(x) nicheCompositionMap(x, "nicheWeighted"),
        function(x) nicheBarplot(x, "nicheWeighted", "niche1")
    )
    for (f in plots) {
        # a plain data frame is rejected, with the same message everywhere
        expect_error(f(cells),
            "`spe` must be a SummarizedExperiment, for example a SpatialExperiment")
        expect_error(f("spe"), "`spe` must be a SummarizedExperiment")
        expect_no_error(quiet(f(spe)))
    }
})

test_that("highlighted cells have size 1 and the grey background 0.7 by default", {
    for (nm in names(cellPlots)) {
        sizes = layerSizes(cellPlots[[nm]]())
        expect_equal(sizes[1], 0.7, info = nm)
        expect_true(all(sizes[-1] == 1), info = nm)
        expect_gt(length(sizes), 1)
    }
})

test_that("cellSize sets the size of the highlighted cells only", {
    for (nm in names(cellPlots)) {
        sizes = layerSizes(cellPlots[[nm]](cellSize = 2.5))
        expect_equal(sizes[1], 0.7, info = nm)
        expect_true(all(sizes[-1] == 2.5), info = nm)
    }
})

test_that("cellSize must be a single positive number", {
    for (nm in names(cellPlots)) {
        for (bad in list(0, -1, NA, "1", c(1, 2), Inf, NULL)) {
            expect_error(cellPlots[[nm]](cellSize = bad),
                "`cellSize` must be a single positive number", info = nm)
        }
    }
})

test_that("plots use the data in colData even when it has its own group column", {
    # a dataset that already has `group` and `lineage` columns of its own
    own = spe
    own$group = "mine"
    own$lineage = "mine"
    expect_no_message(plotOriginalLabelSpatial(own, tree2, "img1", "G"))
    expect_no_message(plotWeightedSpatial(own, tree2, conf, "img1", "G", "A"))
})

test_that("plotWeightedSpatialPair has a separate legend for each cell type", {
    for (enc in c("alpha", "colour")) {
        p = plotWeightedSpatialPair(spe, tree2, conf, "img1", "G", "A", "H", "C",
            encoding = enc)
        scaleNames = vapply(p$scales$scales, function(s) {
            if (is.character(s$name)) s$name[1] else NA_character_
        }, character(1))
        expect_true("A\nconfidence" %in% scaleNames, info = enc)
        expect_true("C\nconfidence" %in% scaleNames, info = enc)
    }
    # alpha is the default encoding
    p = plotWeightedSpatialPair(spe, tree2, conf, "img1", "G", "A", "H", "C")
    expect_true(any(grepl("^alpha", vapply(p$scales$scales,
        function(s) s$aesthetics[1], ""))))
})

## ---- marker plots read the counts from an assay of spe ------------------------
test_that("the marker plots take their counts from the assay of spe", {
    # a second assay with different counts is used when asked for
    other = spe
    SummarizedExperiment::assay(other, "doubled") = 2 * SummarizedExperiment::assay(spe, "counts")
    expect_no_error(plotOneMarkerSpatial(other, "img1", "G", "g1",
        cellTypeTree = tree2, assayName = "doubled"))
    expect_error(plotOneMarkerSpatial(spe, "img1", "G", "g1", cellTypeTree = tree2,
        assayName = "logcounts"),
        "`assayName = \"logcounts\"` is not an assay of spe. Available assays: counts")
    expect_error(plotTwoMarkerSpatial(spe, "img1", "G", "g1", "purple", "H", "g4",
        "green", cellTypeTree = tree2, assayName = 1),
        "`assayName` must be a single string")
})

test_that("the marker plots check the genes", {
    expect_error(plotOneMarkerSpatial(spe, "img1", "G", "nogene", cellTypeTree = tree2),
        "`marker = \"nogene\"` is not a gene \\(a row\\) of the 'counts' assay")
    expect_error(plotOneMarkerSpatial(spe, "img1", "G", c("g1", "g2"), cellTypeTree = tree2),
        "`marker` must be a single string naming a gene")
    expect_error(plotTwoMarkerSpatial(spe, "img1", "G", "g1", "purple", "H", "nogene",
        "green", cellTypeTree = tree2), "`markerB = \"nogene\"` is not a gene")
})

test_that("the marker plots need a tree to find the groups", {
    expect_error(plotOneMarkerSpatial(spe, "img1", "G", "g1"),
        "`cellTypeTree` is needed to find the cells of the group")
    expect_error(plotTwoMarkerSpatial(spe, "img1", "G", "g1", "purple", "H", "g4", "green"),
        "`cellTypeTree` is needed to find the cells of the groups")
})

## ---- niche plots ----------------------------------------------------------------
test_that("the niche plots check their columns and niche", {
    expect_error(plotNicheSpatial(spe, "img1", "niche1", nicheCol = "nope"),
        "`nicheCol = \"nope\"` is not a column of the data")
    expect_error(nicheCompositionMap(spe, "nope"),
        "`niche = \"nope\"` is not a column of the data")
    expect_error(nicheCompositionMap(spe, "nicheWeighted", cellType = "nope"),
        "`cellType = \"nope\"` is not a column of the data")
    expect_error(nicheBarplot(spe, "nicheWeighted", "niche99"),
        "`target` must be one of the niches in `nicheWeighted`: niche1, niche2, niche3")
})

## ---- plotMetricBoxplot ----------------------------------------------------------
test_that("plotMetricBoxplot takes lfunction or kontextual", {
    skip_if_not_installed("ggiraph")
    res = quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome"))
    for (metric in c("lfunction", "kontextual")) {
        p = plotMetricBoxplot(res, d$spe, "A", "B", "outcome", metric = metric)
        expect_s3_class(p, "girafe")

        # boxplots and points, with no lines joining the points
        gg = scBoot:::metricBoxplotGG(res, cells, "A", "B", "outcome", metric)
        geoms = vapply(gg$layers, function(l) class(l$geom)[1], "")
        expect_gt(length(geoms), 0)
        expect_true("GeomBoxplot" %in% geoms)
        expect_false("GeomLine" %in% geoms)
    }
    # the default is lfunction
    expect_equal(eval(formals(plotMetricBoxplot)$metric)[1], "lfunction")
    # other names are rejected
    expect_error(plotMetricBoxplot(res, d$spe, "A", "B", "outcome",
        metric = "colocalisation"), "'arg' should be one of")
    # and so is anything that is not an spe
    expect_error(plotMetricBoxplot(res, cells, "A", "B", "outcome"),
        "`spe` must be a SummarizedExperiment")
})

test_that("the removed plot functions are no longer exported", {
    exports = getNamespaceExports("scBoot")
    expect_false("plotConfidenceMarkers" %in% exports)
    expect_false("plotPValComparison" %in% exports)
    expect_true("testNicheOutcome" %in% exports)
})

## ---- plotProportionBoxplot ------------------------------------------------------
test_that("plotProportionBoxplot handles both proportion pipelines", {
    skip_if_not_installed("ggiraph")
    pipelines = list(
        all = quiet(runProportionsPipeline(d$spe, d$tree, d$conf, "outcome")),
        parent = quiet(runProportionsRelativeToParentPipeline(d$spe, d$tree,
            d$conf, "outcome"))
    )
    for (nm in names(pipelines)) {
        res = pipelines[[nm]]
        p = plotProportionBoxplot(res, "A", d$spe, "outcome")
        expect_s3_class(p, "girafe")

        gg = scBoot:::proportionBoxplotGG(res, "A", cells, "outcome")
        geoms = vapply(gg$layers, function(l) class(l$geom)[1], "")
        expect_true("GeomBoxplot" %in% geoms, info = nm)
        expect_false("GeomLine" %in% geoms, info = nm)
        # one point per image and weighting
        expect_equal(nrow(gg$data), 2 * 6, info = nm)
        expect_equal(levels(gg$data$type), c("Unweighted", "Weighted"), info = nm)
    }
    # the axis says what the proportion is of
    expect_equal(scBoot:::proportionBoxplotGG(pipelines$all, "A", cells,
        "outcome")$labels$y, "Proportion of all cells")
    expect_equal(scBoot:::proportionBoxplotGG(pipelines$parent, "A", cells,
        "outcome")$labels$y, "Proportion of G cells")
    # whole groups and ungrouped cell types are available in the main pipeline
    expect_no_error(plotProportionBoxplot(pipelines$all, "G", d$spe, "outcome"))
    expect_no_error(plotProportionBoxplot(pipelines$all, "C", d$spe, "outcome"))
})

test_that("plotProportionBoxplot checks its arguments", {
    skip_if_not_installed("ggiraph")
    res = quiet(runProportionsPipeline(d$spe, d$tree, d$conf, "outcome"))
    expect_error(plotProportionBoxplot(list(), "A", d$spe, "outcome"),
        "`results` must be the output of runProportionsPipeline")
    expect_error(plotProportionBoxplot(res$proportions, "A", d$spe, "outcome"),
        "`results` must be the output of runProportionsPipeline")
    expect_error(plotProportionBoxplot(res, "Typo", d$spe, "outcome"),
        "`cellType = \"Typo\"` is not in the results. Available: A, B, C, G")
    expect_error(plotProportionBoxplot(res, c("A", "B"), d$spe, "outcome"),
        "`cellType` must be a single string")
    expect_error(plotProportionBoxplot(res, "A", d$spe, "nope"),
        "`outcomeCol = \"nope\"` is not a column of the data")
    expect_error(plotProportionBoxplot(res, "A", cells, "outcome"),
        "`spe` must be a SummarizedExperiment")
})

## ---- nicheProportionBoxplot -----------------------------------------------------
test_that("nicheProportionBoxplot plots the table it is given", {
    skip_if_not_installed("ggiraph")
    tests = quiet(testNicheOutcome(spe, niches$nicheUnweighted, niches$nicheWeighted,
        "outcome"))
    for (type in c("weighted", "unweighted")) {
        props = tests[[paste0("prop", tools::toTitleCase(type))]]
        expect_equal(attr(props, "nicheType"), type)
        expect_s3_class(nicheProportionBoxplot(props, "niche1", spe, "outcome"), "girafe")

        # the table decides which niches are plotted, and the title says which
        gg = scBoot:::nicheProportionBoxplotGG(props, "niche1", cells, "outcome")
        expect_equal(gg$labels$title, paste0("niche1 (", type, ")"))
        expect_equal(nrow(gg$data), sum(props$niche == "niche1"))
        geoms = vapply(gg$layers, function(l) class(l$geom)[1], "")
        expect_true("GeomBoxplot" %in% geoms)
    }
    # there is no weighted argument any more
    expect_false("weighted" %in% names(formals(nicheProportionBoxplot)))

    # the two tables really are different plots
    w = scBoot:::nicheProportionBoxplotGG(tests$propWeighted, "niche1", cells, "outcome")
    u = scBoot:::nicheProportionBoxplotGG(tests$propUnweighted, "niche1", cells, "outcome")
    expect_false(isTRUE(all.equal(w$data$proportion, u$data$proportion)) &&
        w$labels$title == u$labels$title)
})

test_that("nicheProportionBoxplot checks its arguments", {
    skip_if_not_installed("ggiraph")
    tests = quiet(testNicheOutcome(spe, niches$nicheUnweighted, niches$nicheWeighted,
        "outcome"))
    msg = "`proportions` must be the `propWeighted` or `propUnweighted` element of the results of testNicheOutcome\\(\\)"
    # the whole result list is no longer accepted
    expect_error(nicheProportionBoxplot(tests, "niche1", spe, "outcome"), msg)
    expect_error(nicheProportionBoxplot(NULL, "niche1", spe, "outcome"), msg)
    expect_error(nicheProportionBoxplot(data.frame(a = 1), "niche1", spe, "outcome"), msg)
    expect_error(nicheProportionBoxplot(tests$propWeighted, "niche99", spe, "outcome"),
        "`niche` must be one of the niches in `proportions`: niche1, niche2, niche3")
    expect_error(nicheProportionBoxplot(tests$propWeighted, "niche1", cells, "outcome"),
        "`spe` must be a SummarizedExperiment")
    expect_error(nicheProportionBoxplot(tests$propWeighted, "niche1", spe, "nope"),
        "`outcomeCol = \"nope\"` is not a column of the data")
})

## ---- the image, groups, subtypes, niche and pair are checked ---------------------
test_that("every single-image plot checks that the image exists", {
    msg = "`imageID = \"img99\"` is not an image in the data. Images: img1, img2, img3, img4, img5, img6"
    expect_error(plotOriginalLabelSpatial(spe, tree2, "img99", "G"), msg)
    expect_error(plotHardSpatialPair(spe, "img99", "A", "B", tree2), msg)
    expect_error(plotWeightedSpatial(spe, tree2, conf, "img99", "G", "A"), msg)
    expect_error(plotWeightedSpatialPair(spe, tree2, conf, "img99", "G", "A", "H", "C"), msg)
    expect_error(plotOneMarkerSpatial(spe, "img99", "G", "g1", cellTypeTree = tree2), msg)
    expect_error(plotTwoMarkerSpatial(spe, "img99", "G", "g1", "purple", "H", "g4",
        "green", cellTypeTree = tree2), msg)
    expect_error(plotNicheSpatial(spe, "img99", "niche1"), msg)
    expect_error(plotOriginalLabelSpatial(spe, tree2, c("img1", "img2"), "G"),
        "`imageID` must be a single string")
    expect_error(plotNicheSpatial(spe, NA_character_, "niche1"),
        "`imageID` must be a single string")
})

test_that("groups are checked against the data", {
    msg = function(arg, value) paste0("`", arg, " = \"", value,
        "\"` is not a group in the data. Groups: G, H")
    expect_error(plotOriginalLabelSpatial(spe, tree2, "img1", "GG"), msg("targetGroup", "GG"))
    expect_error(plotWeightedSpatial(spe, tree2, conf, "img1", "GG", "A"), msg("targetGroup", "GG"))
    expect_error(plotWeightedSpatialPair(spe, tree2, conf, "img1", "G", "A", "HH", "C"),
        msg("groupB", "HH"))
    expect_error(plotWeightedSpatialPair(spe, tree2, conf, "img1", "Gx", "A", "HH", "C"),
        msg("groupA", "Gx"))
    expect_error(plotOneMarkerSpatial(spe, "img1", "GG", "g1", cellTypeTree = tree2),
        msg("group", "GG"))
    expect_error(plotTwoMarkerSpatial(spe, "img1", "G", "g1", "purple", "HH", "g4",
        "green", cellTypeTree = tree2), msg("groupB", "HH"))
    # a cell type is not a group
    expect_error(plotOriginalLabelSpatial(spe, tree2, "img1", "A"), msg("targetGroup", "A"))
    expect_error(plotOriginalLabelSpatial(spe, tree2, "img1", c("G", "H")),
        "`targetGroup` must be a single group name")
    # a group of the tree without any cells in the data is not a group of the data
    emptyTree = list(G = list(children = c("A", "B")), H = list(children = c("X", "Y")))
    expect_error(plotOriginalLabelSpatial(spe, emptyTree, "img1", "H"),
        "`targetGroup = \"H\"` is not a group in the data. Groups: G")
})

test_that("labels of plotHardSpatialPair can be cell types or groups, and must exist", {
    expect_no_error(plotHardSpatialPair(spe, "img1", "A", "H", tree2))
    expect_no_error(plotHardSpatialPair(spe, "img1", "G", "C", tree2))
    expect_error(plotHardSpatialPair(spe, "img1", "A", "Typo", tree2),
        "`subtypeB = \"Typo\"` is neither a cell type nor a group in the data")
    expect_error(plotHardSpatialPair(spe, "img1", "Typo", "B", tree2),
        "`subtypeA = \"Typo\"` is neither a cell type nor a group in the data")
    expect_error(plotHardSpatialPair(spe, "img1", c("A", "B"), "B", tree2),
        "`subtypeA` must be a single cell type or group name")
})

test_that("subtypes need confidence scores in confList", {
    expect_error(plotWeightedSpatial(spe, tree2, conf, "img1", "G", "Typo"),
        "`targetSubtype = \"Typo\"` is not a subtype of 'G' in `confList`. Subtypes: A, B")
    expect_error(plotWeightedSpatialPair(spe, tree2, conf, "img1", "G", "A", "H", "Typo"),
        "`subtypeB = \"Typo\"` is not a subtype of 'H' in `confList`. Subtypes: C, D")
    expect_error(plotWeightedSpatialPair(spe, tree2, conf, "img1", "G", "Typo", "H", "C"),
        "`subtypeA = \"Typo\"` is not a subtype of 'G'")
    expect_error(plotWeightedSpatial(spe, tree2, conf["G"], "img1", "H", "C"),
        "`confList` has no confidence scores for the group 'H'")
    # the group itself stands for the whole group
    expect_no_error(plotWeightedSpatial(spe, tree2, conf, "img1", "G", "G"))
})

test_that("plotNicheSpatial checks that the niche exists", {
    expect_error(plotNicheSpatial(spe, "img1", "niche9"),
        "`niche = \"niche9\"` is not a niche in `nicheWeighted`. Niches: niche1, niche2, niche3")
    expect_error(plotNicheSpatial(spe, "img1", c("niche1", "niche2")),
        "`niche` must be a single niche name")
    expect_no_error(plotNicheSpatial(spe, "img1", "niche2"))
})

test_that("plotMetricBoxplot checks the pair, the results and the outcome", {
    skip_if_not_installed("ggiraph")
    res = quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome"))
    expect_error(plotMetricBoxplot(res, d$spe, "A", "Zzz", "outcome"),
        "The pair 'A' -> 'Zzz' is not in the results$")
    expect_error(plotMetricBoxplot(res, d$spe, "Zzz", "B", "outcome"),
        "The pair 'Zzz' -> 'B' is not in the results")
    # relative co-localisation only has pairs whose `to` has a parent
    expect_error(plotMetricBoxplot(res, d$spe, "A", "G", "outcome", metric = "kontextual"),
        "The pair 'A' -> 'G' is not in the results for metric = \"kontextual\". Relative co-localisation only has pairs whose `toType` has a parent other than itself")
    expect_error(plotMetricBoxplot(list(), d$spe, "A", "B", "outcome"),
        "`results` must be the output of runColocalisationModels\\(\\)")
    expect_error(plotMetricBoxplot(res$modelWeighted, d$spe, "A", "B", "outcome"),
        "`results` must be the output of runColocalisationModels\\(\\)")
    expect_error(plotMetricBoxplot(res, d$spe, c("A", "B"), "B", "outcome"),
        "`fromType` must be a single cell type or group name")
    expect_error(plotMetricBoxplot(res, d$spe, "A", "B", "nope"),
        "`outcomeCol = \"nope\"` is not a column of the data")

    # outcomeCol has no default any more
    expect_true("outcomeCol" %in% names(formals(plotMetricBoxplot)))
    expect_true(identical(formals(plotMetricBoxplot)$outcomeCol, quote(expr = )))
    expect_error(plotMetricBoxplot(res, d$spe, "A", "B"), "outcomeCol")
})
