## End-to-end run of the spatial functions on a small hand-made dataset, so
## they are exercised even when the full-data tests are skipped.
d = miniPipelineData()
serial = BiocParallel::SerialParam()
quiet = function(x) suppressWarnings(suppressMessages(x))

test_that("runColocalisationModels runs end to end", {
    res = quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome"))
    expect_named(res, c("weighted", "unweighted", "modelWeighted",
        "modelUnweighted", "relativeWeighted", "relativeUnweighted",
        "modelRelativeWeighted", "modelRelativeUnweighted"))
    expect_setequal(unique(res$weighted$imageID), paste0("img", 1:6))
    expect_true(all(c("from", "to", "estimate", "p.adj") %in%
        names(res$modelWeighted)))
})

test_that("niche identification runs end to end", {
    set.seed(1)
    res = quiet(runNicheIdentification(d$spe, d$tree, d$conf, r = c(20, 50),
        k = 3))
    expect_length(res$nicheWeighted, ncol(d$spe))
    expect_lte(length(unique(res$nicheWeighted)), 3)

    tst = quiet(testNicheOutcome(d$spe, res$nicheUnweighted,
        res$nicheWeighted, "outcome"))
    expect_named(tst, c("testUnweighted", "testWeighted", "propUnweighted",
        "propWeighted"))
})

test_that("proportion pipelines run end to end", {
    res = quiet(runProportionsPipeline(d$spe, d$tree, d$conf, "outcome"))
    expect_true(all(res$weighted$proportionWeighted >= 0))
    rel = quiet(runProportionsRelativeToParentPipeline(d$spe, d$tree,
        d$conf, "outcome"))
    expect_true(all(rel$proportions$proportionUnweighted <= 1, na.rm = TRUE))
})

test_that("results do not depend on the BiocParallel backend", {
    skip_on_os("windows")
    multi = BiocParallel::MulticoreParam(2)
    expect_equal(
        quiet(bootstrapLineageConfidence(d$spe, d$tree, nBoot = 3, seed = 1,
            BPPARAM = serial)),
        quiet(bootstrapLineageConfidence(d$spe, d$tree, nBoot = 3, seed = 1,
            BPPARAM = multi)))
    expect_equal(
        quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome",
            BPPARAM = serial)),
        quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome",
            BPPARAM = multi)))
})
