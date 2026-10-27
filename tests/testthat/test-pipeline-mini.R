## End-to-end run of the spatial functions on a small hand-made dataset, so
## they are exercised even when the full-data tests are skipped.
d = miniPipelineData()
quiet = function(x) suppressWarnings(suppressMessages(x))

test_that("runColocalisationModels runs end to end", {
    res = quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome"))
    expect_named(res, c("weighted", "unweighted", "modelWeighted",
        "modelUnweighted", "relativeWeighted", "relativeUnweighted",
        "modelRelativeWeighted", "modelRelativeUnweighted", "spicy"))
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

test_that("results do not depend on the number of cores", {
    skip_on_os("windows")
    expect_equal(
        quiet(bootstrapLineageConfidence(d$spe, d$tree, nBoot = 3, seed = 1,
            cores = 1)),
        quiet(bootstrapLineageConfidence(d$spe, d$tree, nBoot = 3, seed = 1,
            cores = 2)))
    expect_equal(
        quiet(crossfitLineageProbabilities(d$spe, d$tree, nFolds = 3, seed = 1, cores = 1)),
        quiet(crossfitLineageProbabilities(d$spe, d$tree, nFolds = 3, seed = 1, cores = 2)))
    expect_equal(
        quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome", cores = 1)),
        quiet(runColocalisationModels(d$spe, d$tree, d$conf, "outcome", cores = 2)))
})

test_that("the random forest does not depend on the number of cores", {
    skip_on_os("windows")
    expect_equal(
        quiet(bootstrapLineageConfidence(d$spe, d$tree, nBoot = 3, seed = 1, classifier = "rf", ntree = 50,
            cores = 1)),
        quiet(bootstrapLineageConfidence(d$spe, d$tree, nBoot = 3, seed = 1, classifier = "rf", ntree = 50,
            cores = 3)))
})

test_that("parallel workers run BLAS on one thread, and the session's setting is restored", {
    skip_on_os("windows")
    before = RhpcBLASctl::blas_get_num_procs()
    on.exit(RhpcBLASctl::blas_set_num_threads(before))
    RhpcBLASctl::blas_set_num_threads(2L)
    inWorkers = unlist(parLapply(1:4, function(i) RhpcBLASctl::blas_get_num_procs(), parallelParam(2)))
    expect_true(all(inWorkers == 1L))
    expect_equal(RhpcBLASctl::blas_get_num_procs(), 2L)
})
