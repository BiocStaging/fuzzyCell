test_that("runColocalisationModels returns the expected components (full data)", {
    skip_unless_full()
    d = fullTestData()
    res = suppressMessages(runColocalisationModels(d$spe, d$heartCellTypeTree,
        d$confList, outcomeCol = "resolution_broad",
        subject = "subject"))
    expect_named(res, c("weighted", "unweighted", "modelWeighted",
        "modelUnweighted", "relativeWeighted", "relativeUnweighted",
        "modelRelativeWeighted", "modelRelativeUnweighted", "spicy"))
    expect_true(all(c("from", "to", "estimate", "p.adj") %in%
        names(res$modelWeighted)))
})

test_that("weighted co-localisation equals unweighted for hard-labelled pairs (full data)", {
    skip_unless_full()
    d = fullTestData()
    res = suppressMessages(runColocalisationModels(d$spe, d$heartCellTypeTree,
        d$confList, outcomeCol = "resolution_broad", subject = "subject"))
    pick = function(m) m[m$from == "Endothelial" & m$to == "NK", c("estimate", "se", "p.value")]
    expect_equal(pick(res$modelWeighted), pick(res$modelUnweighted))
})
