test_that("runColocalisationModels returns the expected components (full data)", {
    skip_unless_full()
    d = fullTestData()
    res = suppressMessages(runColocalisationModels(d$spe, d$heartCellTypeTree,
        d$confList, outcomeCol = "resolution_broad",
        subject = "subject"))
    expect_named(res, c("weighted", "unweighted", "modelWeighted",
        "modelUnweighted", "relativeWeighted", "relativeUnweighted",
        "modelRelativeWeighted", "modelRelativeUnweighted"))
    expect_true(all(c("from", "to", "estimate", "p.adj") %in%
        names(res$modelWeighted)))
})

test_that("weighted co-localisation equals unweighted for hard-labelled pairs (full data)", {
    skip_unless_full()
    d = fullTestData()
    cells = annotateCellTypeTree(getCellData(d$spe), d$heartCellTypeTree)
    cells$group[is.na(cells$group)] = cells$cellType[is.na(cells$group)]
    pairs = data.frame(from = "Endothelial", to = "NK")
    coloc = stBoot:::calcLfuncPerImage(cells, d$confList, pairs, r = 50)
    expect_equal(coloc$weighted$colocalisation, coloc$unweighted$colocalisation)
})
