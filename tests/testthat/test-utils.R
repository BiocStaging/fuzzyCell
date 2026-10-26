test_that("getCellData adds coordinates from a SpatialExperiment", {
    cells = getCellData(miniSpe())
    expect_equal(nrow(cells), 4)
    expect_equal(cells$x, c(1, 2, 3, 4))
    expect_equal(cells$y, c(4, 3, 2, 1))
    expect_true(all(c("cellID", "imageID", "cellType") %in% names(cells)))
})

test_that("generateComparisons returns all pair types", {
    data(heartCellTypeTree)
    confList = list("T cells" = NULL, "B lineage" = NULL)
    comp = fuzzyCell:::generateComparisons(heartCellTypeTree, confList)
    expect_named(comp, c("groupGroup", "groupLeaf", "leafLeaf"))
    nLeaves = nrow(buildCellTypeMap(heartCellTypeTree)$cellMap)
    expect_equal(nrow(comp$leafLeaf), nLeaves^2)
    expect_equal(nrow(comp$groupGroup), 4)
})
