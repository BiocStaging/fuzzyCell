test_that("getCellData adds coordinates from a SpatialExperiment", {
    cells = getCellData(miniSpe())
    expect_equal(nrow(cells), 4)
    expect_equal(cells$x, c(1, 2, 3, 4))
    expect_equal(cells$y, c(4, 3, 2, 1))
    expect_true(all(c("cellID", "imageID", "cellType") %in% names(cells)))
})

test_that("Lfunction is zero under complete spatial randomness", {
    n1 = 10
    n2 = 100
    area = 1e4
    r = 50
    numerator = n1 * (n2 / area) * pi * r^2
    expect_equal(scBoot:::Lfunction(numerator, r, area, n1, n2), 0)
})

test_that("borderEdge is 1 in the interior and < 1 near the edge", {
    X = spatstat.geom::ppp(x = c(100, 1), y = c(100, 1),
        window = spatstat.geom::owin(c(0, 200), c(0, 200)))
    e = scBoot:::borderEdge(X, 20)
    expect_equal(e[1], 1)
    expect_lt(e[2], 1)
})

test_that("generateComparisons returns all pair types", {
    data(heartCellTypeTree)
    confList = list("T cells" = NULL, "B lineage" = NULL)
    comp = scBoot:::generateComparisons(heartCellTypeTree, confList)
    expect_named(comp, c("groupGroup", "groupLeaf", "leafLeaf"))
    nLeaves = nrow(buildCellTypeMap(heartCellTypeTree)$cellMap)
    expect_equal(nrow(comp$leafLeaf), nLeaves^2)
    expect_equal(nrow(comp$groupGroup), 4)
})
