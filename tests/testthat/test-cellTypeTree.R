data(heartCellTypeTree)

test_that("buildCellTypeMap flattens the tree", {
    m = buildCellTypeMap(heartCellTypeTree)
    expect_true(all(c("CD8+ T cells", "Macrophages", "NK", "Endothelial") %in%
        m$cellMap$cellType))
    expect_true("T cells" %in% m$groupsWithChildren)
    expect_false("NK" %in% m$groupsWithChildren)
    expect_true(is.na(m$cellMap$group[m$cellMap$cellType == "NK"]))
    expect_equal(m$cellMap$group[m$cellMap$cellType == "CD8+ T cells"],
        "T cells")
    expect_equal(m$cellMap$broad[m$cellMap$cellType == "CD8+ T cells"],
        "Immune")
})

test_that("annotateCellTypeTree adds group and lineage", {
    cells = data.frame(cellType = c("CD8+ T cells", "NK", "Endothelial"))
    cells = annotateCellTypeTree(cells, heartCellTypeTree)
    expect_equal(cells$group, c("T cells", NA, NA))
    expect_equal(cells$lineage, c("Immune", "Immune", "Endothelial lineage"))
})

test_that("annotateCellTypeTree errors on cell types that are not in the tree", {
    cells = data.frame(cellType = c("Treg", "TYPO", "TYPO", "NK", NA))
    expect_error(annotateCellTypeTree(cells, heartCellTypeTree),
        "2 cell types \\(3 cells\\) in the data are not in cellTypeTree: 'TYPO' \\(2 cells\\), 'NA' \\(1 cell\\)")
    expect_error(
        annotateCellTypeTree(data.frame(cellType = c("Treg", "tReg")),
            heartCellTypeTree),
        "1 cell type \\(1 cell\\) in the data is not in cellTypeTree: 'tReg' \\(1 cell\\)")
    # factors are fine
    ok = annotateCellTypeTree(data.frame(cellType = factor(c("Treg", "NK"))),
        heartCellTypeTree)
    expect_equal(ok$group, c("T cells", NA))
})

test_that("resolveTreeLeaves finds leaves at any depth", {
    mesenchymal = resolveTreeLeaves(heartCellTypeTree, "Mesenchymal")
    expect_true(all(c("Fibroblasts", "Pericytes", "Cardiomyocytes") %in%
        mesenchymal))
    expect_equal(resolveTreeLeaves(heartCellTypeTree, "Treg"), "Treg")
    expect_error(resolveTreeLeaves(heartCellTypeTree, "nonexistent"),
        "No node named")
})

test_that("autoParentOverrides maps groups to named ancestors", {
    ov = autoParentOverrides(heartCellTypeTree)
    expect_equal(unname(ov["T cells"]), "Immune")
    expect_equal(unname(ov["Fibroblast lineage"]), "Mesenchymal")
})

test_that("buildCellTypeMap rejects cell types that appear more than once", {
    dup = list(Immune = list(
        "T cells" = list(children = c("CD8", "Shared")),
        "B lineage" = list(children = c("B", "Shared"))))
    expect_error(buildCellTypeMap(dup),
        "'Shared' \\(under 'T cells' and 'B lineage'\\) appears more than once")

    # an ungrouped leaf repeated inside a group, and a repeat within one group
    dup2 = list(Immune = list("NK", "T cells" = list(children = c("NK", "x", "x"))))
    expect_error(buildCellTypeMap(dup2), "'NK'.*'x'.*appears more than once")

    # annotateCellTypeTree and resolveTreeLeaves inherit the check
    expect_error(annotateCellTypeTree(data.frame(cellType = "B"), dup),
        "exactly once")
    expect_error(resolveTreeLeaves(dup, "Immune"), "exactly once")

    # the shipped tree and a plain character-vector lineage are fine
    expect_no_error(buildCellTypeMap(heartCellTypeTree))
    expect_no_error(buildCellTypeMap(list(E = c("a", "b"))))
})

test_that("buildCellTypeMap rejects lineage/group names that clash with cell types", {
    # group named like a leaf elsewhere in the tree
    clash = list(Immune = list("Macrophage" = list(children = c("M1", "M2"))),
        Other = c("Macrophage"))
    expect_error(buildCellTypeMap(clash),
        "'Macrophage' is used both as a lineage/group name and as a cell type")

    # organisational wrapper named like a leaf
    wrapper = list(Mesenchymal = list("Vascular" = list(
        "Pericytes", "x" = list(children = c("a", "b")))),
        Other = c("Vascular"))
    expect_error(buildCellTypeMap(wrapper), "'Vascular' is used both")

    # a lineage named like one of its own leaves
    lineage = list(Endothelial = c("Endothelial", "Lymphatic"))
    expect_error(buildCellTypeMap(lineage), "'Endothelial' is used both")

    # the same group name used twice
    twice = list(A = list("G" = list(children = c("a", "b"))),
        B = list("G" = list(children = c("c", "d"))))
    expect_error(buildCellTypeMap(twice),
        "Lineage and group names must be unique.*'G'")

    # a lineage given directly as list(children = ...) is fine
    expect_no_error(buildCellTypeMap(list(Epi = list(children = c("a", "b")))))
    expect_no_error(buildCellTypeMap(heartCellTypeTree))
})

test_that("buildCellTypeMap rejects trees that are not non-empty named lists", {
    msg = "cellTypeTree must be a non-empty named list"
    expect_error(buildCellTypeMap(list()), msg)
    expect_error(buildCellTypeMap(c("A", "B")), msg)
    expect_error(buildCellTypeMap(NULL), msg)
    expect_error(buildCellTypeMap("Immune"), msg)
    expect_error(buildCellTypeMap(list("A", "B")), msg)
    expect_no_error(buildCellTypeMap(list(Immune = c("A", "B"))))
})

test_that("annotateCellTypeTree takes a custom cell type column and checks it", {
    cells = data.frame(ct = c("CD8+ T cells", "NK"), other = 1:2)
    ann = annotateCellTypeTree(cells, heartCellTypeTree, cellType = "ct")
    expect_equal(ann$group, c("T cells", NA))
    expect_equal(ann$lineage, c("Immune", "Immune"))

    expect_error(annotateCellTypeTree(cells, heartCellTypeTree),
        "`cellType = \"cellType\"` is not a column of the data. Available columns: ct, other")
    expect_error(annotateCellTypeTree(cells, heartCellTypeTree, cellType = 5),
        "`cellType` must be a single string naming a column")

    # also works on a SummarizedExperiment, whose colData is checked
    se = SummarizedExperiment::SummarizedExperiment(
        list(counts = matrix(1, 1, 2)),
        colData = S4Vectors::DataFrame(ct = c("Treg", "NK")))
    expect_equal(annotateCellTypeTree(se, heartCellTypeTree,
        cellType = "ct")$group, c("T cells", NA))
    expect_error(annotateCellTypeTree(se, heartCellTypeTree),
        "is not a column of colData\\(spe\\)")
})

test_that("annotateCellTypeTree reports replaced group/lineage columns", {
    cells = data.frame(cellType = c("Treg", "NK"), lineage = c("mine", "mine"))
    expect_message(ann <- annotateCellTypeTree(cells, heartCellTypeTree),
        "Replacing the existing 'lineage' column in the data")
    expect_equal(ann$lineage, c("Immune", "Immune"))

    # re-annotating already annotated cells changes nothing, so stays silent
    expect_no_message(annotateCellTypeTree(ann, heartCellTypeTree))

    both = data.frame(cellType = "Treg", group = "g", lineage = "l")
    expect_message(annotateCellTypeTree(both, heartCellTypeTree),
        "'group' and 'lineage' columns")
})

test_that("buildCellTypeMap requires every top-level entry to be named", {
    expect_error(buildCellTypeMap(list(G = list(children = c("a", "b")), "C")),
        "entry 2 is unnamed. Give it a name")
    expect_error(buildCellTypeMap(list(A = "a", "B", "C")),
        "entry 2, 3 are unnamed")
    expect_no_error(buildCellTypeMap(list(G = list(children = c("a", "b")),
        Other = "C")))
})

test_that("a group may only contain `children`", {
    extra = list(Imm = list("G" = list(children = c("a", "b"),
        "W" = list(children = c("c", "d")))))
    expect_error(buildCellTypeMap(extra),
        "group 'Imm > G' may only contain `children = c\\(...\\)`, but it also has 1 other entry \\('W'\\)")
    unnamedExtra = list(Imm = list("G" = list(children = c("a", "b"), "x")))
    expect_error(buildCellTypeMap(unnamedExtra), "'\\[\\[2\\]\\]'")
    # a lineage that is itself a group is validated the same way
    expect_error(buildCellTypeMap(list(Epi = list(children = c("a"), "x"))),
        "group 'Epi' may only contain")
})

test_that("cell type names must be non-empty strings", {
    expect_error(buildCellTypeMap(list(Imm = list("G" = list(children = c(1, 2))))),
        "`children` of 'Imm > G' must be a non-empty character vector")
    expect_error(buildCellTypeMap(list(Imm = list("G" = list(children = character(0))))),
        "must be a non-empty character vector")
    expect_error(buildCellTypeMap(list(Imm = c("a", ""))),
        "'Imm' has a missing or empty cell type name")
    expect_error(buildCellTypeMap(list(Imm = c("a", NA))),
        "'Imm' has a missing or empty cell type name")
    expect_error(buildCellTypeMap(list(Imm = list(W = list("a", 5)))),
        "'Imm > W > \\[\\[2\\]\\]' must be a character vector of cell types or a list, not a numeric")
    expect_error(buildCellTypeMap(list(Imm = NULL)), "not NULL")
})

test_that("unnamed nested lists are rejected", {
    expect_error(buildCellTypeMap(list(Imm = list("a", list(children = c("b", "c"))))),
        "Entry 2 of 'Imm' is an unnamed list")
})

test_that("only an exact `children` element makes a group", {
    # a wrapper whose name merely starts with "children" is not a group
    tree = list(Imm = list(childrenX = c("a", "b")))
    m = buildCellTypeMap(tree)
    expect_true(all(is.na(m$cellMap$group)))
    expect_length(m$groupsWithChildren, 0)
})

test_that("resolveTreeLeaves needs a single non-missing string", {
    msg = "`nodeName` must be a single, non-missing string"
    expect_error(resolveTreeLeaves(heartCellTypeTree, c("Treg", "NK")), msg)
    expect_error(resolveTreeLeaves(heartCellTypeTree, NA_character_), msg)
    expect_error(resolveTreeLeaves(heartCellTypeTree, 5), msg)
    expect_error(resolveTreeLeaves(heartCellTypeTree, character(0)), msg)
})
