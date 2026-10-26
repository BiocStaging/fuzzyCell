## Heavy tests use the real heart transplant data (SpatialDatasets) and are
## skipped unless FUZZYCELL_FULL_TESTS is set, e.g. Sys.setenv(FUZZYCELL_FULL_TESTS =
## "true"). They are far too slow for Bioconductor's build system.
skip_unless_full = function() {
    skip_on_bioc()
    skip_if(!nzchar(Sys.getenv("FUZZYCELL_FULL_TESTS")),
        "set FUZZYCELL_FULL_TESTS=true to run the full-data tests")
    skip_if_not_installed("SpatialDatasets")
}

## load the example data and confidence scores once per session
.fullData = new.env()
fullTestData = function() {
    if (is.null(.fullData$spe)) {
        data(heartCellTypeTree, package = "fuzzyCell", envir = .fullData)
        .fullData$spe = suppressMessages(loadAmancherla("pre"))
        .fullData$confList = suppressMessages(suppressWarnings(
            bootstrapLineageConfidence(.fullData$spe,
                .fullData$heartCellTypeTree, nBoot = 5, seed = 1)))
    }
    .fullData
}

## tiny hand-built SpatialExperiment with coordinates in spatialCoords()
miniSpe = function() {
    counts = matrix(1:12, nrow = 3,
        dimnames = list(paste0("g", 1:3), paste0("c", 1:4)))
    SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(cellID = paste0("c", 1:4),
            imageID = c("a", "a", "b", "b"),
            cellType = c("T", "T", "B", "B")),
        spatialCoords = cbind(x = c(1, 2, 3, 4), y = c(4, 3, 2, 1)))
}

## small hand-made dataset (6 images, 3 cell types, one group) with
## bootstrap confidences, used to exercise the pipelines end to end
miniPipelineData = function() {
    set.seed(3)
    nImg = 6
    nCell = 70
    cd = do.call(rbind, lapply(seq_len(nImg), function(i) data.frame(
        cellID = paste0("i", i, "_", seq_len(nCell)),
        imageID = paste0("img", i),
        cellType = sample(c("A", "B", "C"), nCell, TRUE),
        outcome = rep(c("x", "y"), each = 3)[i],
        x = runif(nCell, 0, 300), y = runif(nCell, 0, 300))))
    mu = matrix(1, 8, nrow(cd))
    mu[1:3, cd$cellType == "A"] = 5
    mu[4:6, cd$cellType == "B"] = 5
    counts = matrix(rpois(length(mu), mu), 8, nrow(cd),
        dimnames = list(paste0("g", 1:8), cd$cellID))
    spe = SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        colData = cd[, c("cellID", "imageID", "cellType", "outcome")],
        spatialCoords = as.matrix(cd[, c("x", "y")]))
    tree = list(G = list(children = c("A", "B")), Other = "C")
    conf = suppressWarnings(suppressMessages(
        bootstrapLineageConfidence(spe, tree, nBoot = 3, seed = 1)))
    list(spe = spe, tree = tree, conf = conf)
}

## a copy of a SpatialExperiment whose cell ID, cell type, image and coordinate
## columns have non-standard names (barcode, ct, sample, px, py)
renamedSpe = function(spe) {
    spe$barcode = spe$cellID
    spe$ct = spe$cellType
    spe$sample = spe$imageID
    spe$px = SpatialExperiment::spatialCoords(spe)[, 1]
    spe$py = SpatialExperiment::spatialCoords(spe)[, 2]
    spe$cellID = NULL
    spe$cellType = NULL
    spe$imageID = NULL
    spe
}

## wrap a per-cell data frame into a SummarizedExperiment (one column per cell),
## for the functions that take an `spe`
asSpe = function(cells) {
    ids = paste0("c", seq_len(nrow(cells)))
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(0, 1, nrow(cells), dimnames = list("g", ids))),
        colData = S4Vectors::DataFrame(cells, row.names = ids, check.names = FALSE))
}
