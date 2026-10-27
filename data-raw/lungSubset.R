## A small subset of the Vannan et al. (2025) pulmonary fibrosis lung Xenium
## data (SpatialDatasets::spe_Vannan_2025()), small enough to ship with the
## package: five control and five disease biopsies (the five with the fewest
## cells of each), each cropped to the ~3,000 cells in a square window at its
## centre, with the raw counts only. Cell types with fewer than 20 cells in the
## crops are removed, as the classifier needs enough cells of each subtype.
## Re-create with `source("data-raw/lungSubset.R")` from the package root.

library(SpatialExperiment)
library(Matrix)

spe = SpatialDatasets::spe_Vannan_2025()

images = c(
    "VUHD038", "VUHD095", "VUHD113", "VUHD090", "VUHD069",            # control
    "VUILD142MA", "TILD049MA", "TILD167LA", "VUILD105MA2", "TILD113LA") # disease
spe = spe[, as.character(spe$sample) %in% images]

# the ~3,000 cells nearest the centre of each image, in a square window
nPerImage = 3000
xy = spatialCoords(spe)
keep = unlist(lapply(split(seq_len(ncol(spe)), as.character(spe$sample)), function(i) {
    centre = apply(xy[i, , drop = FALSE], 2, stats::median)
    distance = pmax(abs(xy[i, 1] - centre[1]), abs(xy[i, 2] - centre[2]))
    i[distance <= sort(distance)[nPerImage]]
}))
spe = spe[, sort(keep)]

counts = as(assay(spe, "counts"), "CsparseMatrix")
colnames(counts) = spe$full_cell_id
lungSubset = SpatialExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
        cellID = spe$full_cell_id,
        imageID = as.character(spe$sample),
        subject = as.character(spe$patient),
        cellType = as.character(spe$final_CT),
        disease_status = as.character(spe$disease_status),
        row.names = spe$full_cell_id),
    spatialCoords = spatialCoords(spe),
    sample_id = as.character(spe$sample))

rare = names(which(table(lungSubset$cellType) < 20))
message("Removing cell types with fewer than 20 cells: ", paste(rare, collapse = ", "))
lungSubset = lungSubset[, !lungSubset$cellType %in% rare]

save(lungSubset, file = file.path("data", "lungSubset.rda"), compress = "xz", version = 2)
