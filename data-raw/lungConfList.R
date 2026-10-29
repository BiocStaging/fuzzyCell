## Bootstrap confidences for `lungSubset`, so the examples can skip the
## bootstrap: the same call as in the vignette.
## Re-create with `source("data-raw/lungConfList.R")` from the package root,
## after data-raw/lungSubset.R and data-raw/lungCellTypeTree.R.

devtools::load_all()
load(file.path("data", "lungSubset.rda"))
load(file.path("data", "lungCellTypeTree.rda"))

lungConfList = bootstrapLineageConfidence(lungSubset, lungCellTypeTree,
    nBoot = 20, classifier = "pcaLogistic", seed = 1)

save(lungConfList, file = file.path("data", "lungConfList.rda"),
    compress = "xz", version = 2)
