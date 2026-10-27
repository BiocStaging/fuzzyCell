#' Heart transplant example data
#'
#' Downloads the Amancherla et al. (2025) heart transplant Xenium dataset from
#' the SpatialDatasets package (via ExperimentHub) and prepares it for fuzzyCell:
#' renames the columns fuzzyCell expects, keeps patients with cellular rejection,
#' and keeps one randomly chosen image per patient for each biopsy timing.
#' The first call downloads the data, which is large. With `"both"` every patient
#' has an image before and after transplant, so those images are not independent
#' and cannot be passed together to [runColocalisationModels()] with `subject = "subject"`.
#'
#' The prepared object has the columns fuzzyCell functions use: `cellID`,
#' `imageID`, `cellType`, `x`, `y` and `subject`, and the per-patient outcome
#' `resolution_broad` (`"persistence"` or `"resolved"`). The cell type
#' hierarchy for these cell types is in [heartCellTypeTree].
#'
#' @param timing Which biopsies to keep: `"pre"` or `"post"` transplant, or
#'   `"both"`.
#' @param seed Random seed used to choose one image per patient.
#'
#' @return A `SpatialExperiment`.
#'
#' @examples
#' \donttest{
#' spe = loadAmancherla("pre")
#' spe
#' table(spe$resolution_broad)
#' }
#' @export
loadAmancherla = function(timing = c("pre", "post", "both"), seed = 51773) {
    timing = match.arg(timing)
    if (!requireNamespace("SpatialDatasets", quietly = TRUE))
        stop("Package 'SpatialDatasets' is required; install it with ",
            "BiocManager::install(\"SpatialDatasets\")", call. = FALSE)

    spe = SpatialDatasets::spe_Amancherla_2025()
    spe$cellType = as.character(spe$ct_second_pass)
    spe$x = spe$adj_x_centroid
    spe$y = spe$adj_y_centroid
    spe$imageID = as.character(spe$Sample)
    spe$subject = as.character(spe$patient_id)
    spe$cellID = spe$cell_id

    spe = spe[, !is.na(spe$resolution_broad)]
    cellular = spe$patient_rejection_type == "cellular_rejection"

    # one image per patient: a random cell is drawn per patient and its image
    # kept. Pre is drawn before post so each timing's choice is reproducible.
    set.seed(seed)
    pickImages = function(speTiming) {
        colData(speTiming) |>
            as.data.frame() |>
            dplyr::group_by(subject) |>
            dplyr::slice_sample(n = 1) |>
            dplyr::pull(imageID)
    }
    spePre = spe[, cellular & spe$biopsy_timing == "pre"]
    spePre = spePre[, spePre$imageID %in% pickImages(spePre)]
    spePost = spe[, cellular & spe$biopsy_timing == "post"]
    spePost = spePost[, spePost$imageID %in% pickImages(spePost)]

    switch(timing, pre = spePre, post = spePost, both = BiocGenerics::cbind(spePre, spePost))
}

#' Cell type tree for the heart transplant data
#'
#' The cell type hierarchy used with [loadAmancherla()]: broad lineages, groups
#' of closely related subtypes (for example the T cell subtypes), and ungrouped
#' leaf cell types. See [cellTypeTree] for the format.
#'
#' @format A nested list; see [cellTypeTree].
#' @source Defined in `data-raw/heartCellTypeTree.R`.
#' @examples
#' data(heartCellTypeTree)
#' buildCellTypeMap(heartCellTypeTree)$groupsWithChildren
#' @name heartCellTypeTree
#' @docType data
#' @keywords datasets
"heartCellTypeTree"

#' A small cropped subset of the pulmonary fibrosis lung data
#'
#' Xenium data of lung biopsies from Vannan et al. (2025), available in full as
#' `SpatialDatasets::spe_Vannan_2025()`, cut down to ship with the package. It
#' has five control and five pulmonary fibrosis biopsies (the five with the
#' fewest cells of each), each cropped to the roughly 3,000 cells in a square
#' window at its centre, so each image is a window of a larger biopsy. Every
#' image comes from a different subject. Cell types with fewer than 20 cells in
#' the windows are removed. The hierarchy for these cell types is in
#' [lungCellTypeTree].
#'
#' @format A `SpatialExperiment` with the raw counts of 343 genes in 29,971
#'   cells (assay `counts`), the cell coordinates in `spatialCoords()`, and in
#'   `colData`:
#' \describe{
#'   \item{cellID}{A unique cell ID.}
#'   \item{imageID}{The image (the window of a biopsy); also `sample_id`.}
#'   \item{subject}{The patient.}
#'   \item{cellType}{The cell type label.}
#'   \item{disease_status}{`"Control"` or `"Disease"`.}
#' }
#' @source Built from `SpatialDatasets::spe_Vannan_2025()` by
#'   `data-raw/lungSubset.R`.
#' @examples
#' data(lungSubset)
#' lungSubset
#' table(lungSubset$imageID, lungSubset$disease_status)
#' @name lungSubset
#' @docType data
#' @keywords datasets
"lungSubset"

#' Cell type tree for the lung data
#'
#' The cell type hierarchy for [lungSubset] (and the full
#' `SpatialDatasets::spe_Vannan_2025()` data): broad lineages, groups of closely
#' related subtypes (for example the fibroblast subtypes), and ungrouped leaf
#' cell types. See [cellTypeTree] for the format.
#'
#' @format A nested list; see [cellTypeTree].
#' @source Defined in `data-raw/lungCellTypeTree.R`.
#' @examples
#' data(lungCellTypeTree)
#' buildCellTypeMap(lungCellTypeTree)$groupsWithChildren
#' @name lungCellTypeTree
#' @docType data
#' @keywords datasets
"lungCellTypeTree"
