# fuzzyCell <img src="man/figures/logo.png" align="right" height="139" alt="fuzzyCell logo" />

Cell type labels in spatial single-cell data are uncertain, particularly between
closely related subtypes (for example CD4+ and CD8+ T cells), but most
downstream analyses treat them as known. fuzzyCell estimates a per-cell
classification **confidence** by bootstrapping (resampling each cell's
transcript counts and re-classifying it), and uses the confidences as **soft
cell type weights** in downstream spatial analyses:

- cell type proportions, among all cells and relative to the parent group;
- co-localisation of pairs of cell types, with spicyR's L-function, and relative
  co-localisation with Statial's Kontextual test (via spicyR);
- niche identification, from lisaClust's local indicators of spatial
  association.

Each analysis is run with the original hard labels and with the confidence
weights, so the two can be compared, and each can be tested against a
sample-level outcome. Cell types are organised with a user-supplied hierarchy
of lineages, groups and subtypes.

## Installation

### Development version

The development version is on [GitHub](https://github.com/shreyarajeshrao/fuzzyCell).
It depends on several Bioconductor packages (for example `SpatialExperiment`
and `Statial`), so install it with `BiocManager`, which installs those
dependencies from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("shreyarajeshrao/fuzzyCell")
```

Or, with `remotes`, which also needs the Bioconductor repositories set:

```r
# install.packages("remotes")
remotes::install_github("shreyarajeshrao/fuzzyCell",
    repos = BiocManager::repositories())
```

fuzzyCell requires R >= 4.5.0. To also build the vignette, add
`build_vignettes = TRUE`. This needs the [Quarto](https://quarto.org) command
line tool, the `quarto` R package, and the example data package
`SpatialDatasets`:

```r
BiocManager::install(c("quarto", "SpatialDatasets"))
BiocManager::install("shreyarajeshrao/fuzzyCell", build_vignettes = TRUE)
```

### Bioconductor release

fuzzyCell is being prepared for Bioconductor. Once it is accepted it will be
installed with:

```r
BiocManager::install("fuzzyCell")
```

## Quick start

```r
library(fuzzyCell)

# heart transplant Xenium biopsies (downloads the data on first use)
spe = loadAmancherla("pre")
data(heartCellTypeTree)

# per-cell confidence for each group of related subtypes
confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 50)

# weighted vs unweighted analyses, tested against the outcome
props = runProportionsPipeline(spe, heartCellTypeTree, confList,
    outcomeCol = "resolution_broad", subject = "subject")
coloc = runColocalisationModels(spe, heartCellTypeTree, confList,
    outcomeCol = "resolution_broad", subject = "subject")

# the weights themselves, to pass to spicyR or lisaClust directly
w = confidenceWeights(spe, heartCellTypeTree, confList)
```

The input is any `SummarizedExperiment` (such as a `SpatialExperiment`) with a
`"counts"` assay and `cellID`, `imageID`, `cellType`, `x` and `y` in its
`colData`.

See `vignette("fuzzyCell")` for a full walkthrough.

## Issues

Please report bugs and requests at
<https://github.com/shreyarajeshrao/fuzzyCell/issues>.
