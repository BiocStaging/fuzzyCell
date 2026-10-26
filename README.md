# fuzzyCell

Bootstrap classification confidence for spatial single-cell cell types, and
confidence-weighted downstream spatial analyses: proportions, co-localisation
(spicyR and Statial's Kontextual test) and niches (lisaClust), with the
confidences as soft cell type weights.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("shreyarajeshrao/fuzzyCell")
```

## Quick start

```r
library(fuzzyCell)
spe = loadAmancherla("pre")      # downloads the heart transplant Xenium data
data(heartCellTypeTree)
confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 50)
res = runColocalisationModels(spe, heartCellTypeTree, confList, "resolution_broad",
    subject = "subject")
w = confidenceWeights(spe, heartCellTypeTree, confList)   # the weights, for spicyR / lisaClust directly
```

See `vignette("fuzzyCell")` for a full walkthrough.
