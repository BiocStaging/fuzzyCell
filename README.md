# stBoot

Bootstrap classification confidence for spatial single-cell cell types, and
confidence-weighted downstream spatial analyses (proportions, co-localisation,
niches).

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("shreyarajeshrao/scBoot")
```

## Quick start

```r
library(stBoot)
spe = loadAmancherla("pre")      # downloads the heart transplant Xenium data
data(heartCellTypeTree)
confList = bootstrapLineageConfidence(spe, heartCellTypeTree, nBoot = 50)
res = runColocalisationModels(spe, heartCellTypeTree, confList, "resolution_broad")
```

See `vignette("stBoot")` for a full walkthrough.
