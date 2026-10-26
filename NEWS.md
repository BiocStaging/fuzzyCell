# fuzzyCell 0.99.1

* Co-localisation (`runColocalisationModels()`) is spicyR's cell-level test with the confidences as soft cell type
  weights (`spicyR::spicy(cellWeights =)`), and relative co-localisation is Statial's Kontextual test with the
  same weights (`spicyR::kontextualEngine()`); images of a subject are combined, so several images per subject are
  allowed. New `effect` argument (`"allocation"` or `"count"`).
* Niches (`runNicheIdentification()`) cluster lisaClust's local curves with the weights (`lisaClust::lisa()`).
* `confidenceWeights()` turns a confidence list into the cells x labels weight matrix spicyR and lisaClust take.
* The default classifier is ridge multinomial logistic regression on the first 50 principal components
  (`classifier = "pcaLogistic"`, glmnet and irlba); `"logistic"` fits it on all genes.
* `priors = "observed"` weights each class by its share of the group's cells instead of equally.
* `crossfitLineageProbabilities()`: cross-fitted class probabilities, an alternative to the bootstrap.
* Checkpointing is removed: `bootstrapLineageConfidence()` no longer takes `checkpointDir`, and
  `loadCheckpointedConfList()` is gone. Running each group in a fresh R session, which keeps forked workers small,
  no longer needs checkpointing: it is the new argument `subprocess = TRUE`.
* The bootstrap resamples all cells' counts in one pass, grows a random forest on the workers, and accepts sparse
  counts.
* The bootstrap resamples the cells in blocks of about 4 million transcripts, so its memory no longer grows with the
  lineage's total count (the counts drawn are unchanged), and rounds counts stored with floating-point error.
* `classifier = "pcaLogistic"` fits the genes themselves when fewer than three vary, and the logistic classifiers
  accept a count matrix without gene names.
* `crossfitLineageProbabilities()` stratifies its folds by cell type, so every fold's classifier sees every child
  cell type; a child cell type needs at least 2 cells.
* When no pair has an eligible parent, the relative co-localisation tables are empty but keep their columns.

# fuzzyCell 0.99.0

* Initial package, prepared for Bioconductor submission.
* Bootstrap classification confidence (`bootstrapLineageConfidence()`).
* Weighted cell type proportions, co-localisation and niche analyses.
* Example data `simSpe`, `simCellTypeTree`, `simConfList`.
