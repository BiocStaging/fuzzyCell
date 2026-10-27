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
* The bootstrap resamples all cells' counts in one pass, grows a random forest on the workers, and accepts sparse
  counts.
* The bootstrap resamples the cells in blocks of about 4 million transcripts, so its memory no longer grows with the
  lineage's total count. Each transcript is drawn from gene j with probability count_j / total using the counts as
  given, so whole numbers stored with floating-point error are resampled correctly (whole-number counts give the
  same draws as before).
* Parallel loops (the bootstrap resamples, the random forest and the cross-fitting folds) run as in spicyR and
  lisaClust: the workers are forked once per call and run BLAS on one thread, and tasks are handed out in about
  eight chunks per worker. Results are unchanged.
* `correctedCounts = TRUE` in `bootstrapClassify()` and `bootstrapLineageConfidence()` accepts counts that are not
  whole numbers, such as counts corrected for ambient RNA, and resamples them as they are; they must still be on
  the scale of molecule counts. By default counts must be whole numbers, as before.
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
