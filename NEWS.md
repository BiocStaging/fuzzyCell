# stBoot 0.99.1

* Co-localisation (`runColocalisationModels()`) is spicyR's cell-level test with the confidences as soft cell type
  weights (`spicyR::spicy(cellTypeWeights =)`), and relative co-localisation is Statial's Kontextual test with the
  same weights (`spicyR::kontextualEngine()`); images of a subject are combined, so several images per subject are
  allowed. New `effect` argument (`"allocation"` or `"count"`).
* Niches (`runNicheIdentification()`) cluster lisaClust's local curves with the weights (`lisaClust::lisa()`).
* `confidenceWeights()` turns a confidence list into the cells x labels weight matrix spicyR and lisaClust take.
* The default classifier is ridge multinomial logistic regression (`classifier = "logistic"`, glmnet);
  `"pcaLogistic"` fits it on the first 50 principal components.
* `priors = "observed"` weights each class by its share of the group's cells instead of equally.
* `crossfitLineageProbabilities()`: cross-fitted class probabilities, an alternative to the bootstrap.
* The bootstrap resamples all cells' counts in one pass, grows a random forest on the workers, and accepts sparse
  counts.

# stBoot 0.99.0

* Initial package, prepared for Bioconductor submission.
* Bootstrap classification confidence (`bootstrapLineageConfidence()`).
* Weighted cell type proportions, co-localisation and niche analyses.
* Example data `simSpe`, `simCellTypeTree`, `simConfList`.
