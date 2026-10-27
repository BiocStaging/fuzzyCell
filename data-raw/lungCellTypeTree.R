## Cell type tree for the Vannan et al. (2025) pulmonary fibrosis lung Xenium
## data (SpatialDatasets::spe_Vannan_2025()), and so for `lungSubset`.
## Re-create with `source("data-raw/lungCellTypeTree.R")` from the package root.

lungCellTypeTree = list(
    "Immune" = list(
        "Myeloid" = list(
            "Macrophage lineage" = list(children = c("Alveolar Macrophages",
                "Interstitial Macrophages", "SPP1+ Macrophages",
                "Macrophages - IFN-activated", "Monocytes/MDMs")),
            "Dendritic cells" = list(children = c("cDCs", "Migratory DCs",
                "Langerhans cells", "pDCs")),
            "Granulocytes" = list(children = c("Mast", "Basophils", "Neutrophils")),
            "Proliferating Myeloid"),
        "Lymphoid" = list(
            "T cells" = list(children = c("CD4+ T-cells", "CD8+ T-cells", "Tregs",
                "Proliferating T-cells")),
            "NK lineage" = list(children = c("NK/NKT", "Proliferating NK/NKT")),
            "B lineage" = list(children = c("B cells", "Plasma",
                "Proliferating B cells")))),
    "Mesenchymal" = list(
        "Fibroblast lineage" = list(children = c("Alveolar FBs", "Adventitial FBs",
            "Activated Fibrotic FBs", "Inflammatory FBs", "Subpleural FBs",
            "Proliferating FBs", "Myofibroblasts")),
        "SMCs/Pericytes",
        "Mesothelial"),
    "Endothelial lineage" = list(children = c("Arteriole", "Venous", "Capillary",
        "Lymphatic")),
    "Epithelial lineage" = list(
        "Alveolar" = list(children = c("AT1", "AT2", "Transitional AT2",
            "Proliferating AT2")),
        "Airway" = list(children = c("Basal", "Secretory", "Goblet", "Multiciliated",
            "RASC", "PNEC", "KRT5-/KRT17+", "Proliferating Airway"))))

save(lungCellTypeTree, file = file.path("data", "lungCellTypeTree.rda"),
    compress = "xz", version = 2)
