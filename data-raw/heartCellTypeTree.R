## Cell type tree for the Amancherla et al. (2025) heart transplant Xenium
## data (SpatialDatasets::spe_Amancherla_2025()). Re-create with
## `source("data-raw/heartCellTypeTree.R")` from the package root.

heartCellTypeTree = list(
    "Immune" = list(
        "T cells" = list(children = c("CD8+ T cells", "CD4+ T cells",
            "Proliferating T cells", "Treg")),
        "B lineage" = list(children = c("B cells", "Plasma")),
        "DC lineage" = list(children = c("cDC1", "cDC2", "pDC",
            "Proliferating DCs", "mDC")),
        "Macrophage lineage" = list(children = c("Macrophages",
            "SPP1+ Macrophages")),
        "NK",
        "Mast"),
    "Mesenchymal" = list(
        "Fibroblast lineage" = list(children = c("Fibroblasts",
            "POSTN+ Fibroblasts", "Myofibroblasts")),
        "Vascular" = list(children = c("Pericytes", "Proliferating pericytes",
            "vSMCs")),
        "Cardiomyocytes",
        "Adipocytes"),
    "Endothelial lineage" = c("Endothelial", "Activated endothelial",
        "BMX+ Activated endothelial", "Proliferating endothelial",
        "Lymphatic endothelial"))

save(heartCellTypeTree, file = file.path("data", "heartCellTypeTree.rda"),
    compress = "xz", version = 2)
