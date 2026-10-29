#' Cell type tree utilities
#'
#' Helpers for the cell type tree: a nested list describing how cell types are
#' organised into broad lineages and (optionally) groups of closely related
#' subtypes. [bootstrapLineageConfidence()] refines the cells within each group,
#' and the analysis functions compare cell types, groups and lineages.
#'
#' @details
#' A tree is a named list with one entry per broad lineage (these become the
#' `lineage` column). Each lineage is either a character vector of cell types,
#' or a list whose entries are:
#'
#' * a bare string: a cell type that is not in any group;
#' * a named `list(children = c(...))`: a **group** of closely related cell
#'   types, the unit that is bootstrap-refined. A group may contain only
#'   `children`;
#' * a named list of further entries: purely organisational nesting (a
#'   "wrapper") that is not a group itself, to any depth.
#'
#' A lineage can also itself be `list(children = c(...))`, in which case the
#' lineage name doubles as the group name. A cell type's group is its nearest
#' `list(children = ...)` ancestor, however many wrappers sit above it. In the
#' example below `CD8_Resting` is in the group `CD8+ T cells`, not `T cells`.
#'
#' Every cell type must appear exactly once, and lineage, wrapper and group
#' names must be unique and different from every cell type name.
#'
#' ```
#' list(
#'     "Immune" = list(
#'         "T cells" = list(
#'             "CD8+ T cells" = list(children = c("CD8_Resting", "CD8_Memory")),
#'             "CD4+ T cells" = list(children = c("CD4_naive", "CD4_Memory"))),
#'         "NK",
#'         "Mast"),
#'     "Endothelial lineage" = c("Endothelial", "Activated endothelial"))
#' ```
#'
#' `resolveTreeLeaves()` expands any node, however far up the tree, to the cell
#' types beneath it. It is used by relative co-localisation when it conditions
#' on a parent that is broader than a cell type's own group (see
#' `parentOverrides` in [runColocalisationModels()]), and `autoParentOverrides()` finds
#' those parents automatically.
#'
#' @param cellTypeTree A nested list describing the cell type hierarchy (see
#'   Details and [heartCellTypeTree]).
#' @param spe A `SummarizedExperiment` (e.g. a `SpatialExperiment`), or a
#'   `data.frame`/`DataFrame` of per-cell data.
#' @param cellType Name of the column of `colData(spe)` (or of the data frame)
#'   holding each cell's cell type label. An error is raised if it is not a
#'   column.
#' @param nodeName A single string: the name of any node in `cellTypeTree`, a
#'   lineage, an organisational wrapper, a group or a cell type.
#'
#' @return
#' * `buildCellTypeMap()`: a list with `cellMap` (a `data.frame` with columns
#'   `cellType`, `group` and `broad`, one row per cell type) and
#'   `groupsWithChildren` (names of the groups to be bootstrap-refined). An
#'   error is raised if the tree is not a valid tree as described in Details:
#'   it must be a non-empty named list, every top-level entry must be named,
#'   groups may only contain `children`, cell type names must be non-empty
#'   strings, and names must not be repeated or shared between cell types and
#'   lineages, wrappers or groups.
#' * `annotateCellTypeTree()`: `spe` with `group` and `lineage` columns added
#'   (a message reports any existing `group` or `lineage` column whose values
#'   are replaced). An error is raised if any cell's cell type is not in the
#'   tree (including `NA`).
#' * `resolveTreeLeaves()`: character vector of every cell type beneath
#'   `nodeName`.
#' * `autoParentOverrides()`: named character vector mapping each group that
#'   has a named ancestor to that ancestor, for use as `parentOverrides` in
#'   [runColocalisationModels()].
#'
#' @examples
#' data(lungCellTypeTree)
#' buildCellTypeMap(lungCellTypeTree)
#' resolveTreeLeaves(lungCellTypeTree, "Immune")
#' autoParentOverrides(lungCellTypeTree)
#'
#' data(lungSubset)
#' cells = annotateCellTypeTree(getCellData(lungSubset), lungCellTypeTree)
#' table(cells$group, useNA = "ifany")
#' @name cellTypeTree
#' @export
buildCellTypeMap = function(cellTypeTree) {
    if (!is.list(cellTypeTree) || length(cellTypeTree) == 0 || is.null(names(cellTypeTree))) {
        stop("cellTypeTree must be a non-empty named list, with one named entry ",
            "per lineage (see ?cellTypeTree)", call. = FALSE)
    }
    unnamed = which(is.na(names(cellTypeTree)) | !nzchar(names(cellTypeTree)))
    if (length(unnamed) > 0) {
        stop("Every top-level entry of cellTypeTree must be named (one per lineage), ",
            "but entry ", paste(unnamed, collapse = ", "), " ",
            if (length(unnamed) == 1) "is" else "are",
            " unnamed. Give it a name, for example list(Other = \"NK\").",
            call. = FALSE)
    }
    for (i in seq_along(cellTypeTree)) {
        .validateTreeNode(cellTypeTree[[i]], names(cellTypeTree)[i])
    }

    rows = list()

    walk = function(node, broad) {
        if (is.character(node)) {
            # a plain character vector: these leaves have no group of their own
            for (leaf in node) rows[[length(rows) + 1]] <<- list(broad = broad, group = NA_character_, child = leaf)
            return(invisible())
        }

        if (.isGroup(node)) {
            # a group passed in directly, with no organisational wrapper name
            # above it (e.g. a broad category that is immediately
            # list(children = ...)): the broad category name doubles as the
            # group name since there is nothing more specific
            for (leaf in node[["children"]]) rows[[length(rows) + 1]] <<- list(broad = broad, group = broad, child = leaf)
            return(invisible())
        }

        nms = names(node)
        if (is.null(nms)) nms = rep("", length(node))

        for (i in seq_along(node)) {
            entry = node[[i]]
            nm = nms[i]

            if (nm == "") {
                # bare leaf string sitting directly among siblings - no group of
                # its own
                rows[[length(rows) + 1]] <<- list(broad = broad, group = NA_character_, child = entry)
            } else if (.isGroup(entry)) {
                # a refinable group: nm is the group, its children are the
                # leaves
                for (leaf in entry[["children"]]) rows[[length(rows) + 1]] <<- list(broad = broad, group = nm, child = leaf)
            } else {
                # purely organisational nesting (e.g. "T cells" wrapping further
                # groups): recurse without introducing a group at this level
                walk(entry, broad)
            }
        }
    }

    for (i in seq_along(cellTypeTree)) walk(cellTypeTree[[i]], broad = names(cellTypeTree)[i])

    treeDF = dplyr::bind_rows(lapply(rows, as.data.frame, stringsAsFactors = FALSE))

    # every cell type must appear exactly once: annotateCellTypeTree() matches
    # cells to the tree by cell type, so a repeated name would silently put all
    # its cells in whichever group comes first
    repeated = unique(treeDF$child[duplicated(treeDF$child)])
    if (length(repeated) > 0) {
        describe = vapply(repeated, function(ct) {
            hits = treeDF[treeDF$child == ct, ]
            where = ifelse(is.na(hits$group),
                paste0("'", hits$broad, "' (no group)"), paste0("'", hits$group, "'"))
            sprintf("'%s' (under %s)", ct, paste(where, collapse = " and "))
        }, character(1))
        stop("Each cell type must appear exactly once in cellTypeTree, but ",
            paste(describe, collapse = "; "), " appears more than once.", call. = FALSE)
    }

    # lineage / wrapper / group names and cell type names share one namespace:
    # a label that is both would be read as the cell type everywhere, so the
    # lineage or group could never be analysed as such
    nodeNames = c(names(cellTypeTree), unlist(lapply(cellTypeTree, .collectNodeNames),
        use.names = FALSE))
    clash = intersect(nodeNames, treeDF$child)
    if (length(clash) > 0) {
        stop(paste0("'", clash, "'", collapse = ", "), " ",
            if (length(clash) == 1) "is" else "are",
            " used both as a lineage/group name and as a cell type in ",
            "cellTypeTree. Names must be unique across the whole tree.",
            call. = FALSE)
    }
    repeatedNodes = unique(nodeNames[duplicated(nodeNames)])
    if (length(repeatedNodes) > 0) {
        stop("Lineage and group names must be unique in cellTypeTree, but ",
            paste0("'", repeatedNodes, "'", collapse = ", "), " ",
            if (length(repeatedNodes) == 1) "is" else "are", " used more than once.",
            call. = FALSE)
    }

    cellMap = treeDF |> dplyr::transmute(cellType = child, group, broad)

    list(
        cellMap = cellMap,
        # groups that have children are the ones that need bootstrap refinement;
        # ungrouped leaf types (NK, Mast, Cardiomyocytes, ...) have group == NA
        # and are skipped.
        groupsWithChildren = unique(treeDF$group[!is.na(treeDF$group)])
    )
}

#' @rdname cellTypeTree
#' @export
annotateCellTypeTree = function(spe, cellTypeTree, cellType = "cellType") {
    checkColumn(spe, cellType, "cellType")
    treeInfo = buildCellTypeMap(cellTypeTree)
    cellMap = treeInfo$cellMap

    # every cell must belong to a cell type in the tree: otherwise it would
    # quietly get NA for group and lineage and drop out of the analyses
    labels = as.character(spe[[cellType]])
    notInTree = setdiff(unique(labels), cellMap$cellType)
    if (length(notInTree) > 0) {
        nCells = vapply(notInTree, function(ct) {
            if (is.na(ct)) sum(is.na(labels)) else sum(labels == ct, na.rm = TRUE)
        }, integer(1))
        stop(length(notInTree), " cell type", if (length(notInTree) > 1) "s",
            " (", sum(nCells), " cell", if (sum(nCells) != 1) "s", ") in the data ",
            if (length(notInTree) > 1) "are" else "is", " not in cellTypeTree: ",
            paste0("'", notInTree, "' (", nCells, " cell", ifelse(nCells == 1, "", "s"), ")",
                collapse = ", "),
            ". Add them to the tree, or remove these cells, before continuing.",
            call. = FALSE)
    }

    newValues = list(
        group = cellMap$group[match(spe[[cellType]], cellMap$cellType)],
        lineage = cellMap$broad[match(spe[[cellType]], cellMap$cellType)]
    )

    # tell the user if a column of theirs is being replaced with different
    # values (re-annotating already annotated cells is silent, as nothing changes)
    replaced = names(newValues)[vapply(names(newValues), function(col) {
        col %in% cellDataNames(spe) &&
            !isTRUE(all.equal(as.character(spe[[col]]), newValues[[col]]))
    }, logical(1))]
    if (length(replaced) > 0) {
        message("Replacing the existing ", paste0("'", replaced, "'", collapse = " and "),
            " column", if (length(replaced) > 1) "s", " in the data with ",
            if (length(replaced) > 1) "values" else "a column", " derived from cellTypeTree")
    }

    spe$group = newValues$group
    spe$lineage = newValues$lineage
    spe
}

# TRUE if `node` is a group: a list with an element named exactly "children"
# (exact matching, unlike `$`, which would also match e.g. "childrenX")
.isGroup = function(node) {
    is.list(node) && "children" %in% names(node)
}

# checks that the part of the tree below one lineage is well formed, and stops
# with a message saying where if not. `path` is the vector of names leading to
# `node`, used only for the message.
.validateTreeNode = function(node, path) {
    where = paste(path, collapse = " > ")

    if (is.character(node)) {
        if (length(node) == 0 || anyNA(node) || any(!nzchar(node))) {
            stop("'", where, "' has a missing or empty cell type name; every cell ",
                "type must be a non-empty string", call. = FALSE)
        }
        return(invisible())
    }
    if (!is.list(node)) {
        stop("'", where, "' must be a character vector of cell types or a list, ",
            "not ", if (is.null(node)) "NULL" else paste0("a ", class(node)[1]),
            call. = FALSE)
    }

    nms = names(node)
    if (is.null(nms)) nms = rep("", length(node))

    if (.isGroup(node)) {
        # a group is just a set of children, nothing nested inside it
        extra = setdiff(seq_along(node), which(nms == "children"))
        if (length(extra) > 0 || sum(nms == "children") > 1) {
            stop("The group '", where, "' may only contain `children = c(...)`",
                if (length(extra) > 0) paste0(", but it also has ", length(extra),
                    " other entr", if (length(extra) == 1) "y" else "ies",
                    " (", paste0("'", ifelse(nzchar(nms[extra]), nms[extra], paste0("[[", extra, "]]")),
                        "'", collapse = ", "), ")"),
                ". Put other cell types next to the group, not inside it.", call. = FALSE)
        }
        children = node[["children"]]
        if (!is.character(children) || length(children) == 0 || anyNA(children) || any(!nzchar(children))) {
            stop("`children` of '", where, "' must be a non-empty character vector of ",
                "cell type names, without missing or empty names", call. = FALSE)
        }
        return(invisible())
    }

    for (i in seq_along(node)) {
        entry = node[[i]]
        if (is.list(entry) && !nzchar(nms[i])) {
            stop("Entry ", i, " of '", where, "' is an unnamed list; nested ",
                "entries need a name", call. = FALSE)
        }
        .validateTreeNode(entry, c(path, if (nzchar(nms[i])) nms[i] else paste0("[[", i, "]]")))
    }
    invisible()
}

# names of every non-leaf node beneath `node` (organisational wrappers and
# list(children = ...) groups), at any depth. Leaves are character values and
# are not included.
.collectNodeNames = function(node) {
    if (!is.list(node)) {
        return(character(0))
    }
    nms = names(node)
    if (is.null(nms)) nms = rep("", length(node))
    out = character(0)
    for (i in seq_along(node)) {
        if (is.list(node[[i]])) {
            if (nzchar(nms[i])) out = c(out, nms[i])
            out = c(out, .collectNodeNames(node[[i]]))
        }
    }
    out
}

# collects every leaf beneath ANY node value - a bare character vector, a
# list(children = ...) group, or plain organisational nesting - regardless of
# depth
.collectAllLeaves = function(node) {
    if (is.character(node)) {
        return(node)
    }
    if (.isGroup(node)) {
        return(node[["children"]])
    }
    unlist(lapply(node, .collectAllLeaves), use.names = FALSE)
}

# finds the named node `target` anywhere inside `node` and collects every leaf
# beneath it - NULL if `target` doesn't appear anywhere in this subtree
.findAndCollectLeaves = function(node, target) {
    if (is.character(node)) {
        return(NULL)
    }
    nms = names(node)
    if (is.null(nms)) nms = rep("", length(node))
    for (i in seq_along(node)) {
        if (identical(nms[i], target)) {
            return(.collectAllLeaves(node[[i]]))
        }
        if (is.list(node[[i]])) {
            hit = .findAndCollectLeaves(node[[i]], target)
            if (!is.null(hit)) {
                return(hit)
            }
        }
    }
    NULL
}

#' @rdname cellTypeTree
#' @export
resolveTreeLeaves = function(cellTypeTree, nodeName) {
    if (!is.character(nodeName) || length(nodeName) != 1 || is.na(nodeName)) {
        stop("`nodeName` must be a single, non-missing string", call. = FALSE)
    }
    allLeaves = buildCellTypeMap(cellTypeTree)$cellMap$cellType
    if (nodeName %in% allLeaves) {
        return(nodeName)
    }

    for (broadName in names(cellTypeTree)) {
        if (identical(broadName, nodeName)) {
            return(unique(.collectAllLeaves(cellTypeTree[[broadName]])))
        }
        hit = .findAndCollectLeaves(cellTypeTree[[broadName]], nodeName)
        if (!is.null(hit)) {
            return(unique(hit))
        }
    }
    stop("No node named '", nodeName, "' found in cellTypeTree", call. = FALSE)
}

# finds the name of `target`'s immediate container - the node whose own named
# entries include `target` - anywhere inside `node`. `containerName` is what
# `node` itself is currently known as (the name that would be returned if
# `target` turned up as one of node's direct entries). NULL if `target` isn't
# found in this subtree.
.findImmediateAncestor = function(node, containerName, target) {
    if (is.character(node)) {
        return(NULL)
    }
    nms = names(node)
    if (is.null(nms)) nms = rep("", length(node))
    for (i in seq_along(node)) {
        if (identical(nms[i], target)) {
            return(containerName)
        }
        if (is.list(node[[i]])) {
            childContainerName = if (nzchar(nms[i])) nms[i] else containerName
            hit = .findImmediateAncestor(node[[i]], childContainerName, target)
            if (!is.null(hit)) {
                return(hit)
            }
        }
    }
    NULL
}

#' @rdname cellTypeTree
#' @export
autoParentOverrides = function(cellTypeTree) {
    groupNames = buildCellTypeMap(cellTypeTree)$groupsWithChildren

    parents = vapply(groupNames, function(g) {
        for (broadName in names(cellTypeTree)) {
            if (identical(broadName, g)) {
                # already top-level - no ancestor above it
                return(NA_character_)
            }
            hit = .findImmediateAncestor(cellTypeTree[[broadName]], broadName, g)
            if (!is.null(hit)) {
                return(hit)
            }
        }
        NA_character_
    }, character(1))

    parents[!is.na(parents)]
}
