#' Validate input parameters for ratio score computation
#'
#' Validates and processes input parameters for compute_ratio_score
#' function, checking counts matrix, cells, genes, weights, and method
#' specifications.
#'
#' @param counts Matrix, sparse matrix, or data frame containing gene
#'   expression counts.
#' @param cells Character vector of cell IDs (optional).
#' @param genes_1 Character vector of gene names for first gene set.
#' @param genes_2 Character vector of gene names for second gene set (optional).
#' @param method Character vector specifying scoring method.
#' @param weights_1 Named numeric vector of weights for genes_1
#'   (optional).
#' @param weights_2 Named numeric vector of weights for genes_2
#'   (optional).
#' @details Modifies variables in the parent frame using assign().
#' @return Returns 0 on success, stops with error if validation fails.
#' @importFrom methods is
validate_input <- function(
    counts,
    cells,
    genes_1,
    genes_2,
    method,
    weights_1,
    weights_2) {
  #-------------------------------------------------------------------------#
  # check counts
  #-------------------------------------------------------------------------#
  if (!is(counts, "matrix") &&
      !is(counts, "CsparseMatrix") &&
      !is(counts, "RsparseMatrix") &&
      !is(counts, "TsparseMatrix") &&
      !is(counts, "data.frame")) {
    stop("Counts is not a matrix or sparse matrix!")
  }
  #-------------------------------------------------------------------------#
  # check cells
  #-------------------------------------------------------------------------#
  if (!is.null(cells)) {
    c_1 <- cells[cells %in% colnames(counts)]
    if (length(c_1) == 0) {
      stop("No cells provided are present in count matrix")
    }
    if (length(c_1) < length(cells)) {
      warning(paste0(
        paste0(!cells %in% c_1, sep = " ", collapse = " "),
        " are not present in count matrix"))
    }
  } else {
    c_1 <- colnames(counts)
  }

  #-------------------------------------------------------------------------#
  # check genes
  #-------------------------------------------------------------------------#
  if (!is.null(genes_1)) {
    g_1 <- genes_1[genes_1 %in% rownames(counts)]
    if (length(g_1) == 0) {
      stop("No genes provided are present in count matrix")
    }
    if (length(g_1) < length(genes_1)) {
      warning(paste0(
        paste0(!genes_1 %in% g_1, sep = " ", collapse = " "),
        " are not present in count matrix"))
    }
  } else {
    stop("Please provide at least one gene group - parse to genes_1")
  }
  if (!is.null(genes_2)) {
    g_2 <- genes_2[genes_2 %in% rownames(counts)]
    if (length(g_2) == 0) {
      stop("No genes provided are present in count matrix")
    }
    if (length(g_2) < length(genes_2)) {
      warning(paste0(
        paste0(!genes_2 %in% g_2, sep = " ", collapse = " "),
        " are not present in count matrix"))
    }
  } else {
    g_2 <- g_1
  }
  #-------------------------------------------------------------------------#
  # check weights
  #-------------------------------------------------------------------------#
  if (!is.null(weights_1)) {
    w_1 <- weights_1[names(weights_1) %in% g_1]
    if (length(w_1) == 0) {
      stop("Weight names do not match gene provided")
    }
    if (length(w_1) < length(weights_1)) {
      warning(paste0(
        paste0(!names(weights_1) %in% g_1, sep = " ", collapse = " "),
        " are not present gene list"))
    }
    weights_1_pad <- rep(1, length(g_1))
    names(weights_1_pad) <- g_1
    weights_1_pad[names(w_1)] <- w_1
  } else {
    weights_1_pad <- rep(1, length(g_1))
  }
  if (!is.null(weights_2)) {
    w_2 <- weights_2[names(weights_2) %in% g_2]
    if (length(w_2) == 0) {
      stop("Weight names do not match gene provided")
    }
    if (length(w_2) < length(weights_2)) {
      warning(paste0(
        paste0(!names(weights_2) %in% g_2, sep = " ", collapse = " "),
        " are not present gene list"))
    }
    weights_2_pad <- rep(1, length(g_2))
    names(weights_2_pad) <- g_2
    weights_2_pad[names(w_2)] <- w_2
  } else {
    weights_2_pad <- rep(1, length(g_2))
  }
  #-------------------------------------------------------------------------#
  # check methods
  #-------------------------------------------------------------------------#
  if (!method[1L] %in% c("mean_ratio", "mean_sub", "sum_ratio",
      "sum_sub")) {
    stop("Method provided not in available methods: ",
      "mean_ratio, mean_sub, sum_ratio, sum_sub")
  }
  assign("counts", counts, envir = parent.frame())
  assign("cells", c_1, envir = parent.frame())
  assign("genes_1", g_1, envir = parent.frame())
  assign("genes_2", g_2, envir = parent.frame())
  assign("weights_1", weights_1_pad, envir = parent.frame())
  assign("weights_2", weights_2_pad, envir = parent.frame())
  assign("method", method[1L], envir = parent.frame())
  return(0)
}

#' Check and format grouping structure
#'
#' Validates and formats group assignments into a standardized data
#'   frame format.
#'
#' @param group_id Data frame, matrix, or named vector containing group
#'   assignments.
#' @param group_name Character string specifying column name if group_id
#'   is a data frame.
#' @importFrom methods is
#' @return Data frame with columns "cell_id" and "group_id".
#' @importFrom methods is
check_grouping <- function(group_id, group_name) {
  if (is(group_id, "data.frame") || is(group_id,"matrix")) {
    if (is.null(group_name)) {
      stop("No group name specified for group_id data frame")
    }
    if (!group_name %in% colnames(group_id)) {
      stop("Grouping name is not present in grouping data frame.")
    }
    group_id <- data.frame(
      "cell_id" = rownames(group_id),
      "group_id" = group_id[, group_name])
  } else if (is(group_id, "vector")) {
    if (is.null(names(group_id))) {
      stop("No cell ids in group_id vector - add cell ids as names")
    }
    group_id <- data.frame(
      "cell_id" = names(group_id),
      "group_id" = group_id)
  } else {
    stop("Unsupported group_id format - use ",
      "data.frame/matrix/vector")
  }
  return(group_id)
}

#' Check and filter gene set against count matrix
#'
#' Validates that genes in the gene set are present in the count matrix
#' and returns filtered gene list.
#'
#' @param counts Matrix or sparse matrix with genes as row names.
#' @param genes Character vector of gene names to check.
#' @return Character vector of genes present in counts matrix, or NULL
#'   if input is NULL.
check_gene_set <- function(counts, genes){
  if (!is.null(genes)) {
    gs <- genes[genes %in% rownames(counts)]
    if (length(gs) == 0) {
      stop("No genes provided are present in count matrix")
    }
    if (length(gs) < length(genes)) {
      missing <- genes[!genes %in% gs]
      warning(paste0(
        paste0(missing, sep = " ", collapse = " "),
        " are not present in count matrix"))
    }
    return(gs)
  } else {
    return(genes)
  }
}
