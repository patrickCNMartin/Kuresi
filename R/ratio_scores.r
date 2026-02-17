#' Compute ratio scores by group
#'
#' Computes ratio scores for each group in a grouping structure, calculating
#' competitive fitness metrics between two gene sets.
#'
#' @param counts A matrix or sparse matrix containing gene expression
#'   counts.
#' @param groups A data frame, matrix, or named vector containing group
#'   assignments.
#' @param group_name Character string specifying the column name in
#'   groups if groups is a data frame.
#' @param genes_1 Character vector of gene names for the first gene set.
#' @param genes_2 Character vector of gene names for the second gene set
#'   (optional).
#' @param method Character vector specifying the scoring method:
#'   "mean_ratio", "mean_sub", "sum_ratio", or "sum_sub".
#' @param weights_1 Named numeric vector of weights for genes_1
#'   (optional).
#' @param weights_2 Named numeric vector of weights for genes_2
#'   (optional).
#' @param scale Logical whether to scale counts (default: FALSE).
#' @param center Logical whether to center scores (default: FALSE).
#' @param rank logical wether to use score rank (default: FALSE).
#' @param add_name Character string to add to output column name
#'   (optional).
#' @param verbose Logical whether to print progress messages (default: TRUE).
#' @return A data frame with group assignments and computed scores.
#' @export
compute_ratio_bygroup <- function(
    counts,
    groups,
    group_name = NULL,
    genes_1 = NULL,
    genes_2 = NULL,
    method = c("mean_ratio", "mean_sub", "sum_ratio", "sum_sub"),
    weights_1 = NULL,
    weights_2 = NULL,
    scale = FALSE,
    center = FALSE,
    raw = FALSE,
    rank = FALSE,
    add_name = NULL,
    verbose = TRUE) {
  #-------------------------------------------------------------------------#
  # Get groupings
  #-------------------------------------------------------------------------#
  group_id <- check_grouping(groups, group_name)
  group_id <- split(group_id, group_id[, "group_id"])

  for (g in seq_along(group_id)) {
    scores <- compute_ratio_score(
      counts,
      cells = group_id[[g]]$cell_id,
      genes_1 = genes_1,
      genes_2 = genes_2,
      method = method[1L],
      weights_1 = weights_1,
      weights_2 = weights_2,
      scale = scale,
      center = FALSE,
      raw = TRUE,
      collapse = TRUE,
      verbose = verbose)
    rownames(group_id[[g]]) <- group_id[[g]]$cell_id
    group_id[[g]]$score <- scores
  }
  #-------------------------------------------------------------------------#
  # Compute global metrics after computing per grouping scores
  #-------------------------------------------------------------------------#
  group_id <- do.call("rbind", group_id)
  if (center && grepl("mean", method[1L])) {
    group_id$score <- group_id$score - 1
  } else if (center && grepl("sub", method[1L])) {
    group_id$score <- group_id$score - min(group_id$score)
  } else if (!center && !raw) {
    group_id$score <- (group_id$score - min(group_id$score)) / (
      max(group_id$score) - min(group_id$score))
  }
  if (rank) {
    ord <- sort(unique(group_id$score), decreasing = TRUE)
    group_id$score <- match(group_id$score, ord)
  }
  group_id <- group_id[match(rownames(groups), group_id$cell_id), "score"]
  groups$score <- group_id
  colnames(groups) <- gsub("score", method[1L], colnames(groups))
  return(groups)
}




#' Compute ratio score between two gene sets
#'
#' Computes competitive fitness scores based on the ratio or difference
#' between two gene sets for individual cells or groups of cells.
#'
#' @param counts A matrix or sparse matrix containing gene expression
#'   counts.
#' @param cells Character vector of cell IDs to include (optional, uses
#'   all cells if NULL).
#' @param genes_1 Character vector of gene names for the first gene set.
#' @param genes_2 Character vector of gene names for the second gene set
#'   (optional, uses genes_1 if NULL).
#' @param method Character vector specifying the scoring method:
#'   "mean_ratio", "mean_sub", "sum_ratio", or "sum_sub".
#' @param weights_1 Named numeric vector of weights for genes_1
#'   (optional).
#' @param weights_2 Named numeric vector of weights for genes_2
#'   (optional).
#' @param scale Logical whether to scale counts before computing scores
#'   (default: FALSE).
#' @param center Logical whether to center scores (default: FALSE).
#' @param raw Logical whether to return raw uncentered scores (default: FALSE).
#' @param collapse Logical whether to return a single aggregated value
#'   instead of per-cell scores (default: FALSE).
#' @param verbose Logical whether to print progress messages
#'   (default: TRUE).
#' @return Numeric vector of scores (one per cell) or a single numeric
#'   value if collapse is TRUE.
#' @export
compute_ratio_score <- function(
    counts,
    cells = NULL,
    genes_1 = NULL,
    genes_2 = NULL,
    method = c("mean_ratio", "mean_sub", "sum_ratio", "sum_sub"),
    weights_1 = NULL,
    weights_2 = NULL,
    scale = FALSE,
    center = FALSE,
    raw = FALSE,
    collapse = FALSE,
    verbose = TRUE) {
  #-------------------------------------------------------------------------#
  # checks
  #-------------------------------------------------------------------------#
  validate_input(
    counts,
    cells,
    genes_1,
    genes_2,
    method,
    weights_1,
    weights_2
  )
  #-------------------------------------------------------------------------#
  # subset counts
  #-------------------------------------------------------------------------#
  counts_1 <- counts[genes_1, cells] * weights_1
  counts_2 <- counts[genes_2, cells] * weights_2
  if (scale) {
    counts_1 <- t(scale(t(as.matrix(counts_1))))
    counts_2 <- t(scale(t(as.matrix(counts_2))))
  }
  #-------------------------------------------------------------------------#
  # Compute score ratio scores
  #-------------------------------------------------------------------------#
  score <- compute_score(counts_1, counts_2, method = method)
  if (center && grepl("mean", method[1L])) {
    score <- score - 1
  } else if (center && grepl("sub", method[1L])) {
    score <- score - min(score)
  } else if (!center && !raw){
    score <- (score - min(score)) / (
      max(score) - min(score))
  }
  
  #-------------------------------------------------------------------------#
  # If collapse is TRUE we take the mean or sum across all cells and
  # return a single value, otherwise by individual cells
  #-------------------------------------------------------------------------#
  if (collapse && grepl("mean", method[1L])) {
    return(mean(score))
  } else if (collapse && grepl("sub", method[1L])) {
    return(sum(score))
  } else {
    return(score)
  }
}




#' Compute score from two count matrices
#'
#' Computes ratio or subtraction scores between two count matrices.
#'
#' @param counts_1 Matrix or vector representing counts for first gene set.
#' @param counts_2 Matrix or vector representing counts for second gene set.
#' @param method Character string specifying "ratio" or "sub"
#'   (subtraction) method.
#' @return Numeric vector of computed scores.
compute_score <- function(counts_1, counts_2, method) {
  counts_1 <- collapse_score(counts_1, method)
  if (all(rownames(counts_1) == rownames(counts_2))) {
    counts_2 <- rep(0, length(counts_1))
  } else {
    counts_2 <- collapse_score(counts_2, method)
  }
  if (grepl("ratio", method)) {
    score <- (counts_1 + 1) / (counts_2 + 1)
  } else if (grepl("sub", method)) {
    score <- counts_1 - counts_2
  }
  return(score)
}

#' Collapse counts matrix to vector
#'
#' Collapses a counts matrix to a single vector using mean or sum aggregation.
#'
#' @param counts Matrix or numeric vector of counts.
#' @param method Character string containing "mean" or "sum" to specify
#'   aggregation method.
#' @return Numeric vector of aggregated values.
collapse_score <- function(counts, method) {
  if (nrow(counts) < 2) {
    return(counts)
  }
  if (grepl("mean", method)) {
    return(colMeans(counts, na.rm = TRUE))
  }
  if (grepl("sum", method)) {
    return(colSums(counts, na.rm = TRUE))
  }
}
