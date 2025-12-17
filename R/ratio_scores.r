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
    add_name = NULL,
    verbose = TRUE) {
  #-------------------------------------------------------------------------#
  # Get groupings
  #-------------------------------------------------------------------------#
  group_id <- check_grouping(groups, group_name)
  group_id <- split(group_id, group_id[, "group_id"])
 
  for (g in seq_along(group_id)) {
    scores <- compute_ratio_score(counts,
                                  cells = group_id[[g]]$cell_id,
                                  genes_1 = genes_1,
                                  genes_2 = genes_2,
                                  method = method,
                                  weights_1 = weights_1,
                                  weights_2 = weights_2,
                                  scale = scale,
                                  center = center,
                                  collapse = TRUE,
                                  verbose = TRUE)
    groups[[g]]$score <- scores
  }
  groups <- do.call("rbind", groups)
  return(groups)
}




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
    collapse = FALSE,
    verbose = TRUE) {
  simple_bar(verbose)
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
  # Compute score ration scores
  #-------------------------------------------------------------------------#
  score <- compute_score(counts_1, counts_2, method = method)
  if (center && grepl("mean", method[1L])) {
    score <- score - 1
  } else if (center && grepl("sub", method[1L])) {
    score <- score - min(score)
  } else {
    score <- (score - min(score)) / (max(score) - min(score))
  }
  #-------------------------------------------------------------------------#
  # If collapse is TRUE we take the mean or sum across all cells and
  # return a single value, otherwise by indiv cells
  #-------------------------------------------------------------------------#
  if (collapse && grepl("mean", method[1L])) {
    return(mean(score))
  } else if (collapse && grepl("sub", method[1L])) {
    return(sum(score))
  } else {
    return(score)
  }
}




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
