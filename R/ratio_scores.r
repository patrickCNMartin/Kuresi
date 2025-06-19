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
    add_name = NULL,
    verbose = TRUE){
    # simple_bar(verbose)
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
        weights_2)
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
    } else if (center && grepl("sub", method[1L])){
        score <- score - min(score)
    } else {
        score <- (score - min(score)) / (max(score) - min(score))
    }
    return(score)
}




compute_score <- function(counts_1, counts_2, method) {
    counts_1 <- collapse_score(counts_1, method)
    if (all(rownames(counts_1) == rownames(counts_2))){
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
