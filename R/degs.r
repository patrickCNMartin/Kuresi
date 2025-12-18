#' Perform Wilcoxon rank-sum test for differential expression
#'
#' Computes differential expression between two count matrices using
#' Wilcoxon rank-sum tests for each gene.
#'
#' @param s_counts Matrix or sparse matrix of counts for seed/reference
#'   group.
#' @param q_counts Matrix or sparse matrix of counts for query group.
#' @param genes Character vector of gene names (optional, uses rownames
#'   if not provided).
#' @return Data frame with columns: genes, p_value, p_value_adj,
#'   fold_change, effect_size.
#' @importFrom Matrix rowMeans
#' @importFrom stats wilcox.test p.adjust
k_wilcox <- function(s_counts, q_counts, genes) {
  fc <- rowMeans(s_counts) - rowMeans(q_counts)
  genes <- rownames(s_counts)
  p_val <- sapply(
    seq(1, nrow(s_counts)), function(idx, s_counts, q_counts) {
      return(
        suppressWarnings(
          wilcox.test(
            as.vector(s_counts[idx, ]),
            as.vector(q_counts[idx, ]))$p.value))
    },
    s_counts = s_counts,
    q_counts = q_counts
  )
  effect_size <- sapply(p_val, effect_size,
    s_counts = ncol(s_counts),
    q_counts = ncol(q_counts)
  )
  degs <- data.frame(
    "genes" = genes,
    "p_value" = p_val,
    "p_value_adj" = p.adjust(p_val, "bonferroni"),
    "fold_change" = fc,
    "effect_size" = effect_size
  )
  return(degs)
}

#' Perform t-test for differential expression
#'
#' Computes differential expression between two count matrices using
#' Student's t-tests for each gene.
#'
#' @param s_counts Matrix or sparse matrix of counts for seed/reference
#'   group.
#' @param q_counts Matrix or sparse matrix of counts for query group.
#' @param genes Character vector of gene names (optional, uses rownames
#'   if not provided).
#' @return Data frame with columns: genes, p_value, p_value_adj,
#'   fold_change, effect_size.
#' @importFrom Matrix rowMeans
#' @importFrom stats t.test p.adjust
k_ttest <- function(s_counts, q_counts, genes) {
  fc <- rowMeans(s_counts) - rowMeans(q_counts)
  genes <- rownames(s_counts)
  p_val <- sapply(
    seq(1, nrow(s_counts)), function(idx, s_counts, q_counts) {
      return(
        suppressWarnings(
          t.test(s_counts[idx, ], q_counts[idx, ])$p.value))
    },
    s_counts = s_counts,
    q_counts = q_counts
  )
  effect_size <- sapply(p_val, effect_size,
    s_counts = ncol(s_counts),
    q_counts = ncol(q_counts)
  )
  degs <- data.frame(
    "genes" = genes,
    "p_value" = p_val,
    "p_value_adj" = p.adjust(p_val, "bonferroni"),
    "fold_change" = fc,
    "effect_size" = effect_size
  )
  return(degs)
}


#' Calculate effect size from p-value
#'
#' Estimates effect size (Cohen's h) from p-value and sample sizes using
#' power analysis methods.
#'
#' @param pval Numeric p-value.
#' @param s_counts Integer number of samples in seed/reference group.
#' @param q_counts Integer number of samples in query group.
#' @return Numeric effect size value or NA if pval is NA.
#' @importFrom pwr pwr.2p2n.test
effect_size <- function(pval, s_counts, q_counts) {
  if (is.na(pval)) {
    return(NA)
  }
  pval <- max(c(pval, 1e-100))
  pval <- min(c(pval, 0.8))
  effect_size <- pwr::pwr.2p2n.test(
    h = NULL,
    n1 = s_counts,
    n2 = q_counts,
    sig.level = pval,
    power = 0.8
  )$h
  return(effect_size)
}
