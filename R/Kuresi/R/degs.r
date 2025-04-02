k_wilcox <- function(s_counts, q_counts) {
    log2_fc <- log2(rowMeans(s_counts) - rowMeans(q_counts))
    genes <- rownames(s_counts)
    p_val <- sapply(
        seq(1, nrow(s_counts)), function(idx, s_counts, q_counts){
            return(
                suppressWarnings(
                wilcox.test(s_counts[idx, ], q_counts[idx, ])$p.value))
    })
    effect_size <- sapply(p_val, effect_size,
        s_counts = ncol(s_counts),
        q_counts = ncol(q_counts))
    degs <- data.frame(
        "genes" = genes,
        "p_value" = p_val,
        "p_value_adj" = p.adjust(p_val, "bonferroni"),
        "log2_fc" = log2_fc,
        "effect_size" = effect_size
    )
    return(degs)
}

k_ttest <- function(s_counts, q_counts) {
    log2_fc <- log2(rowMeans(s_counts) - rowMeans(q_counts))
    genes <- rownames(s_counts)
    p_val <- sapply(
        seq(1, nrow(s_counts)), function(idx, s_counts, q_counts){
            return(
                suppressWarnings(
                t.test(s_counts[idx, ], q_counts[idx, ])$p.value))
    })
    effect_size <- sapply(p_val, effect_size,
        s_counts = ncol(s_counts),
        q_counts = ncol(q_counts))
    degs <- data.frame(
        "genes" = genes,
        "p_value" = p_val,
        "p_value_adj" = p.adjust(p_val, "bonferroni"),
        "log2_fc" = log2_fc,
        "effect_size" = effect_size
    )
    return(degs)
}

# k_lm <- function(s_counts, q_counts) {
#     log2_fc <- log2(rowMeans(s_counts) - rowMeans(q_counts))
#     genes <- rownames(s_counts)

# }

effect_size <- function(pval, s_counts, q_counts) {
  if (is.na(pval)) {
      return(NA)
  }
  pval <- max(c(pval, 1e-100))
  pval <- min(c(pval, 0.8))
  effect_size <- pwr::pwr.2p2n.test(h = NULL,
    n1 = seed,
    n2 = query,
    sig.level = pval,
    power = 0.8)$h
  return(effect_size)
}