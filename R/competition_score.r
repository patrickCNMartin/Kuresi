#-----------------------------------------------------------------------------#
################################ ELO COST #####################################
#-----------------------------------------------------------------------------#
#' @export
compute_competition_outcomes <- function(
    counts,
    groups,
    gene_set1,
    group_name = NULL,
    gene_set2 = NULL,
    method = "wilcox",
    log_fc = 0,
    pval = 0.05) {
  #-------------------------------------------------------------------------#
  # Get groupings
  #-------------------------------------------------------------------------#
  group_id <- check_grouping(groups, group_name)
  group_id <- split(group_id, group_id[, "group_id"]) 
  #-------------------------------------------------------------------------#
  # Verify genes sets
  #-------------------------------------------------------------------------#
  gene_set1 <- check_gene_set(counts, gene_set1)
  gene_set2 <- check_gene_set(counts, gene_set2)
  gene_set <- c(gene_set1, gene_set2)
  if (is.null(gene_set)) {
    stop("Please provide at least one gene set (gene_set1 or gene_set2)")
  }
  #-------------------------------------------------------------------------#
  # Looping over grouping subsets
  #-------------------------------------------------------------------------#
  degs <- vector("list", length(group_id))
  names(degs) <- names(group_id)
  degs <- lapply(degs, function(x, group_id) {
    tmp <- vector("list", length(group_id))
    names(tmp) <- names(group_id)
    return(tmp)
  }, group_id = group_id)
  
  for (s in seq_along(group_id)) {
    for (q in seq_along(group_id)) {
      s_counts <- as.matrix(counts[gene_set, group_id[[s]]$cell_id])
      q_counts <- as.matrix(counts[gene_set, group_id[[q]]$cell_id])
      tmp <- switch(EXPR = method,
        "wilcox" = k_wilcox(s_counts, q_counts)
      )
      degs[[s]][[q]] <- get_outcomes(tmp,
                                     gene_set1,
                                     gene_set2,
                                     fold_change = log_fc,
                                     pval = pval)
    }
  }
  return(degs)
}




get_outcomes <- function(
    deg,
    gene_set1,
    gene_set2,
    fold_change = 0.1,
    pval = 0.05) {
  total_genes <- length(c(gene_set1, gene_set2))
  gene_set1 <- deg$genes %in% gene_set1
  gene_set2 <- deg$genes %in% gene_set2
  outcome <- rep(0, nrow(deg))
  outcome[gene_set1] <- as.numeric(deg$fold_change[gene_set1] > fold_change)
  outcome[gene_set2] <- as.numeric(deg$fold_change[gene_set2] < fold_change)
  outcome[deg$p_value_adj > pval] <- 0.5
  outcome <- c(outcome, rep(0.5, times = total_genes - length(outcome)))
  return(outcome)
}

#' @export
outcomes_as_score_matrix <- function(scores) {
    cols <- length(scores)
    rows <- length(scores[[1]])
    score_matrix <- matrix(0, nrow = rows, ncol = cols)
    colnames(score_matrix) <- names(scores)
    rownames(score_matrix) <- names(scores[[1]])
    for (i in seq_along(scores)) {
      for (j in seq_along(scores[[i]])) {
        score_matrix[j, i] <- (sum(scores[[i]][[j]]) / length(scores[[i]][[j]]))
      }
    }
  return(score_matrix)
}

#' @export
outcomes_as_cost <- function(scores) {
    cost_matrix <- matrix(0, nrow = length(scores[[1]]), ncol = length(scores))
    colnames(cost_matrix) <- names(scores)
    rownames(cost_matrix) <- names(scores[[1]])
    for (i in seq_along(scores)) {
      for (j in seq_along(scores[[i]])) {
        local_score <- (sum(scores[[i]][[j]]) / length(scores[[i]][[j]])) - 0.5
        local_score <- 1 - abs(local_score)
        cost_matrix[j, i] <- local_score
      }
    }
  return(cost_matrix)
}

#' @export
outcomes_as_elo <- function(
    scores,
    initial_elo = 1000,
    k = 32,
    n_tournaments = 100) {
  scores <- outcomes_as_score_matrix(scores)
  elo_seed <- rep(initial_elo, ncol(scores))
  names(elo_seed) <- colnames(scores)
  elo_query <- rep(initial_elo, nrow(scores))
  names(elo_query) <- rownames(scores)
  for (n in seq_len(n_tournaments)) {
    for (i in seq_len(ncol(scores))) {
      for (j in seq_len(nrow(scores))) {
        local_seed_rating <- elo_seed[i]
        local_query_rating <- elo_query[j]
        local_scores <- ifelse(scores[j, i] == 0.5, 0.5, round(scores[j, i]))
        compute_elo(local_seed_rating, local_query_rating, local_scores)
        elo_seed[i] <- local_seed_rating
        elo_query[j] <- local_query_rating
      }
    }
  }

  return(list("elo_seed" = elo_seed, "elo_query" = elo_query))
}



compute_elo <- function(
    local_seed_rating,
    local_query_rating,
    outcome_seed,
    k = 32) {
  expected_seed <- 1 / (1 + 10^((local_query_rating - local_seed_rating) / 400))
  new_local_seed_rating <- local_seed_rating +
    k * (outcome_seed - expected_seed)
  outcome_query <- 1 - outcome_seed
  expected_query <- 1 /
    (1 + 10^((local_seed_rating - local_query_rating) / 400))
  new_local_query_rating <- local_query_rating +
    k * (outcome_query - expected_query)
  assign("local_seed_rating",
         round(new_local_seed_rating),
         env = parent.frame())
  assign("local_query_rating",
         round(new_local_query_rating),
         env = parent.frame())
}
