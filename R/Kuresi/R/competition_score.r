#-----------------------------------------------------------------------------#
################################ ELO COST #####################################
#-----------------------------------------------------------------------------#
compute_competition_outcomes <- function(
    counts,
    grouping,
    expressed,
    grouping_name = NULL,
    repressed = NULL,
    method = "wilcox",
    log_fc = 0,
    pval = 0.05) {
    #-------------------------------------------------------------------------#
    # Get groupings
    #-------------------------------------------------------------------------#
    grouping_name <- check_grouping_name(grouping, grouping_name)
    grouping <- split(grouping, grouping[, grouping_name])
    #-------------------------------------------------------------------------#
    # Looping over grouping subsets
    #-------------------------------------------------------------------------#
    for (s in seq_along(grouping)){
        for (q in seq_along(grouping)){
            s_counts <- counts[, rownames(grouping[[s]])]
            q_counts <- counts[, rownames(grouping[[q]])]
            degs <- switch(EXPR = method,
                "wilcox" = k_wilcox(s_counts, q_counts),
                "edgeR" = k_edgeR(s_counts, q_counts))
        }
    }
}

#'@importFrom vesalius identify_markers
#' @export
compute_competition_outcomes_old <- function(seed_assay,
    expressed,
    repressed = NULL,
    query_assay = NULL,
    integration_method = "none",
    nfeatures = 2000,
    dimensions = 30,
    use_counts = "raw",
    method = "wilcox",
    log_fc = 0,
    pval = 0.05,
    min_pct = 0,
    min_spatial_index = 10,
    norm_method = "last",
    trial_seed = "last",
    labels_seed = NULL,
    trial_query = "last",
    labels_query = NULL,
    verbose = TRUE) {
    #-------------------------------------------------------------------------#
    # initial sanity check and gets
    #-------------------------------------------------------------------------#
    expressed <- check_comp_genes(expressed, seed_assay, norm_method, method)
    repressed <- check_comp_genes(repressed, seed_assay, norm_method, method)
    signal <- c(expressed, repressed)
    trial <- trial_seed
    seed_trial <- vesalius:::check_territory_trial(seed_assay, trial_seed)
    if (is.null(labels_seed)) {
        labels_seed <- unique(seed_trial$trial)
    }
    seed <- vesalius:::check_group_value(seed_trial, labels_seed)
    seed <- filter_territories(seed, seed_trial, min_spatial_index)
    query <- seed
    sample <- FALSE
    #-------------------------------------------------------------------------#
    # first we check integration 
    #-------------------------------------------------------------------------#
    if (!is.null(query_assay)) {
        trial_seed <- vesalius:::check_territory_trial(seed_assay,
            trial = trial_seed,
            return_label = TRUE)
        trial_query <- vesalius:::check_territory_trial(query_assay,
            trial = trial_query,
            return_label = TRUE)
        seed_counts <- get_counts(seed_assay, type = use_counts)
        query_counts <- get_counts(query_assay, type = use_counts)
        integrated <- count_integration(seed_counts,
            query_counts,
            integration_method = integration_method,
            nfeatures = nfeatures,
            dimensions = dimensions,
            features = signal)
        trial <- paste0(trial_seed, "_", trial_query)
        query_trial <- vesalius:::check_territory_trial(query_assay, trial_query)
        if (is.null(labels_query)){
            labels_query <- unique(query_trial$trial)
        }
        query <- vesalius:::check_group_value(query_trial, labels_query)
        query <- filter_territories(query, query_trial, min_spatial_index)
        sample <- TRUE

        buffer <- quick_assay_build(seed_assay,
            query_assay,
            integrated,
            trial_seed,
            trial_query,
            seed,
            query,
            min_spatial_index)
        seed_assay <- buffer$assay
        seed <- buffer$seed
        query <- buffer$query
    }
    #-------------------------------------------------------------------------#
    # initialize, Loop, and clean
    # checking genes to make sure they are still present after integration
    # if integration was used
    #-------------------------------------------------------------------------#

    signal <- check_comp_genes(signal, seed_assay, norm_method, method)
    scores <- vector("list", length(seed))
    names(scores) <- as.character(seed)
    scores <- lapply(scores, function(i, query){
            tmp <- vector("list", length(query))
            names(tmp) <- as.character(query)
            return(tmp)
        }, query = query)
    for (s in seq_along(scores)) {
        for (q in seq_along(scores[[s]])) {
            message_switch(verbose, "kuresi",seed = seed[s], query = query[q])
            tmp <- vesalius::identify_markers(seed_assay,
                trial = trial,
                norm_method = norm_method,
                seed = seed[s],
                query = query[q],
                sample = sample,
                pval = 1,
                log_fc = log_fc,
                min_pct = min_pct,
                min_spatial_index = 10,
                genes = signal,
                verbose = FALSE)
            tmp <- vesalius::get_markers(tmp)
            
            scores[[s]][[q]] <- get_outcomes(tmp, expressed, repressed, pval = pval)
        }
    }
    return(scores)
}



get_outcomes <- function(deg,
    expressed,
    repressed,
    pval = 0.05) {
    total_genes <- length(c(expressed, repressed))
    expressed <- deg$genes %in% expressed
    repressed <- deg$genes %in% repressed
    outcome <- rep(0, nrow(deg))
    outcome[expressed] <- as.numeric(deg$fold_change[expressed] > 0)
    outcome[repressed] <- as.numeric(deg$fold_change[repressed] < 0)
    outcome[deg$p_value_adj > pval] <- 0.5
    outcome <- c(outcome, rep(0.5, times = total_genes - length(outcome)))
    return(outcome)
}

#' @export
outcomes_as_score <- function(scores,
    seed_assay = NULL,
    query_assay = NULL,
    trial_seed = "last",
    trial_query = "last") {
    if (!is.null(seed_assay) && !is.null(query_assay)){
        seed <- vesalius:::check_territory_trial(seed_assay, trial_seed)
        seed <- seed %>% filter(trial %in% names(scores))
        query <- vesalius:::check_territory_trial(query_assay, trial_query)
        query <- query %>% filter(trial %in% names(scores[[1]]))
        cost_matrix <- matrix(0, nrow = nrow(query), ncol = nrow(seed))
        colnames(cost_matrix) <- seed$barcodes
        rownames(cost_matrix) <- query$barcodes
        for (i in seq_along(scores)){
            for (j in seq_along(scores[[i]])){
                local_score <- (sum(scores[[i]][[j]]) / length(scores[[i]][[j]]))
                query_loc <- query$barcodes[query$trial == names(scores[[i]])[j]]
                seed_loc <- seed$barcodes[seed$trial == names(scores)[i]]
                cost_matrix[query_loc, seed_loc] <- local_score 
            }
        }
        
    } else {
        cols <- length(scores)
        rows <- length(scores[[1]])
        cost_matrix <- matrix(0, nrow = rows, ncol = cols)
        colnames(cost_matrix) <- names(scores)
        rownames(cost_matrix) <- names(scores[[1]])
        for (i in seq_along(scores)){
            for (j in seq_along(scores[[i]])){
                cost_matrix[j,i] <- (sum(scores[[i]][[j]]) / length(scores[[i]][[j]]))
            }
        }
    }
    return(cost_matrix)
}

#' @export
outcomes_as_cost <- function(scores,
    seed_assay = NULL,
    query_assay = NULL,
    trial_seed = "last",
    trial_query = "last") {
    if (!is.null(seed_assay) && !is.null(query_assay)){
        seed <- vesalius:::check_territory_trial(seed_assay, trial_seed)
        seed <- seed %>% filter(trial %in% names(scores))
        query <- vesalius:::check_territory_trial(query_assay, trial_query)
        query <- query %>% filter(trial %in% names(scores[[1]]))
        cost_matrix <- matrix(0, nrow = nrow(query), ncol = nrow(seed))
        colnames(cost_matrix) <- seed$barcodes
        rownames(cost_matrix) <- query$barcodes
        for (i in seq_along(scores)){
            for (j in seq_along(scores[[i]])){
                local_score <- (sum(scores[[i]][[j]]) / length(scores[[i]][[j]])) - 0.5
                local_score <- 1 - abs(local_score)
                query_loc <- query$barcodes[query$trial == names(scores[[i]])[j]]
                seed_loc <- seed$barcodes[seed$trial == names(scores)[i]]
                cost_matrix[query_loc, seed_loc] <- local_score 
            }
        }
        
    } else {
        cost_matrix <- matrix(0, nrow = length(scores[[1]]), ncol = length(scores))
        colnames(cost_matrix) <- names(scores)
        rownames(cost_matrix) <- names(scores[[1]])
        for (i in seq_along(scores)){
            for (j in seq_along(scores[[i]])){
                local_score <- (sum(scores[[i]][[j]]) / length(scores[[i]][[j]])) - 0.5
                local_score <- 1 - abs(local_score)
                cost_matrix[j, i] <- local_score 
            }
        }
        
    }
    return(cost_matrix)
}

#'@export
outcomes_as_elo <- function(scores,
    seed_assay = NULL,
    query_assay = NULL,
    trial_seed = "last",
    trial_query = "last",
    initial_elo = 1000,
    k = 32,
    n_tournaments = 100) {
    scores <- outcomes_as_score(scores)
    elo_seed <- rep(initial_elo, ncol(scores))
    names(elo_seed) <- colnames(scores)
    elo_query <- rep(initial_elo, nrow(scores))
    names(elo_query) <- rownames(scores)
    for (n in seq_len(n_tournaments)){
        for (i in seq_len(ncol(scores))) {
            for (j in seq_len(nrow(scores))) {
                local_seed_rating <- elo_seed[i]
                local_query_rating <- elo_query[j]
                local_scores <- ifelse(scores[j,i] == 0.5, 0.5 ,round(scores[j,i]))
                compute_elo(local_seed_rating, local_query_rating, local_scores)
                elo_seed[i] <- local_seed_rating
                elo_query[j] <- local_query_rating
            }
        }
    }
    
    return(list("elo_seed" = elo_seed, "elo_query" = elo_query))
}


#' @export
assign_score <- function(vesalius_assay,
    score,
    trial= "last",
    signif = 5,
    as_query = FALSE,
    add_name = NULL) {
    if (!is.null(add_name)){
        new_trial <- vesalius:::create_trial_tag(
            names(get_territories(vesalius_assay)), add_name) %>%
            tail(1)
    } else {
        new_trial <- vesalius:::create_trial_tag(
            names(get_territories(vesalius_assay)), 
            "Score") %>%
            tail(1)
    }
    if(as_query){
        scored_territories <- rownames(score)
        score <- 1 - apply(score, MARGIN = 1, mean)
    } else {
        scored_territories <- colnames(score)
        score <- apply(score, MARGIN = 2, mean)
    }
    trial <- vesalius:::check_territory_trial(vesalius_assay, trial)
    locs <- match(trial$trial,scored_territories)
    trial$trial <- "Not Selected"
    trial$trial[!is.na(locs)] <- signif(score[locs[!is.na(locs)]],digits = signif)
    trial <- trial[, c("barcodes", "x","y", "trial")]
    
    colnames(trial) <- gsub("trial", new_trial, colnames(trial))
    vesalius_assay <- vesalius:::update_vesalius_assay(vesalius_assay = vesalius_assay,
        data = trial,
        slot = "territories",
        append = TRUE)
    return(vesalius_assay)
}

#' @importFrom vesalius add_cells
#' @export
assign_elo <- function(vesalius_assay,
    elo,
    trial= "last",
    dim = 1,
    add_name = NULL){
    if (!is.null(add_name)){
        new_trial <- vesalius:::create_trial_tag(
            names(get_territories(vesalius_assay)), add_name) %>%
            tail(1)
    } else {
        new_trial <- vesalius:::create_trial_tag(
            names(get_territories(vesalius_assay)), 
            "ELO") %>%
            tail(1)
    }
    scored_territories <- names(elo[[dim]])
    trial <- vesalius:::check_territory_trial(vesalius_assay, trial)
    score <- elo[[dim]]
    locs <- match(trial$trial,scored_territories)
    trial$trial <- "Not Selected"
    trial$trial[!is.na(locs)] <- score[locs[!is.na(locs)]]
    trial <- trial[, c("barcodes", "x","y", "trial")]
    colnames(trial) <- gsub("trial", new_trial, colnames(trial))
    vesalius_assay <- vesalius:::update_vesalius_assay(vesalius_assay = vesalius_assay,
        data = trial,
        slot = "territories",
        append = TRUE)
    return(vesalius_assay)

   
}


compute_elo <- function(local_seed_rating,
    local_query_rating,
    outcome_seed,
    k = 32) {
    expected_seed <- 1 / (1 + 10^((local_query_rating - local_seed_rating) / 400))
    new_local_seed_rating <- local_seed_rating + k * (outcome_seed - expected_seed)
    outcome_query <- 1 - outcome_seed
    expected_query <- 1 / (1 + 10^((local_seed_rating - local_query_rating) / 400))
    new_local_query_rating <- local_query_rating + k * (outcome_query - expected_query)
    assign("local_seed_rating", round(new_local_seed_rating), env = parent.frame())
    assign("local_query_rating", round(new_local_query_rating), env = parent.frame())
}


