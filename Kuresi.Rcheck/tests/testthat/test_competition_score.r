test_that("outcomes_as_score_matrix works", {
  # Create mock outcomes structure
  scores <- list(
    "Group1" = list(
      "Group1" = c(0.5, 0.5, 0.5),
      "Group2" = c(1, 1, 0.5),
      "Group3" = c(0, 0.5, 1)),
    "Group2" = list(
      "Group1" = c(0, 0.5, 1),
      "Group2" = c(0.5, 0.5, 0.5),
      "Group3" = c(1, 1, 0.5)),
    "Group3" = list(
      "Group1" = c(1, 0.5, 0),
      "Group2" = c(0.5, 0, 0),
      "Group3" = c(0.5, 0.5, 0.5)))
  
  result <- outcomes_as_score_matrix(scores)
  
  expect_type(result, "double")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("outcomes_as_cost works", {
  scores <- list(
    "Group1" = list(
      "Group1" = c(0.5, 0.5, 0.5),
      "Group2" = c(1, 1, 0.5)),
    "Group2" = list(
      "Group1" = c(0, 0.5, 1),
      "Group2" = c(0.5, 0.5, 0.5)))
  
  result <- outcomes_as_cost(scores)
  
  expect_type(result, "double")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("outcomes_as_elo works", {
  scores <- list(
    "Group1" = list(
      "Group1" = c(0.5, 0.5, 0.5),
      "Group2" = c(1, 1, 0.5)),
    "Group2" = list(
      "Group1" = c(0, 0.5, 1),
      "Group2" = c(0.5, 0.5, 0.5)))
  
  result <- outcomes_as_elo(scores, n_tournaments = 10)
  
  expect_type(result, "list")
  expect_named(result, c("elo_seed", "elo_query"))
  expect_type(result$elo_seed, "double")
  expect_type(result$elo_query, "double")
  expect_length(result$elo_seed, 2)
  expect_length(result$elo_query, 2)
  expect_true(all(result$elo_seed > 0))
  expect_true(all(result$elo_query > 0))
})

test_that("compute_competition_outcomes works", {
  # Create test data
  counts <- matrix(
    rpois(200, 5),
    nrow = 20,
    ncol = 20,
    dimnames = list(
      paste0("Gene", 1:20),
      paste0("Cell", 1:20)))
  
  groups <- data.frame(
    cell_id = paste0("Cell", 1:20),
    group_id = rep(c("A", "B"), each = 10),
    row.names = paste0("Cell", 1:20))
  
  gene_set1 <- paste0("Gene", 1:5)
  gene_set2 <- paste0("Gene", 6:10)
  
  result <- compute_competition_outcomes(
    counts = counts,
    groups = groups,
    group_name = "group_id",
    gene_set1 = gene_set1,
    gene_set2 = gene_set2,
    method = "wilcox")
  
  expect_type(result, "list")
  expect_length(result, 2)
  expect_named(result, c("A", "B"))
  expect_true(all(sapply(result, is.list)))
})

