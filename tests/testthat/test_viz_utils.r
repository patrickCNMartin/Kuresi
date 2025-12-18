test_that("view_scores works", {
  # Create simple score matrix
  scores <- matrix(
    runif(9, 0, 1),
    nrow = 3,
    ncol = 3,
    dimnames = list(
      paste0("Query", 1:3),
      paste0("Seed", 1:3)))
  
  result <- view_scores(scores)
  
  expect_s3_class(result, "ggplot")
  # ggplot objects should not be NULL
  expect_true(!is.null(result))
})

test_that("view_elo works", {
  elo_scores <- list(
    elo_seed = c("Group1" = 1200, "Group2" = 1000, "Group3" = 1100),
    elo_query = c("Group1" = 1150, "Group2" = 1050))
  
  result <- view_elo(elo_scores)
  
  expect_s3_class(result, "ggplot")
  expect_true(!is.null(result))
})

test_that("generate_palette works", {
  palette <- c("red", "blue", "green")
  seed <- c("A", "B", "C")
  query <- c("B", "C", "D")
  
  result <- generate_palette(palette, seed, query)
  
  expect_type(result, "list")
  expect_named(result, c("seed", "query"))
  expect_type(result$seed, "character")
  expect_type(result$query, "character")
  expect_length(result$seed, 3)
  expect_length(result$query, 3)
})

test_that("score_plot works", {
  score_data <- data.frame(
    x = 1:10,
    y = 1:10,
    score = runif(10, 0, 1))
  
  result <- score_plot(score_data)
  
  expect_s3_class(result, "ggplot")
  expect_true(!is.null(result))
})

test_that("score_plot works with bins", {
  score_data <- data.frame(
    x = 1:10,
    y = 1:10,
    score = runif(10, 0, 1))
  
  result <- score_plot(score_data, bins = 3)
  
  expect_s3_class(result, "ggplot")
  expect_true(!is.null(result))
})

