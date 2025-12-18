test_that("compute_ratio_score works with simple input", {
  # Create simple test data
  counts <- matrix(
    rpois(100, 5),
    nrow = 10,
    ncol = 10,
    dimnames = list(
      paste0("Gene", 1:10),
      paste0("Cell", 1:10)))
  
  result <- compute_ratio_score(
    counts = counts,
    genes_1 = c("Gene1", "Gene2", "Gene3"),
    genes_2 = c("Gene4", "Gene5"),
    method = "mean_ratio",
    collapse = FALSE)
  
  expect_type(result, "double")
  expect_length(result, 10) # one score per cell
  expect_true(all(!is.na(result)))
  expect_true(all(result >= 0))
})

test_that("compute_ratio_score works with collapse = TRUE", {
  counts <- matrix(
    rpois(100, 5),
    nrow = 10,
    ncol = 10,
    dimnames = list(
      paste0("Gene", 1:10),
      paste0("Cell", 1:10)))
  
  result <- compute_ratio_score(
    counts = counts,
    genes_1 = c("Gene1", "Gene2"),
    method = "mean_ratio",
    collapse = TRUE)
  
  expect_type(result, "double")
  expect_length(result, 1) # single aggregated value
  expect_true(!is.na(result))
})

test_that("compute_ratio_score works with different methods", {
  counts <- matrix(
    rpois(100, 5),
    nrow = 10,
    ncol = 10,
    dimnames = list(
      paste0("Gene", 1:10),
      paste0("Cell", 1:10)))
  
  methods <- c("mean_ratio", "mean_sub", "sum_ratio", "sum_sub")
  
  for (method in methods) {
    result <- compute_ratio_score(
      counts = counts,
      genes_1 = c("Gene1", "Gene2"),
      genes_2 = c("Gene3", "Gene4"),
      method = method,
      collapse = FALSE)
    
    expect_type(result, "double")
    expect_length(result, 10)
    expect_true(all(!is.na(result)))
  }
})

test_that("compute_ratio_bygroup works", {
  counts <- matrix(
    rpois(100, 5),
    nrow = 10,
    ncol = 10,
    dimnames = list(
      paste0("Gene", 1:10),
      paste0("Cell", 1:10)))
  
  groups <- data.frame(
    cell_id = paste0("Cell", 1:10),
    group = rep(c("A", "B"), each = 5),
    row.names = paste0("Cell", 1:10))
  
  result <- compute_ratio_bygroup(
    counts = counts,
    groups = groups,
    group_name = "group",
    genes_1 = c("Gene1", "Gene2"),
    genes_2 = c("Gene3", "Gene4"),
    method = "mean_ratio")
  
  expect_s3_class(result, "data.frame")
  expect_true("mean_ratio" %in% colnames(result))
  expect_equal(nrow(result), 10)
})
