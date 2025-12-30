test_that("win_lose_genes returns expected structure", {
  result <- win_lose_genes()
  
  expect_type(result, "list")
  expect_named(result, c("win", "lose"))
  expect_type(result$win, "character")
  expect_type(result$lose, "character")
  expect_true(length(result$win) > 0)
  expect_true(length(result$lose) > 0)
})

test_that("cancer_maker_list returns genes for valid types", {
  pancreas <- cancer_maker_list("pancreas")
  breast <- cancer_maker_list("breast")
  ovary <- cancer_maker_list("ovary")
  
  expect_type(pancreas, "character")
  expect_type(breast, "character")
  expect_type(ovary, "character")
  expect_true(length(pancreas) > 0)
  expect_true(length(breast) > 0)
  expect_true(length(ovary) > 0)
})

