data(Kuresi)



test_that("Compute Simple Ratio score", {
    out <- compute_ratio_score(counts,
        genes_1 = c("Gene1", "Gene2"))

})