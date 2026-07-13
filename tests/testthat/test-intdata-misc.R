
test_that("V_genes_with_varying_fwrcdr_boundaries()", {
    human_intdata0 <- load_intdata("human")
    ## The results below might change when NCBI updates the internal data
    ## (e.g. when they release the next IgBLAST release). In that case the
    ## examples in man/intdata-misc.Rd might also need some adjustments.
    var_genes <- V_genes_with_varying_fwrcdr_boundaries(human_intdata0)
    expected <- c("IGHV4-31", "IGHV4-4", "IGHV5-a", "TRBV30")
    expect_identical(var_genes, expected)
    var_genes <- V_genes_with_varying_fwrcdr_boundaries(human_intdata0,
                                                        V_segment="cdr1")
    expect_identical(var_genes, expected[c(2L, 4L)])
})

