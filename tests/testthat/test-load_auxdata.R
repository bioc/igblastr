
.AUXDATA_COLNAMES <- names(igblastr:::AUXDATA_COL2CLASS)

test_that("load_auxdata()", {
    organisms <- list_igblast_organisms()
    for (organism in organisms) {
        auxdata <- load_auxdata(organism, which="original")
        expect_true(is.data.frame(auxdata))
        expect_identical(colnames(auxdata), .AUXDATA_COLNAMES)
        ## 1 row is repeated in human_gl.aux (the row for TRAJ13*02)
        if (organism == "human") {
            ok <- igblastr:::rows_with_same_key_are_identical(auxdata,
                                                              "allele_name")
            expect_true(ok)
        } else {
            expect_identical(anyDuplicated(auxdata[ , "allele_name"]), 0L)
        }
    }
})

