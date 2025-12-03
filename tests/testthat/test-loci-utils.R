test_that("normalize_user_supplied_loci()", {
    normalize_user_supplied_loci <- igblastr:::normalize_user_supplied_loci
    expect_identical(normalize_user_supplied_loci(),
                     c("IGH", "IGK", "IGL"))
    expect_identical(normalize_user_supplied_loci(tcr.db=TRUE),
                     c("TRA", "TRB", "TRG", "TRD"))
    expect_identical(normalize_user_supplied_loci("IGL+IGK"),
                     c("IGK", "IGL"))
    expect_identical(normalize_user_supplied_loci("TRG+TRB+TRD+TRA"),
                     c("TRA", "TRB", "TRG", "TRD"))
    expect_identical(normalize_user_supplied_loci("TRG+TRD+TRA",
                                                  stop.if.missing.regions=TRUE),
                     c("TRA", "TRG", "TRD"))

    ## 'loci' must be a non-empty character vector.
    expect_error(normalize_user_supplied_loci(character(0)))
    ## 'loci' cannot contain NAs or duplicates
    expect_error(normalize_user_supplied_loci(c("IGH", "IGH")))
    expect_error(normalize_user_supplied_loci(c("IGH", NA)))
    ## 'tcr.db' should not be used when 'loci' is supplied.
    expect_error(normalize_user_supplied_loci("IGH", tcr.db=TRUE))
    ## Invalid locus.
    expect_error(normalize_user_supplied_loci("TRC"))
    ## Invalid locus selection.
    expect_error(normalize_user_supplied_loci("IGH+TRA"))
    ## No D regions.
    expect_error(
        normalize_user_supplied_loci("IGL", stop.if.missing.regions=TRUE)
    )
    expect_error(
        normalize_user_supplied_loci("IGK+IGL", stop.if.missing.regions=TRUE)
    )
    expect_error(
        normalize_user_supplied_loci("TRA", stop.if.missing.regions=TRUE)
    )
    expect_error(
        normalize_user_supplied_loci("TRA+TRG", stop.if.missing.regions=TRUE)
    )
})

