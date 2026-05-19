
test_that("make_unique_allele_names()", {
    make_unique_allele_names <- igblastr:::make_unique_allele_names

    allele_names <- c(A="ID1", B="ID2")
    expect_identical(make_unique_allele_names(allele_names), allele_names)
    expect_identical(make_unique_allele_names(allele_names, suffixes.only=TRUE),
                     setNames(character(2), names(allele_names)))

    allele_names <- c(A="ID1", B="ID2", C="ID1", D="ID1", E="ID3", F="ID3")
    expected <- c(A="a", B="", C="b", D="c", E="a", F="b")
    expect_identical(make_unique_allele_names(allele_names, suffixes.only=TRUE),
                     expected)
    expect_identical(make_unique_allele_names(allele_names),
                     add_suffix(allele_names, expected))

    allele_names <- c(allele_names, G="ID1c")
    regexp <- "disambiguation scheme"
    expect_error(make_unique_allele_names(allele_names), regexp=regexp)
    expect_error(make_unique_allele_names(allele_names, suffixes.only=TRUE),
                 regexp=regexp)
})

