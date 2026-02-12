test_that("drop_repeated_sequences()", {
    drop_repeated_sequences <- igblastr:::drop_repeated_sequences

    x <- c(X="ACTT", Y="ACTT", X="CC", Y="TT", Y="GGGGGGGG",
           X="ACTT", X="ACTT", Z="", Y="TT", Y="GgGgGgGg")

    expect_identical(drop_repeated_sequences(x), x[-c(6:7, 9)])

    dna <- DNAStringSet(x)
    object <- drop_repeated_sequences(dna)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), x[-c(6:7, 9:10)])

    y <- rev(x)
    expect_identical(drop_repeated_sequences(y), y[-c(5L, 7L, 10L)])

    dna <- DNAStringSet(y)
    object <- drop_repeated_sequences(dna)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), toupper(y[-c(5:7, 10L)]))
})

