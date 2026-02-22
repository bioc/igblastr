test_that(".drop_repeated_alleles()", {
    drop_repeated_alleles <- igblastr:::.drop_repeated_alleles

    x <- c(X="ACTT", Y="ACTT", X="CC", Y="TT", Y="GGGGGGGG",
           X="ACTT", X="ACTT", Z="", Y="TT", Y="GgGgGgGg")

    expect_identical(drop_repeated_alleles(x), x[-c(6:7, 9)])

    dna <- DNAStringSet(x)
    object <- drop_repeated_alleles(dna)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), x[-c(6:7, 9:10)])

    y <- rev(x)
    expect_identical(drop_repeated_alleles(y), y[-c(5L, 7L, 10L)])

    dna <- DNAStringSet(y)
    object <- drop_repeated_alleles(dna)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), toupper(y[-c(5:7, 10L)]))
})

.check_object_mcols <- function(object)
{
    object_mcols <- mcols(object)
    expect_true(is(object_mcols, "DataFrame"))
    expected_colnames <- c(
        "locus",
        setdiff(names(igblastr:::IGBLAST_INTDATA_COL2CLASS), "chain_type"),
        "starting_gap", "all_gaps_in_frame", "all_gaps_contained",
        "chain_type"
    )
    expect_identical(colnames(object_mcols), expected_colnames)
    expect_identical(object_mcols[ , "allele_name"], names(object))
}

test_that("clean_allele_set()", {
    clean_allele_set <- igblastr:::clean_allele_set

    ## ======== WITHOUT GAPS ========

    dna0 <- DNAStringSet(c(A1="GG", A2="AAAAACG", A3="T", A4="GG", A5="CTAATA",
                           A6="GG", A7="TTT", A8="GG", A9="TTTTTTTTTT"))
    object <- clean_allele_set(dna0)
    expect_identical(as.character(object), as.character(dna0))

    dna <- dna0
    names(dna) <- paste0("   ", names(dna), "   \t ")
    object <- clean_allele_set(dna)
    expect_identical(as.character(object), as.character(dna0))

    dna <- dna0
    names(dna) <- paste0("\t stuff |", names(dna), "\t   |more|| stuff  ")
    object <- clean_allele_set(dna)
    expect_identical(as.character(object), as.character(dna0))

    dna <- dna0
    names(dna)[2L] <- NA
    expect_error(clean_allele_set(dna), regexp="NAs")

    dna <- dna0
    names(dna)[2L] <- "   \t\t"
    expect_error(clean_allele_set(dna), regexp="empty")

    dna <- dna0
    names(dna)[2L] <- "Z"
    expect_error(clean_allele_set(dna), regexp="less than 2-character long")

    ## With "repeated" alleles (i.e. alleles with identical DNA
    ## sequences **and** names).
    dna <- dna0
    names(dna)[c(1:2, 9L)] <- "A1"
    names(dna)[c(3L, 5:6, 8L)] <- "A3"
    expected_names <- c("A1a", "A1b", "A3a", "A4", "A3b", "A3c", "A7", "A1c")
    expected_object <- setNames(dna0[-8L], expected_names)
    object <- clean_allele_set(dna)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(expected_object))
    expect_true(is.null(mcols(object)))

    regexp <- "V allele sequences have no gaps"
    expect_warning(
        object <- clean_allele_set(dna, gapped=TRUE),
        regexp=regexp
    )
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(expected_object))
    expect_true(is.null(mcols(object)))

    ## ======== WITH GAPS ========

    dna1 <- DNAStringSet(c(A1="AC.G..T", A2="AAA...AACGT", A3="G...G"))
    errmsg <- "Some allele sequences have gaps"
    expect_error2(clean_allele_set(dna1), errmsg)

    object <- clean_allele_set(dna1, gapped=TRUE)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), c(A1="ACGT", A2="AAAAACGT", A3="GG"))
    expect_true(is.null(mcols(object)))

    errmsg <- paste0("'with.intdata=TRUE' can only be used ",
                     "when 'gapped' is TRUE")
    expect_error2(clean_allele_set(dna1, with.intdata=TRUE), errmsg)

    errmsg <- paste0("'dna' must have a \"locus\" metadata column ",
                     "when 'gapped' is TRUE")
    expect_error2(clean_allele_set(dna1, gapped=TRUE, with.intdata=TRUE),
                  errmsg)

    mcols(dna1)$locus <- "IGH"
    dna <- dna1
    names(dna) <- paste0("\t stuff |", names(dna), "\t   |more|| stuff  ")
    object <- clean_allele_set(dna, gapped=TRUE, with.intdata=TRUE)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), c(A1="ACGT", A2="AAAAACGT", A3="GG"))
    .check_object_mcols(object)

    dna2 <- DNAStringSet(c(A1="ACGT", A2="AAAAACGT", A3="G...G"))
    errmsg <- "Some allele sequences have gaps"
    expect_error2(clean_allele_set(dna2), errmsg)

    dna <- dna2
    names(dna) <- paste0("\t ", names(dna), "\t    ")
    regexp <- "V allele sequences have no gaps"
    expect_warning(
        object <- clean_allele_set(dna, gapped=TRUE),
        regexp=regexp
    )
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), c(A1="ACGT", A2="AAAAACGT", A3="GG"))
    expect_true(is.null(mcols(object)))

    mcols(dna2)$locus <- "IGK"
    dna <- dna2
    names(dna) <- paste0("\t ", names(dna), "\t    ")
    regexp <- "V allele sequences have no gaps"
    expect_warning(
        object <- clean_allele_set(dna, gapped=TRUE, with.intdata=TRUE),
        regexp=regexp
    )
    expect_true(is(object, "DNAStringSet"))
    .check_object_mcols(object)

    dna3 <- DNAStringSet(c(A1="G.G.", A2="AAA........AACG",
                           A3="T..", A4="....GG..", A5="CT..AATA",
                           A6="G....G", A7="T..T..T", A8=".G...G",
                           A9=".........TTT...T.T.T.T..TT.....T"))
    mcols(dna3)$locus <- "TRA"

    object <- clean_allele_set(dna3, gapped=TRUE)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(dna0))
    expect_identical(colnames(mcols(object)), colnames(mcols(dna3)))

    object <- clean_allele_set(dna3, gapped=TRUE, with.intdata=TRUE)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(dna0))
    .check_object_mcols(object)

    ## With "repeated" alleles (i.e. alleles with identical **ungapped** DNA
    ## sequences **and** names).
    dna <- dna3
    names(dna)[c(1:2, 9L)] <- "A1"
    names(dna)[c(3L, 5:6, 8L)] <- "A3"
    expected_names <- c("A1a", "A1b", "A3a", "A4", "A3b", "A3c", "A7", "A1c")
    expected_object <- setNames(dna0[-8L], expected_names)
    object <- clean_allele_set(dna, gapped=TRUE)
    expect_true(is(object, "DNAStringSet"))
    expect_identical(as.character(object), as.character(expected_object))
    expect_identical(colnames(mcols(object)), colnames(mcols(dna3)))

    ## Fails because the gaps in A6 and A8 imply different coding_frame_start.
    errmsg <- paste0("alleles with identical ungapped sequences and names ",
                     "must have gaps that result in identical annotations")
    expect_error2(clean_allele_set(dna, gapped=TRUE, with.intdata=TRUE), errmsg)

    ## Let's put gaps in A6 and A8 that result in the two alleles having
    ## identical annotations.
    dna <- dna3
    dna[c("A6", "A8")] <- c(".......GG", ".G....G..")
    mcols(dna)$locus <- "TRB"
    names(dna)[c(1:2, 9L)] <- "A1"
    names(dna)[c(3L, 5:6, 8L)] <- "A3"
    object <- clean_allele_set(dna, gapped=TRUE, with.intdata=TRUE)
    expect_identical(as.character(object), as.character(expected_object))
    .check_object_mcols(object)
})

