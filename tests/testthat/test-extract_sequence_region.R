
test_that("extract_sequence_region()", {
    use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "heavy_sequences.fasta")
    AIRR_df <- igblastn(query)

    ## Extract regions as nucleotide sequences.

    all_fwrcdr_regions <- lapply(igblastr:::FWRCDR_NAMES,
        function(region_type) extract_sequence_region(AIRR_df, region_type)
    )
    for (regions in all_fwrcdr_regions) {
        expect_true(is.character(regions))
        expect_equal(length(regions), nrow(AIRR_df))
        expect_false(anyNA(regions))
    }
    expected <- as.list(AIRR_df[igblastr:::FWRCDR_NAMES])
    expect_identical(all_fwrcdr_regions, unname(expected))

    concat1 <- do.call(paste0, expected)
    concat2 <- with(AIRR_df, substr(sequence, fwr1_start, fwr4_end))
    expect_identical(concat1, concat2)

    for (region_type in igblastr:::VDJ_REGION_TYPES) {
        regions <- extract_sequence_region(AIRR_df, region_type, as.aa=TRUE)
        expect_true(is.character(regions))
        expect_equal(length(regions), nrow(AIRR_df))
    }

    ## Extract regions as amino acid sequences.

    all_fwrcdr_aa_regions <- lapply(igblastr:::FWRCDR_NAMES,
        function(region_type) extract_sequence_region(AIRR_df, region_type,
                                                      as.aa=TRUE)
    )
    for (regions in all_fwrcdr_aa_regions) {
        expect_true(is(regions, "AAStringSet"))
        expect_equal(length(regions), nrow(AIRR_df))
        expect_false(anyNA(regions))
    }
    expected <- as.list(AIRR_df[paste0(igblastr:::FWRCDR_NAMES, "_aa")])
    expect_identical(lapply(all_fwrcdr_aa_regions, as.character),
                     unname(expected))

    concat1 <- do.call(paste0, expected)
    concat2 <- substr(AIRR_df$sequence_aa, 1L, nchar(concat1))
    expect_identical(concat1, concat2)

    for (region_type in igblastr:::VDJ_REGION_TYPES) {
        regions <- extract_sequence_region(AIRR_df, region_type, as.aa=TRUE)
        expect_true(is(regions, "AAStringSet"))
        expect_equal(length(regions), nrow(AIRR_df))
        if (region_type == "V") {
            expected <- substr(AIRR_df$sequence_aa, 1L, nchar(regions))
            expect_true(all(regions == expected))
        } else {
            counts <- vapply(seq_along(regions),
                function(i) {
                    region <- regions[[i]]
                    if (length(region) == 0L)
                        return(NA_integer_)
                    countPattern(region, AAString(AIRR_df$sequence_aa[[i]]))
                }, integer(1))
            if (region_type == "D") {
                expect_true(all(counts >= 1L, na.rm=TRUE))
            } else {
                expect_true(all(counts == 1L))
            }
        }
    }
})

