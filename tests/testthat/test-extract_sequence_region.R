
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
    reconstructed <- do.call(paste0, all_fwrcdr_regions)
    expect_identical(reconstructed,
                     with(AIRR_df, substr(sequence, fwr1_start, fwr4_end)))

    all_VDJ_regions <- lapply(igblastr:::VDJ_REGION_TYPES,
        function(region_type) extract_sequence_region(AIRR_df, region_type)
    )
    for (regions in all_VDJ_regions) {
        expect_true(is.character(regions))
        expect_equal(length(regions), nrow(AIRR_df))
    }

    ## Extract regions as amino acid sequences.

    all_region_types <- c(igblastr:::FWRCDR_NAMES, igblastr:::VDJ_REGION_TYPES)
    all_aa_regions <- lapply(all_region_types,
        function(region_type) extract_sequence_region(AIRR_df, region_type,
                                                      as.aa=TRUE)
    )
    for (regions in all_aa_regions) {
        expect_true(is(regions, "AAStringSet"))
        expect_equal(length(regions), nrow(AIRR_df))
    }
})

