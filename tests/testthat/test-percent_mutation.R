
test_that("percent_mutation()", {
    use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "light_sequences.fasta")
    query <- head(readDNAStringSet(query), n=10L)
    AIRR_df <- igblastn(query)

    for (for.aa in c(FALSE, TRUE)) {
        percent_mut <- percent_mutation(AIRR_df, for.aa=for.aa)
        expect_true(is_tibble(percent_mut))
        expect_identical(dim(percent_mut), c(10L, 5L))
        expected_colnames <- c("sequence_id", "locus",
                               paste0(c("v", "d", "j"), "_perc_mut"))
        expect_identical(colnames(percent_mut), expected_colnames)
        expect_identical(percent_mut$sequence_id, AIRR_df$sequence_id)
        expect_identical(percent_mut$locus, AIRR_df$locus)
        expect_true(all(is.na(percent_mut$d_perc_mut)))
        m <- as.matrix(percent_mut[ , tail(expected_colnames, n=3L)])
        expect_true(typeof(m) == "double")
    }
})

