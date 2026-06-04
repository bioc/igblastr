
test_that("percent_mutation()", {
    use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "light_sequences.fasta")
    query <- head(readDNAStringSet(query), n=10L)
    AIRR_df <- igblastn(query)

    EXPECTED_COLNAMES <- c("sequence_id", "locus",
                           paste0(c("v", "d", "j"), "_perc_mut"))
    for (for.aa in c(FALSE, TRUE)) {
        percent_mut <- percent_mutation(AIRR_df, for.aa=for.aa)
        expect_true(is_tibble(percent_mut))
        expect_identical(dim(percent_mut), c(10L, 5L))
        expected_colnames <- EXPECTED_COLNAMES
        if (for.aa)
            expected_colnames[3:5] <- paste0(expected_colnames[3:5], "_aa")
        expect_identical(colnames(percent_mut), expected_colnames)
        expect_identical(percent_mut$sequence_id, AIRR_df$sequence_id)
        expect_identical(percent_mut$locus, AIRR_df$locus)
        ## Because all the queries are from the light chain, then
        ## 'AIRR_df$d_call' and all the 'd_*' columns in 'AIRR_df' are
        ## filled with NAs. So the 'd_perc_mut(_aa)' column in 'percent_mut'
        ## should also be filled with NAs.
        expect_true(all(is.na(percent_mut[[4L]])))  # d_perc_mut(_aa) col
        m <- as.matrix(percent_mut[ , tail(expected_colnames, n=3L)])
        expect_true(typeof(m) == "double")
    }
})

