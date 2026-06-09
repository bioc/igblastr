
test_that("AIRR output columns are valid Rearrangement fields", {
    library(airr)
    Rearrangement_fields <- c(RearrangementSchema@required,
                              RearrangementSchema@optional)
    use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")
    use_c_region_db("_IMGT.human.IGH+IGK+IGL.202605")
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "light_sequences.fasta")
    query <- head(readDNAStringSet(query), n=10L)
    AIRR_df <- igblastn(query)
    expect_true(all(colnames(AIRR_df) %in% Rearrangement_fields))
})

