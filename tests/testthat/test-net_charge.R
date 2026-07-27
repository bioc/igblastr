
test_that("net_charge()", {
    use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "heavy_sequences.fasta")
    AIRR_df <- igblastn(query)

    cdr3_aa <- extract_sequence_region(AIRR_df, region_type="cdr3", as.aa=TRUE)

    cdr3_charges <- net_charge(cdr3_aa, with.nTer=FALSE, with.cTer=FALSE)
    expect_true(is.numeric(cdr3_charges))
    expect_equal(length(cdr3_charges), length(cdr3_aa))
    expect_equal(cdr3_charges[[2L]], -1.189)

    pKa_values <- c(nTer=9.60, R=12.48, H=6.0, K=10.50,
                    D=3.86, E=4.25, C=8.33, Y=10.07, cTer=2.34)
    cdr3_charges <- net_charge(cdr3_aa, pKa_values=pKa_values)
    expect_true(is.numeric(cdr3_charges))
    expect_equal(length(cdr3_charges), length(cdr3_aa))
    expect_equal(cdr3_charges[[2L]], -1.181)

    cdr3_charges <- net_charge(cdr3_aa, with.nTer=FALSE, with.cTer=FALSE,
                               pKa_values=pKa_values)
    expect_true(is.numeric(cdr3_charges))
    expect_equal(length(cdr3_charges), length(cdr3_aa))
    expect_equal(cdr3_charges[[2L]], -1.175)

    fwr4_aa <- extract_sequence_region(AIRR_df, region_type="fwr4", as.aa=TRUE)

    fwr4_charges <- net_charge(fwr4_aa, as.matrix=TRUE)
    expect_true(is.matrix(fwr4_charges))
    expect_equal(nrow(fwr4_charges), length(fwr4_aa))
    expect_identical(rowSums(fwr4_charges), net_charge(fwr4_aa))
    expect_identical(colnames(fwr4_charges), names(PKA_VALUES))

    fwr4_charges <- net_charge(fwr4_aa, with.nTer=FALSE, as.matrix=TRUE)
    expect_true(is.matrix(fwr4_charges))
    expect_equal(nrow(fwr4_charges), length(fwr4_aa))
    expect_identical(rowSums(fwr4_charges),
                     net_charge(fwr4_aa, with.nTer=FALSE))
    expect_identical(colnames(fwr4_charges), tail(names(PKA_VALUES), n=-1))

    fwr4_charges <- net_charge(fwr4_aa, with.cTer=FALSE, as.matrix=TRUE)
    expect_true(is.matrix(fwr4_charges))
    expect_equal(nrow(fwr4_charges), length(fwr4_aa))
    expect_identical(rowSums(fwr4_charges),
                     net_charge(fwr4_aa, with.cTer=FALSE))
    expect_identical(colnames(fwr4_charges), head(names(PKA_VALUES), n=-1))
})

test_that("extract_region_net_charge()", {
    use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "light_sequences.fasta")
    AIRR_df <- igblastn(query)

    cdr2_aa <- extract_sequence_region(AIRR_df, region_type="cdr2", as.aa=TRUE)

    cdr2_charges <- extract_region_net_charge(AIRR_df, region_type="cdr2")
    expect_identical(cdr2_charges, net_charge(cdr2_aa,
                                              with.nTer=FALSE, with.cTer=FALSE))

    cdr2_charges <- extract_region_net_charge(AIRR_df, region_type="cdr2",
                                              pH=6.5)
    expect_identical(cdr2_charges, net_charge(cdr2_aa, pH=6.5,
                                              with.nTer=FALSE, with.cTer=FALSE))

    cdr2_charges <- extract_region_net_charge(AIRR_df, region_type="cdr2",
                                              as.matrix=TRUE)
    expect_identical(cdr2_charges, net_charge(cdr2_aa,
                                              with.nTer=FALSE, with.cTer=FALSE,
                                              as.matrix=TRUE))
})

