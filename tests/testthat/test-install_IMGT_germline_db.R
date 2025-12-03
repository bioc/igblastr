test_that("install_IMGT_germline_db()", {
    db_name <- install_IMGT_germline_db("202531-1", "Homo sapiens",
                                        force=TRUE)
    expect_identical(db_name, "IMGT-202531-1.Homo_sapiens.IGH+IGK+IGL")
    use_germline_db(db_name)
    rm_germline_db(db_name)

    db_name <- install_IMGT_germline_db("202531-1", "Homo sapiens",
                                        tcr.db=TRUE, force=TRUE)
    expect_identical(db_name, "IMGT-202531-1.Homo_sapiens.TRA+TRB+TRG+TRD")
    use_germline_db(db_name)
    rm_germline_db(db_name)

    db_name <- install_IMGT_germline_db("202531-1", "Homo sapiens",
                                        loci="IGH", force=TRUE)
    expect_identical(db_name, "IMGT-202531-1.Homo_sapiens.IGH")
    use_germline_db(db_name)
    rm_germline_db(db_name)
 
    db_name <- install_IMGT_germline_db("202531-1", "Homo sapiens",
                                        loci=c("TRB", "TRA"), force=TRUE)
    expect_identical(db_name, "IMGT-202531-1.Homo_sapiens.TRA+TRB")
    use_germline_db(db_name)
    rm_germline_db(db_name)
})

test_that("install_IMGT_c_region_db()", {
    db_name <- install_IMGT_c_region_db("Homo sapiens", "TRB+TRA", force=TRUE)
    expect_identical(db_name, "IMGT.human.TRA+TRB.202509")
    use_c_region_db(db_name)
    rm_c_region_db(db_name)
})

