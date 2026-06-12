
test_that("install_custom_germline_db()", {

    db_name0 <- "_OGRDB.human.IGH+IGK+IGL.202605"
    alleles0 <- load_germline_sequences(db_name0)

    ## --- Use install_custom_germline_db() to create a germline db ---
    ## --- that is identical to _OGRDB.human.IGH+IGK+IGL.202605     ---

    db_name <- "cus.human.IGH+IGK+IGL"
    dir.create(fasta_dir <- tempfile())
    germline_sets <- c(IGH_VDJ=10, IGKappa_VJ=5, IGLambda_VJ=4)
    download_OGRDB_germline_sequences("Homo sapiens", germline_sets,
                                       destdir=fasta_dir)

    ## With auto-generated internal and auxiliary data.
    install_custom_germline_db(db_name, fasta_dir, gapped=TRUE,
                               intdata="auto", auxdata="auto")
    alleles <- load_germline_sequences(db_name)
    expect_identical(names(alleles), names(alleles0))
    expect_true(all(alleles == alleles0))
    expect_identical(load_intdata(db_name), load_intdata(db_name0))
    expect_identical(load_auxdata(db_name), load_auxdata(db_name0))
    rm_germline_db(db_name)

    ## With user-supplied internal and auxiliary data.
    my_intdata <- load_intdata(db_name0)
    my_auxdata <- load_auxdata(db_name0)
    install_custom_germline_db(db_name, fasta_dir, gapped=TRUE,
                               intdata=my_intdata, auxdata=my_auxdata)
    alleles <- load_germline_sequences(db_name)
    expect_identical(names(alleles), names(alleles0))
    expect_true(all(alleles == alleles0))
    expect_identical(load_intdata(db_name), my_intdata)
    expect_identical(load_auxdata(db_name), my_auxdata)
    rm_germline_db(db_name)
})

