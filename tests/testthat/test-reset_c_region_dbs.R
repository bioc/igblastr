test_that(".preinstall_c_region_dbs()", {
    install_dir <- tempfile("preinstalled_c_region_dbs_")
    igblastr:::.preinstall_c_region_dbs(install_dir)
    db_list <- igblastr:::list_dbs(install_dir, what="C-region",
                                   long.listing=TRUE)

    expect_identical(db_list[["_IMGT.human.IGH+IGK+IGL.202605"]],
                     c(IGH=161L, IGK=5L, IGL=14L))

    expect_identical(db_list[["_IMGT.human.TRA+TRB+TRG+TRD.202605"]],
                     c(TRA=1L, TRB=6L, TRG=15L, TRD=1L))

    expect_identical(db_list[["_IMGT.mouse.IGH+IGK+IGL.202605"]],
                     c(IGH=59L, IGK=1L, IGL=4L))

    expect_identical(db_list[["_IMGT.rabbit.TRA+TRB+TRG+TRD.202605"]],
                     c(TRA=1L, TRB=3L, TRG=1L, TRD=1L))

    expect_identical(db_list[["_IMGT.rhesus_monkey.TRA+TRB+TRG+TRD.202605"]],
                     c(TRA=2L, TRB=3L, TRG=2L, TRD=3L))
})

