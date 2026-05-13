test_that(".create_all_builtin_c_region_dbs()", {
    destdir <- tempfile("builtin_c_region_dbs_")
    igblastr:::.create_all_builtin_c_region_dbs(destdir)
    db_list <- igblastr:::list_dbs(destdir, what="C-region", long.listing=TRUE)

    expect_identical(db_list[["_IMGT.human.IGH+IGK+IGL.202412"]],
                     c(IGH=58L, IGK=5L, IGL=13L))

    expect_identical(db_list[["_IMGT.human.TRA+TRB+TRG+TRD.202509"]],
                     c(TRA=1L, TRB=3L, TRG=7L, TRD=1L))

    expect_identical(db_list[["_IMGT.mouse.IGH.202509"]],
                     c(IGH=55L))

    expect_identical(db_list[["_IMGT.mouse.TRA+TRB+TRG+TRD.202509"]],
                     c(TRA=2L, TRB=2L, TRG=4L, TRD=1L))
})

