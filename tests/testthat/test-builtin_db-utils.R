test_that("create_all_builtin_germline_dbs()", {
    destdir <- tempfile("builtin_germline_dbs_")
    igblastr:::create_all_builtin_germline_dbs(destdir)
    db_list <- igblastr:::list_dbs(destdir, what="germline", long.listing=TRUE)

    current <- db_list[["_AIRR.human.IGH+IGK+IGL.202410"]]
    expected <- rbind(IGH=c(V=198L, D=31L,  J=7L),
                      IGK=c(   64L,    0L,    7L),
                      IGL=c(   80L,    0L,    9L))
    expect_identical(current, expected)

    current <- db_list[["_AIRR.mouse.CAST_EiJ.IGH+IGK+IGL.202501"]]
    expected <- rbind(IGH=c(V=87L, D=9L,  J= 4L),
                      IGK=c(  88L,   0L,    11L),
                      IGL=c(   9L,   0L,     7L))
    expect_identical(current, expected)
})

test_that("create_all_builtin_c_region_dbs()", {
    destdir <- tempfile("builtin_c_region_dbs_")
    igblastr:::create_all_builtin_c_region_dbs(destdir)
    db_list <- igblastr:::list_dbs(destdir, what="C-region", long.listing=TRUE)

    expect_identical(db_list[["_IMGT.human.IGH+IGK+IGL.202412"]],
                     c(IGH=58L, IGK=5L, IGL=13L))

    expect_identical(db_list[["_IMGT.human.TRA+TRB+TRG+TRD.202509"]],
                     c(TRA=1L, TRB=3L, TRG=7L, TRD=1L))

    expect_identical(db_list[["_IMGT.mouse.IGH.202509"]],
                     c(IGH=56L))

    expect_identical(db_list[["_IMGT.mouse.TRA+TRB+TRG+TRD.202509"]],
                     c(TRA=2L, TRB=2L, TRG=4L, TRD=1L))
})

