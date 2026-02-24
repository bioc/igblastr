test_that(".install_all_builtin_germline_dbs()", {
    destdir <- tempfile("builtin_germline_dbs_")
    igblastr:::.install_all_builtin_germline_dbs(destdir)
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

