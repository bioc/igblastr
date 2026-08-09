test_that(".preinstall_germline_dbs()", {
    install_dir <- tempfile("preinstalled_germline_dbs_")
    igblastr:::.preinstall_germline_dbs(install_dir)
    db_list <- igblastr:::list_dbs(install_dir, what="germline",
                                   long.listing=TRUE)

    current <- db_list[["_OGRDB.human.IGH+IGK+IGL.202410"]]
    expected <- rbind(IGH=c(V=198L, D=31L,  J=7L),
                      IGK=c(   64L,    0L,    7L),
                      IGL=c(   80L,    0L,    9L))
    expect_identical(current, expected)

    current <- db_list[["_OGRDB.human.IGH+IGK+IGL.202605"]]
    expected <- rbind(IGH=c(V=202L, D=31L,  J=7L),
                      IGK=c(   77L,    0L,    7L),
                      IGL=c(   88L,    0L,    9L))
    expect_identical(current, expected)

    current <- db_list[["_OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410"]]
    expected <- rbind(IGH=c(V=92L, D=10L,  J= 4L),
                      IGK=c(  89L,    0L,    11L),
                      IGL=c(   3L,    0L,     7L))
    expect_identical(current, expected)
})

