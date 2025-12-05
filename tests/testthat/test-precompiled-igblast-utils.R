test_that("inst/extdata/igdata_store/ is not outdated", {
    ## igblastr:::check_igdata_store() is not guaranteed to work on Windows.
    if (.Platform$OS.type != "windows") {
        prev_warn <- getOption("warn")
        options(warn=2)
        on.exit(options(warn=prev_warn))
        igblastr:::check_igdata_store(get_igblast_root())
    }
})

