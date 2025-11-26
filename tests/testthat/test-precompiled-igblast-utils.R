test_that("inst/extdata/ncbi_igblast_data_files/ is not out-of-sync", {
    ## igblastr:::check_ncbi_igblast_data_files() is not guaranteed to work
    ## on Windows.
    if (.Platform$OS.type != "windows") {
        prev_warn <- getOption("warn")
        options(warn=2)
        on.exit(options(warn=prev_warn))
        igblastr:::check_ncbi_igblast_data_files(get_igblast_root())
    }
})

