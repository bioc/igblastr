test_that("igblastn()", {
    use_germline_db("_AIRR.human.IGH+IGK+IGL.202501")
    use_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")
    catnap_bnabs <- system.file(package="igblastr",
                                "extdata", "catnap_bnabs.fasta")
    query <- readDNAStringSet(catnap_bnabs)

    ## Call igblastn() on first 10 sequences.
    out0 <- igblastn(head(query, n=10L))
    expect_true(is_tibble(out0))
    expect_identical(dim(out0), c(10L, 111L))
    expect_identical(head(colnames(out0), n=2L), c("sequence_id", "sequence"))
    expect_identical(out0$sequence_id, head(names(query), n=10L))
    expect_identical(out0$sequence, unname(as.character(head(query, n=10L))))

    ## Call igblastn() on first 10 sequences using 1 thread.
    out1 <- igblastn(head(query, n=10L), num_threads=1)
    expect_identical(out1, out0)

    ## Call igblastn() on first 10 sequences using 5 threads.
    out5 <- igblastn(head(query, n=10L), num_threads=5)
    expect_identical(out5, out0)

    ## mclapply() with 'mc.cores' > 1 is not supported on Windows.
    if (.Platform$OS.type != "windows") {
        ## Call igblastn() on first 10 sequences in parallel using 4 workers.
        library(parallel)
        limit_cores <- isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_")))
        mc.cores <- if (limit_cores) 2L else 4L
        res <- mclapply(1:10, function(i) igblastn(query[i]), mc.cores=mc.cores)
        for (i in 1:10) {
            out_i <- res[[i]]
            expect_true(is_tibble(out_i))
            expect_identical(dim(out_i), c(1L, 111L))
            expect_identical(colnames(out_i), colnames(out0))
        }
        expect_identical(do.call(rbind, res), out0)

        ## Call igblastn() on first 5 sequences in parallel using 2 workers
        ## that try to write to the same output file.
        out <- tempfile()
        expect_warning(
            res <- mclapply(1:5, function(i) igblastn(query[i], out=out),
                            mc.cores=2),
            regexp="all scheduled cores encountered errors in user code"
        )
        is_error <- vapply(res, function(out) inherits(out, "try-error"),
                           logical(1))
        expect_true(all(is_error))
    }
})

