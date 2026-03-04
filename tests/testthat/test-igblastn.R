test_that("igblastn(): basic operations", {
    use_germline_db("_OGRDB.human.IGH+IGK+IGL.202410")
    use_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "heavy_sequences.fasta")
    query <- readDNAStringSet(query)

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

test_that("igblastn(): germline_db_[VDJ]_seqidlist arguments", {
    use_germline_db("_OGRDB.human.IGH+IGK+IGL.202410")
    use_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "heavy_sequences.fasta")
    query10 <- head(readDNAStringSet(query), n=10L)

    temp_files0 <- list.files(tempdir(), all.files=TRUE, recursive=TRUE)
    out0 <- igblastn(query10)
    temp_files <- list.files(tempdir(), all.files=TRUE, recursive=TRUE)
    expect_identical(temp_files, temp_files0)

    ## Not really restricting anything here because these V, D, J alleles
    ## are those that got actually called by 'igblastn(query10)' above.
    V_seqidlist <- c("IGHV1-8*01",
                     "IGHV3-7*01",
                     "IGHV3-7*05",
                     "IGHV3-23*01",
                     "IGHV3-30*02",
                     "IGHV3-30*18",
                     "IGHV4-38-2*01",
                     "IGHV4-59*01")
    D_seqidlist <- c("IGHD2-2*01",
                     "IGHD2-2*02",
                     "IGHD2-8*01",
                     "IGHD2-8*02",
                     "IGHD3-10*01",
                     "IGHD3-10*03",
                     "IGHD3-16*03",
                     "IGHD4-4*01",
                     "IGHD4-23*01",
                     "IGHD6-6*01",
                     "IGHD6-13*01",
                     "IGHD6-19*01",
                     "IGHD6-25*01")
    J_seqidlist <- c("IGHJ2*01",
                     "IGHJ3*02",
                     "IGHJ4*02",
                     "IGHJ5*02")
    out1 <- suppressWarnings(
        igblastn(query10, germline_db_V_seqidlist=V_seqidlist,
                          germline_db_D_seqidlist=D_seqidlist,
                          germline_db_J_seqidlist=J_seqidlist)
    )
    temp_files <- list.files(tempdir(), all.files=TRUE, recursive=TRUE)
    expect_identical(temp_files, temp_files0)
    expect_identical(dim(out1), dim(out0))
    expect_identical(colnames(out1), colnames(out0))
    m <- match(paste0(c("v", "d", "j"), "_support"), colnames(out0))
    expect_identical(out1[ , -m], out0[ , -m])

    ## Let's restrict the search space to one allele per V/D gene.
    V_seqidlist <- c("IGHV1-8*01",
                     "IGHV3-7*01",
                     "IGHV3-23*01",
                     "IGHV3-30*02",
                     "IGHV4-59*01")
    D_seqidlist <- c("IGHD2-2*01",
                     "IGHD2-8*01",
                     "IGHD3-10*01",
                     "IGHD3-16*03",
                     "IGHD4-4*01",
                     "IGHD4-23*01",
                     "IGHD6-6*01",
                     "IGHD6-13*01",
                     "IGHD6-19*01",
                     "IGHD6-25*01")
    out2 <- suppressWarnings(
        igblastn(query10, germline_db_V_seqidlist=V_seqidlist,
                          germline_db_D_seqidlist=D_seqidlist)
    )
    temp_files <- list.files(tempdir(), all.files=TRUE, recursive=TRUE)
    expect_identical(temp_files, temp_files0)
    expect_identical(dim(out2), dim(out0))
    expect_true(all(unlist(strsplit(out2$v_call, ",")) %in% V_seqidlist))
    expect_true(all(unlist(strsplit(out2$d_call, ",")) %in% D_seqidlist))
})

