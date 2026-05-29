test_that("igblastn(): basic operations", {
    use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")
    use_c_region_db("_IMGT.human.IGH+IGK+IGL.202605")
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "heavy_sequences.fasta")
    query <- readDNAStringSet(query)

    ## Call igblastn() on first 10 sequences.
    AIRR_df0 <- igblastn(head(query, n=10L))
    expect_true(is_tibble(AIRR_df0))
    expect_identical(dim(AIRR_df0), c(10L, 111L))
    expect_identical(head(colnames(AIRR_df0), n=2L),
                     c("sequence_id", "sequence"))
    expect_identical(AIRR_df0$sequence_id, head(names(query), n=10L))
    expect_identical(AIRR_df0$sequence,
                     unname(as.character(head(query, n=10L))))

    ## Call igblastn() on first 10 sequences using 1 thread.
    AIRR_df1 <- igblastn(head(query, n=10L), num_threads=1)
    expect_identical(AIRR_df1, AIRR_df0)

    ## Call igblastn() on first 10 sequences using 5 threads.
    AIRR_df5 <- igblastn(head(query, n=10L), num_threads=5)
    expect_identical(AIRR_df5, AIRR_df0)

    ## mclapply() with 'mc.cores' > 1 is not supported on Windows.
    if (.Platform$OS.type != "windows") {
        ## Call igblastn() on first 10 sequences in parallel using 4 workers.
        library(parallel)
        limit_cores <- isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_")))
        mc.cores <- if (limit_cores) 2L else 4L
        AIRR_dfs <- mclapply(1:10,
                             function(i) igblastn(query[i]), mc.cores=mc.cores)
        for (i in 1:10) {
            AIRR_df <- AIRR_dfs[[i]]
            expect_true(is_tibble(AIRR_df))
            expect_identical(dim(AIRR_df), c(1L, 111L))
            expect_identical(colnames(AIRR_df), colnames(AIRR_df0))
        }
        expect_identical(do.call(rbind, AIRR_dfs), AIRR_df0)

        ## Call igblastn() on first 5 sequences in parallel using 2 workers
        ## that try to write to the same output file.
        out <- tempfile()
        expect_warning(
            AIRR_dfs <- mclapply(1:5,
                                 function(i) igblastn(query[i], out=out),
                                 mc.cores=2),
            regexp="all scheduled cores encountered errors in user code"
        )
        is_error <- vapply(AIRR_dfs,
                           function(AIRR_df) inherits(AIRR_df, "try-error"),
                           logical(1))
        expect_true(all(is_error))
    }
})

test_that("igblastn(): 'germline_db_[VDJ]_seqidlist' arguments", {
    use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")
    use_c_region_db("_IMGT.human.IGH+IGK+IGL.202605")
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "heavy_sequences.fasta")
    query10 <- head(readDNAStringSet(query), n=10L)

    temp_files0 <- list.files(tempdir(), all.files=TRUE, recursive=TRUE)
    AIRR_df0 <- igblastn(query10)
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
    AIRR_df1 <- suppressWarnings(
        igblastn(query10, germline_db_V_seqidlist=V_seqidlist,
                          germline_db_D_seqidlist=D_seqidlist,
                          germline_db_J_seqidlist=J_seqidlist)
    )
    temp_files <- list.files(tempdir(), all.files=TRUE, recursive=TRUE)
    expect_identical(temp_files, temp_files0)
    expect_identical(dim(AIRR_df1), dim(AIRR_df0))
    expect_identical(colnames(AIRR_df1), colnames(AIRR_df0))
    m <- match(paste0(c("v", "d", "j"), "_support"), colnames(AIRR_df0))
    expect_identical(AIRR_df1[ , -m], AIRR_df0[ , -m])

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
    AIRR_df2 <- suppressWarnings(
        igblastn(query10, germline_db_V_seqidlist=V_seqidlist,
                          germline_db_D_seqidlist=D_seqidlist)
    )
    temp_files <- list.files(tempdir(), all.files=TRUE, recursive=TRUE)
    expect_identical(temp_files, temp_files0)
    expect_identical(dim(AIRR_df2), dim(AIRR_df0))
    expect_true(all(unlist(strsplit(AIRR_df2$v_call, ",")) %in% V_seqidlist))
    expect_true(all(unlist(strsplit(AIRR_df2$d_call, ",")) %in% D_seqidlist))
})

test_that("igblastn(): 'auxiliary_data' argument", {
    db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"
    use_germline_db(db_name)
    use_c_region_db("")
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "heavy_sequences.fasta")
    query10 <- head(readDNAStringSet(query), n=10L)
    AIRR_df0 <- igblastn(query10, num_alignments_J=1)
    auxdata <- load_auxdata(db_name)
    AIRR_df <- igblastn(query10, auxiliary_data=auxdata, num_alignments_J=1)
    expect_identical(AIRR_df, AIRR_df0)

    ## --- With some J alleles with a missing "cdr3_end" in 'auxdata' ---

    ## Columns Of Interest. These are the columns that will get NAs for
    ## queries that are mapped to a J allele (via v_call) with an unknown
    ## CDR3 end position.
    COI1 <- c("fwr4", "fwr4_aa", "cdr3", "cdr3_aa",
              "junction", "junction_length",
              "junction_aa", "junction_aa_length",
              "fwr4_start", "fwr4_end", "cdr3_start", "cdr3_end")
    expect_false(anyNA(AIRR_df0[ , COI1]))

    not_called <- c("IGHJ6*03", "IGLJ7*02")  # not in 'AIRR_df0$j_call'
    auxdata[auxdata$allele_name %in% not_called, "cdr3_end"] <- NA
    AIRR_df <- igblastn(query10, auxiliary_data=auxdata, num_alignments_J=1)
    expect_identical(AIRR_df, AIRR_df0)

    called <- c("IGHJ2*01", "IGHJ4*02")      # in 'AIRR_df0$j_call'
    auxdata[auxdata$allele_name %in% called, "cdr3_end"] <- NA
    AIRR_df <- igblastn(query10, auxiliary_data=auxdata, num_alignments_J=1)
    expect_identical(AIRR_df[ , "j_call"], AIRR_df0[ , "j_call"])
    ROI <- which(AIRR_df$j_call %in% called)  # Rows Of Interest
    expect_true(all(is.na(AIRR_df[ROI, COI1])))
    AIRR_df[ROI, COI1] <- AIRR_df0[ROI, COI1]
    expect_identical(AIRR_df, AIRR_df0)

    ## --- With some J alleles with missing "coding_frame_start", ---
    ## ---        "cdr3_end", and "extra_bps" in 'auxdata'        ---

    COI2 <- c("vj_in_frame", "v_frameshift", "productive", COI1)
    expect_false(anyNA(AIRR_df0[ , COI2]))

    auxdata <- load_auxdata(db_name)
    cols <- c("coding_frame_start", "cdr3_end", "extra_bps")
    called <- "IGHJ4*02"  # in 'AIRR_df0$j_call'
    auxdata[auxdata$allele_name %in% called, cols] <- NA
    AIRR_df1 <- igblastn(query10, auxiliary_data=auxdata, num_alignments_J=1)
    expect_identical(AIRR_df1[ , "j_call"], AIRR_df0[ , "j_call"])
    ROI <- which(AIRR_df1$j_call %in% called)  # Rows Of Interest
    expect_true(all(is.na(AIRR_df1[ROI, COI2])))
    AIRR_df <- AIRR_df1
    AIRR_df[ROI, COI2] <- AIRR_df0[ROI, COI2]
    expect_identical(AIRR_df, AIRR_df0)

    ## ---        With some J alleles not annotated in 'auxdata'.       ---
    ## --- Note that this is equivalent to having "coding_frame_start", ---
    ## ----   "cdr3_end", and "extra_bps" set to NA for these alleles.  ---

    auxdata <- load_auxdata(db_name)
    not_called <- c("IGHJ6*03", "IGLJ7*02")  # not in 'AIRR_df0$j_call'
    auxdata <- subset(auxdata, allele_name %notin% not_called)
    ## Supplying incomplete auxiliary data triggers a warning.
    expect_warning(
      AIRR_df <- igblastn(query10, auxiliary_data=auxdata, num_alignments_J=1),
      regexp="Incomplete auxiliary data"
    )
    expect_identical(AIRR_df, AIRR_df0)

    called <- "IGHJ4*02"  # in 'AIRR_df0$j_call'
    auxdata <- subset(auxdata, allele_name %notin% called)
    ## Supplying incomplete auxiliary data triggers a warning.
    expect_warning(
      AIRR_df <- igblastn(query10, auxiliary_data=auxdata, num_alignments_J=1),
      regexp="Incomplete auxiliary data"
    )
    expect_identical(AIRR_df, AIRR_df1)
})

