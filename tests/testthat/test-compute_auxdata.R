### Install IMGT germline dbs used in tests below. Note that only the first
### installation actually triggers a download from IMGT. All subsequent
### installations obtain the data from the IMGT local store (located
### in 'igblastr_cache(IMGT_STORE)') so are very fast and work offline.
install_IMGT_germline_db("202614-2", "Homo sapiens",
                         without.auxdata=TRUE, overwrite=TRUE)
install_IMGT_germline_db("202614-2", "Mus musculus",
                         without.auxdata=TRUE, overwrite=TRUE)
install_IMGT_germline_db("202614-2", "Rattus norvegicus",
                         without.auxdata=TRUE, overwrite=TRUE)

.AUXDATA_COLNAMES <- names(igblastr:::AUXDATA_COL2CLASS)

test_that("compute_auxdata()", {

    ## --- for human J alleles (from AIRR and IMGT) ---

    db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- compute_auxdata(J_alleles)
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Now we're going to check that 'computed_auxdata' agrees with the
    ## auxiliary data included in IgBLAST. More precisely, we're going to
    ## check that it's a subset of 'load_auxdata("human", which="original")'.

    ## All the J alleles in _OGRDB.human.IGH+IGK+IGL.202605 are annotated
    ## in human_gl.aux so we expect no NAs in 'm' below.
    orig_auxdata <- load_and_fix_human_auxdata()
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    expect_false(anyNA(m))
    orig_auxdata <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m)
    expect_identical(computed_auxdata, orig_auxdata)

    db_name <- "IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- compute_auxdata(J_alleles)
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL
    ## are annotated in human_gl.aux so we expect a few NAs in 'm' below.
    orig_auxdata <- load_and_fix_human_auxdata()
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    keep_idx <- which(!is.na(m))
    current <- S4Vectors:::extract_data_frame_rows(computed_auxdata, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m[keep_idx])
    expect_identical(current, target)

    ## --- for mouse J alleles (from IMGT) ---

    db_name <- "IMGT-202614-2.Mus_musculus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- suppressWarnings(compute_auxdata(J_alleles))
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202614-2.Mus_musculus.IGH+IGK+IGL
    ## are annotated in mouse_gl.aux so we expect a few NAs in 'm' below.
    ## We will also skip validation for alleles for which no CDR3 end was
    ## found.
    orig_auxdata <- load_auxdata("mouse", which="original")
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    keep_idx <- which(!(is.na(computed_auxdata[ , "cdr3_end"]) | is.na(m)))
    current <- S4Vectors:::extract_data_frame_rows(computed_auxdata, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m[keep_idx])
    expect_identical(current, target)

    ## --- for rat J alleles (from IMGT) ---

    db_name <- "IMGT-202614-2.Rattus_norvegicus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- suppressWarnings(compute_auxdata(J_alleles))
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202614-2.Mus_musculus.IGH+IGK+IGL
    ## are annotated in rat_gl.aux so we expect a few NAs in 'm' below.
    ## We will also skip validation for alleles for which no CDR3 end was
    ## found.
    orig_auxdata <- load_auxdata("rat", which="original")
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    keep_idx <- which(!(is.na(computed_auxdata[ , "cdr3_end"]) | is.na(m)))
    current <- S4Vectors:::extract_data_frame_rows(computed_auxdata, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m[keep_idx])
    expect_identical(current, target)
})

test_that("find_discordant_auxdata()/complete_auxdata()", {
    find_discordant_auxdata <- igblastr:::find_discordant_auxdata
    complete_auxdata <- igblastr:::complete_auxdata

    ## --- trivial case ---

    auxdata1 <- data.frame(
        allele_name=c("A", "B", "C", "D", "E", "F"),
        coding_frame_start=c(1L, 2L, 1L, 0L, 1L, 1L),
        chain_type=c("JH", "JH", "JK", "JL", "JH", "JK"),
        cdr3_end=c(11:13, NA, 15:16),
        extra_bps=1L
    )
    expect_identical(complete_auxdata(auxdata1, auxdata1), auxdata1)

    ## --- discordant data ---

    auxdata2 <- data.frame(
        allele_name=c("E", "F", "G", "H", "B", "D", "I", "A"),
        coding_frame_start=c(NA, 1L, 1L, 2L, 1L, 0L, 0L, 1L),
        chain_type=c("JH", "JK", "JH", "JH", "JH", "JL", "JL", "JH"),
        cdr3_end=c(11L, 11:13, 12L, 15:16, 11L),
        extra_bps=1L
    )

    expect_error(complete_auxdata(auxdata1, auxdata2), regexp="discordant")
    disc_rowpairs <- find_discordant_auxdata(auxdata1, auxdata2)
    expected <- data.frame(rowidx1=c(2L, 5L, 6L), rowidx2=c(5L, 1L, 2L))
    expect_identical(disc_rowpairs, expected)

    nopairs <- data.frame(rowidx1=integer(0), rowidx2=integer(0))
    disc_rowpairs <- find_discordant_auxdata(auxdata1,
                                             auxdata2[-disc_rowpairs[[2L]], ])
    expect_identical(disc_rowpairs, nopairs)
    disc_rowpairs <- find_discordant_auxdata(auxdata1[-disc_rowpairs[[1L]], ],
                                             auxdata2)
    expect_identical(disc_rowpairs, nopairs)

    ## Swapping the two data.frames returns the same pairs but they are
    ## possibly in a different order.
    expect_error(complete_auxdata(auxdata2, auxdata1), regexp="discordant")
    disc_rowpairs <- find_discordant_auxdata(auxdata2, auxdata1)
    expected <- data.frame(rowidx1=c(1L, 2L, 5L), rowidx2=c(5L, 6L, 2L))
    expect_identical(disc_rowpairs, expected)

    ## --- concordant data ---

    auxdata1[c(2L, 5L), "coding_frame_start"] <- NA
    auxdata2[1:2, "cdr3_end"] <- NA

    disc_rowpairs <- find_discordant_auxdata(auxdata1, auxdata2)
    expect_identical(disc_rowpairs, nopairs)
    expected <- auxdata1
    expected[2L, "coding_frame_start"] <- 1L
    expected[4L, "cdr3_end"] <- 15L
    auxdata1a <- complete_auxdata(auxdata1, auxdata2)
    expect_identical(auxdata1a, expected)
    expect_identical(complete_auxdata(auxdata1a, auxdata2), auxdata1a)

    disc_rowpairs <- find_discordant_auxdata(auxdata2, auxdata1)
    expect_identical(disc_rowpairs, nopairs)
    expected <- auxdata2
    expected[c(1L, 2L), "cdr3_end"] <- c(15L, 16L)
    auxdata2a <- complete_auxdata(auxdata2, auxdata1)
    expect_identical(auxdata2a, expected)
    expect_identical(complete_auxdata(auxdata2a, auxdata1), auxdata2a)
})

