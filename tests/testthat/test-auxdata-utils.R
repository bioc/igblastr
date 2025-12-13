### Not the true colnames used in IgBLAST auxdata files.
### Ours are shorter and have underscores instead of spaces.
.IGBLAST_AUXDATA_COLNAMES <- c(
    "allele_name",
    "coding_frame_start",
    "chain_type",
    "cdr3_end",
    "extra_bps"
)

test_that("load_auxdata()", {
    organisms <- list_igblast_organisms()
    for (organism in organisms) {
        auxdata <- load_auxdata(organism, "original")
        expect_true(is.data.frame(auxdata))
        expect_identical(colnames(auxdata), .IGBLAST_AUXDATA_COLNAMES)
        ## human_gl.aux has 2 identical rows for TRAJ13*02 !
        if (organism != "human")
            expect_identical(anyDuplicated(auxdata[ , "allele_name"]), 0L)
    }
})

### Fix human aux data on-the-fly.
### We know that NCBI originally messed up with the 'extra_bps' value
### for alleles IGHJ6*02 and IGHJ6*03 in the original human_gl.aux. They
### corrected this later in the updated human_gl.aux that they released
### in April 2025. We do our own correction here.
.load_human_auxdata <- function()
{
    auxdata <- load_auxdata("human", "original")
    fixme <- auxdata[ , "allele_name"] %in% c("IGHJ6*02", "IGHJ6*03")
    auxdata[fixme, "extra_bps"] <- 1L  # (replace 0L with 1L)
    auxdata
}

### Install germline dbs used in tests below. Note that only the first
### installation actually triggers a download from IMGT. All subsequent
### installations obtain the data from the IMGT local store (located
### in 'igblastr_cache(IMGT_LOCAL_STORE)') so are very fast and work offline.
install_IMGT_germline_db("202531-1", "Homo sapiens", force=TRUE)
install_IMGT_germline_db("202531-1", "Mus musculus", force=TRUE)
install_IMGT_germline_db("202531-1", "Rattus norvegicus", force=TRUE)
install_IMGT_germline_db("202531-1", "Oryctolagus cuniculus", force=TRUE)

test_that("translate_J_alleles()", {

    ## --- for human J alleles (from AIRR and IMGT) ---

    auxdata <- .load_human_auxdata()

    db_name <- "_AIRR.human.IGH+IGK+IGL.202410"
    J_alleles <- load_germline_db(db_name, region_types="J")
    J_aa <- translate_J_alleles(J_alleles, auxdata)
    expect_true(is.character(J_aa))
    expect_identical(names(J_aa), names(J_alleles))
    expect_false(anyNA(J_aa))
    expect_identical(J_aa[["IGHJ1*01"]], "AEYFQHWGQGTLVTVSS")
    expect_true(all(grepl("[WF]G.G", J_aa)))

    db_name <- "IMGT-202531-1.Homo_sapiens.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    ## 2 human J alleles in IMGT release 202531-1 have no entries
    ## in 'auxdata'. Note that this could change in the future.
    allele_is_known <- names(J_alleles) %in% auxdata$allele_name
    expect_true(sum(!allele_is_known) <= 2L)
    J_aa <- translate_J_alleles(J_alleles, auxdata)
    expect_identical(unname(is.na(J_aa)), !allele_is_known)
    expect_identical(J_aa[["IGHJ1*01"]], "AEYFQHWGQGTLVTVSS")
    expect_true(all(grepl("[WF]G.G", J_aa[allele_is_known])))

    ## --- for mouse J alleles (from IMGT) ---

    auxdata <- load_auxdata("mouse", "original")

    db_name <- "IMGT-202531-1.Mus_musculus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    ## 3 mouse J alleles in IMGT release 202531-1 have no entries
    ## in 'auxdata'. Note that this could change in the future.
    allele_is_known <- names(J_alleles) %in% auxdata$allele_name
    expect_true(sum(!allele_is_known) <= 3L)
    ## Get rid of the "unknown" J alleles.
    known_J_alleles <- J_alleles[allele_is_known]
    J_aa <- translate_J_alleles(known_J_alleles, auxdata)
    expect_false(anyNA(J_aa))
    expect_identical(J_aa[["IGHJ1*01"]], "YWYFDVWGAGTTVTVSS")
})

test_that("J_allele_has_stop_codon()", {

    ## --- for human J alleles (from AIRR and IMGT) ---

    auxdata <- .load_human_auxdata()

    db_name <- "_AIRR.human.IGH+IGK+IGL.202410"
    J_alleles <- load_germline_db(db_name, region_types="J")
    has_stop_codon <- J_allele_has_stop_codon(J_alleles, auxdata)
    expect_true(is.logical(has_stop_codon))
    expect_identical(names(has_stop_codon), names(J_alleles))
    expect_false(any(has_stop_codon))

    db_name <- "IMGT-202531-1.Homo_sapiens.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    allele_is_known <- names(J_alleles) %in% auxdata$allele_name
    has_stop_codon <- J_allele_has_stop_codon(J_alleles, auxdata)
    expect_identical(unname(is.na(has_stop_codon)), !allele_is_known)
    expect_false(any(has_stop_codon[allele_is_known]))

    ## --- for rabbit J alleles (from IMGT) ---

    auxdata <- load_auxdata("rabbit", "original")

    db_name <- "IMGT-202531-1.Oryctolagus_cuniculus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    has_stop_codon <- J_allele_has_stop_codon(J_alleles, auxdata)
    expect_false(anyNA(has_stop_codon))
    expect_identical(names(J_alleles)[has_stop_codon], "IGKJ1-2*04")
})

test_that("translate_fwr4()", {

    ## --- for human J alleles (from AIRR and IMGT) ---

    auxdata <- .load_human_auxdata()

    db_name <- "_AIRR.human.IGH+IGK+IGL.202410"
    J_alleles <- load_germline_db(db_name, region_types="J")
    fwr4_aa <- translate_fwr4(J_alleles, auxdata)
    expect_true(is.character(fwr4_aa))
    expect_identical(names(fwr4_aa), names(J_alleles))
    expect_false(anyNA(fwr4_aa))
    fwr4_head <- translate_fwr4(J_alleles, auxdata, max_codons=4L)
    expect_true(is.character(fwr4_head))
    expect_identical(names(fwr4_head), names(J_alleles))
    expect_false(anyNA(fwr4_head))
    expect_true(all(nchar(fwr4_head) == 4L))
    expect_true(all(grepl("^[WF]G.G$", fwr4_head)))

    db_name <- "IMGT-202531-1.Homo_sapiens.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    allele_is_known <- names(J_alleles) %in% auxdata$allele_name
    fwr4_head <- translate_fwr4(J_alleles, auxdata, max_codons=4L)
    expect_identical(unname(is.na(fwr4_head)), !allele_is_known)
    expect_true(all(grepl("^[WF]G.G$", fwr4_head[allele_is_known])))

    ## --- for mouse J alleles (from IMGT) ---

    auxdata <- load_auxdata("mouse", "original")

    db_name <- "IMGT-202531-1.Mus_musculus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    ## 3 mouse J alleles in IMGT release 202531-1 have no entries
    ## in 'auxdata'. Note that this could change in the future.
    allele_is_known <- names(J_alleles) %in% auxdata$allele_name
    expect_true(sum(!allele_is_known) <= 3L)
    ## Get rid of the "unknown" J alleles.
    known_J_alleles <- J_alleles[allele_is_known]
    fwr4_head <- translate_fwr4(known_J_alleles, auxdata, max_codons=4L)
    expect_false(anyNA(fwr4_head))
    ## 3 "known" mouse J alleles in IMGT release 202531-1 don't have
    ## the expected motif at the beginning of their FWR4 region.
    ## Is this expected? Could this change in the future?
    surprise <- fwr4_head[!grepl("^[WF]G.G$", fwr4_head)]
    expect_identical(names(surprise), c("IGKJ3*01", "IGKJ3*02", "IGLJ3P*01"))
    expect_identical(unname(surprise), c("FSDG", "FSDG", "FSSN"))

    ## --- for rat J alleles (from IMGT) ---

    auxdata <- load_auxdata("rat", "original")

    db_name <- "IMGT-202531-1.Rattus_norvegicus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    fwr4_head <- translate_fwr4(J_alleles, auxdata, max_codons=4L)
    ## translate_fwr4() uses 'auxdata$cdr3_end' to get the position of
    ## the first FWR4 codons, but this column has an NA for IGKJ3*01.
    ## This could change in the future.
    ok <- !is.na(fwr4_head)
    expect_identical(names(fwr4_head)[!ok], "IGKJ3*01")
    ## Get rid of IGKJ3*01.
    fwr4_head <- fwr4_head[ok]
    ## 2 "known" rat J alleles in IMGT release 202531-1 don't have
    ## the expected motif at the beginning of their FWR4 region.
    ## Is this expected? Could this change in the future?
    surprise <- fwr4_head[!grepl("^[WF]G.G$", fwr4_head)]
    expect_identical(names(surprise), c("IGLJ2*01", "IGLJ4*01"))
    expect_identical(unname(surprise), c("LGKG", "LGKG"))
})

test_that("compute_auxdata()", {

    ## --- for human J alleles (from AIRR and IMGT) ---

    db_name <- "_AIRR.human.IGH+IGK+IGL.202410"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- compute_auxdata(J_alleles)
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .IGBLAST_AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Now we're going to check that 'computed_auxdata' agrees with the
    ## auxiliary data included in IgBLAST. More precisely, we're going to
    ## check that it's a subset of 'load_auxdata("human", "original")'.

    ## All the J alleles in _AIRR.human.IGH+IGK+IGL.202410 are annotated
    ## in human_gl.aux so we expect no NAs in 'm' below.
    orig_auxdata <- .load_human_auxdata()
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    expect_false(anyNA(m))
    orig_auxdata <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m)
    expect_identical(computed_auxdata, orig_auxdata)

    db_name <- "IMGT-202531-1.Homo_sapiens.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- compute_auxdata(J_alleles)
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .IGBLAST_AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202531-1.Homo_sapiens.IGH+IGK+IGL
    ## are annotated in human_gl.aux so we expect a few NAs in 'm' below.
    orig_auxdata <- .load_human_auxdata()
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    keep_idx <- which(!is.na(m))
    current <- S4Vectors:::extract_data_frame_rows(computed_auxdata, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m[keep_idx])
    expect_identical(current, target)

    ## --- for mouse J alleles (from IMGT) ---

    db_name <- "IMGT-202531-1.Mus_musculus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- suppressWarnings(compute_auxdata(J_alleles))
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .IGBLAST_AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202531-1.Mus_musculus.IGH+IGK+IGL
    ## are annotated in mouse_gl.aux so we expect a few NAs in 'm' below.
    ## We will also skip validation for alleles for which no CDR3 end was
    ## found.
    orig_auxdata <- load_auxdata("mouse", "original")
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    keep_idx <- which(!(is.na(computed_auxdata[ , "cdr3_end"]) | is.na(m)))
    current <- S4Vectors:::extract_data_frame_rows(computed_auxdata, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m[keep_idx])
    expect_identical(current, target)

    ## --- for rat J alleles (from IMGT) ---

    db_name <- "IMGT-202531-1.Rattus_norvegicus.IGH+IGK+IGL"
    J_alleles <- load_germline_db(db_name, region_types="J")
    computed_auxdata <- suppressWarnings(compute_auxdata(J_alleles))
    expect_true(is.data.frame(computed_auxdata))
    expect_identical(colnames(computed_auxdata), .IGBLAST_AUXDATA_COLNAMES)
    expect_identical(computed_auxdata[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202531-1.Mus_musculus.IGH+IGK+IGL
    ## are annotated in rat_gl.aux so we expect a few NAs in 'm' below.
    ## We will also skip validation for alleles for which no CDR3 end was
    ## found.
    orig_auxdata <- load_auxdata("rat", "original")
    m <- match(names(J_alleles), orig_auxdata[ , "allele_name"])
    keep_idx <- which(!(is.na(computed_auxdata[ , "cdr3_end"]) | is.na(m)))
    current <- S4Vectors:::extract_data_frame_rows(computed_auxdata, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_auxdata, m[keep_idx])
    expect_identical(current, target)
})

