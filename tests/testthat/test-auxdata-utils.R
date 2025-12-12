### Note really the true colnames used in IgBLAST aux files.
### Ours are shorter and have underscores instead of spaces.
.IGBLAST_AUXDATA_COLNAMES <- c(
    "allele_name",
    "coding_frame_start",
    "chain_type",
    "CDR3_stop",
    "extra_bps"
)

test_that("load_auxdata()", {
    organisms <- list_igblast_organisms()
    for (organism in organisms) {
        aux <- load_auxdata(organism)
        expect_true(is.data.frame(aux))
        expect_identical(colnames(aux), .IGBLAST_AUXDATA_COLNAMES)
    }
})

### Fix human aux data on-the-fly.
### We know that NCBI originally messed up with the 'extra_bps' value
### for alleles IGHJ6*02 and IGHJ6*03 in the original human_gl.aux. They
### corrected this later in the updated human_gl.aux that they released
### in April 2025. We do our own correction here.
.load_human_auxdata <- function()
{
    human_aux <- load_auxdata("human", "original")
    fixme <- human_aux[ , "allele_name"] %in% c("IGHJ6*02", "IGHJ6*03")
    human_aux[fixme, "extra_bps"] <- 1L  # (replace 0L with 1L)
    human_aux
}

test_that("compute_auxdata()", {
    EXTENDED_AUXDATA_COLNAMES <- c(.IGBLAST_AUXDATA_COLNAMES,
                                   #"FWR4_start_pattern",
                                   "has_stop_codon")

    ## --- on human V alleles ---

    db_name <- "_AIRR.human.IGH+IGK+IGL.202410"
    J_alleles <- load_germline_db(db_name, region_types="J")
    human_aux <- compute_auxdata(J_alleles)
    expect_true(is.data.frame(human_aux))
    expect_identical(colnames(human_aux), EXTENDED_AUXDATA_COLNAMES)
    expect_identical(human_aux[ , "allele_name"], names(J_alleles))

    ## Now we're going to check that 'human_aux' agrees with the auxiliary
    ## data included in IgBLAST. More precisely, we're going to check that
    ## it's a subset of 'load_auxdata("human", "original")'.

    ## All the J alleles in _AIRR.human.IGH+IGK+IGL.202410 are annotated
    ## in human_gl.aux so we expect no NAs in 'm' below.
    orig_human_aux <- .load_human_auxdata()
    m <- match(names(J_alleles), orig_human_aux[ , "allele_name"])
    expect_false(anyNA(m))
    orig_human_aux <- S4Vectors:::extract_data_frame_rows(orig_human_aux, m)
    expect_identical(human_aux[ , .IGBLAST_AUXDATA_COLNAMES],
                     orig_human_aux)

    db_name <- install_IMGT_germline_db("202531-1", "Homo sapiens",
                                        force=TRUE)
    J_alleles <- load_germline_db(db_name, region_types="J")
    human_aux <- compute_auxdata(J_alleles)
    expect_true(is.data.frame(human_aux))
    expect_identical(colnames(human_aux), EXTENDED_AUXDATA_COLNAMES)
    expect_identical(human_aux[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202531-1.Homo_sapiens.IGH+IGK+IGL
    ## are annotated in human_gl.aux so we expect a few NAs in 'm' below.
    orig_human_aux <- .load_human_auxdata()
    m <- match(names(J_alleles), orig_human_aux[ , "allele_name"])
    keep_idx <- which(!is.na(m))
    current <- S4Vectors:::extract_data_frame_rows(human_aux, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_human_aux, m[keep_idx])
    expect_identical(current[ , .IGBLAST_AUXDATA_COLNAMES], target)

    ## --- on mouse V alleles ---

    db_name <- install_IMGT_germline_db("202531-1", "Mus musculus",
                                        force=TRUE)
    J_alleles <- load_germline_db(db_name, region_types="J")
    mouse_aux <- suppressWarnings(compute_auxdata(J_alleles))
    expect_true(is.data.frame(mouse_aux))
    expect_identical(colnames(mouse_aux), EXTENDED_AUXDATA_COLNAMES)
    expect_identical(mouse_aux[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202531-1.Mus_musculus.IGH+IGK+IGL
    ## are annotated in mouse_gl.aux so we expect a few NAs in 'm' below.
    ## We will also skip validation for alleles for which no CDR3 stop was
    ## found.
    orig_mouse_aux <- load_auxdata("mouse", "original")
    m <- match(names(J_alleles), orig_mouse_aux[ , "allele_name"])
    keep_idx <- which(!(is.na(mouse_aux[ , "CDR3_stop"]) | is.na(m)))
    current <- S4Vectors:::extract_data_frame_rows(mouse_aux, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_mouse_aux, m[keep_idx])
    expect_identical(current[ , .IGBLAST_AUXDATA_COLNAMES], target)

    ## --- on rat V alleles ---

    db_name <- install_IMGT_germline_db("202531-1", "Rattus norvegicus",
                                        force=TRUE)
    J_alleles <- load_germline_db(db_name, region_types="J")
    rat_aux <- suppressWarnings(compute_auxdata(J_alleles))
    expect_true(is.data.frame(rat_aux))
    expect_identical(colnames(rat_aux), EXTENDED_AUXDATA_COLNAMES)
    expect_identical(rat_aux[ , "allele_name"], names(J_alleles))

    ## Not all the J alleles in IMGT-202531-1.Mus_musculus.IGH+IGK+IGL
    ## are annotated in rat_gl.aux so we expect a few NAs in 'm' below.
    ## We will also skip validation for alleles for which no CDR3 stop was
    ## found.
    orig_rat_aux <- load_auxdata("rat", "original")
    m <- match(names(J_alleles), orig_rat_aux[ , "allele_name"])
    keep_idx <- which(!(is.na(rat_aux[ , "CDR3_stop"]) | is.na(m)))
    current <- S4Vectors:::extract_data_frame_rows(rat_aux, keep_idx)
    target <- S4Vectors:::extract_data_frame_rows(orig_rat_aux, m[keep_idx])
    expect_identical(current[ , .IGBLAST_AUXDATA_COLNAMES], target)
})

