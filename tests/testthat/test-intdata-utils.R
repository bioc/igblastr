.IGBLAST_INTDATA_COLNAMES <- names(igblastr:::.IGBLAST_INTDATA_COL2CLASS)

.check_intdata <- function(organism, domain_system)
{
    intdata <- load_intdata(organism, domain_system=domain_system,
                            which="original")
    expect_true(is.data.frame(intdata))
    expect_identical(colnames(intdata), .IGBLAST_INTDATA_COLNAMES)

    allele_name <- intdata[ , "allele_name"]
    ## 13 rows are repeated in human.ndm.imgt
    ## 24 rows are repeated in mouse.ndm.imgt
    if (organism %in% c("human", "mouse") && domain_system == "imgt") {
        expect_true(rows_with_same_key_are_identical(intdata, "allele_name"))
    } else {
        expect_identical(anyDuplicated(allele_name), 0L)
    }

    fwr1_start <- intdata[ , "fwr1_start"]
    fwr1_end   <- intdata[ , "fwr1_end"]
    expect_true(all(fwr1_start == 1L))
    expect_true(all(fwr1_end >= fwr1_start + 2L))
    expect_true(all(fwr1_end %% 3L == 0L))

    cdr1_start <- intdata[ , "cdr1_start"]
    cdr1_end   <- intdata[ , "cdr1_end"]
    expect_true(all(cdr1_start == fwr1_end + 1L))
    expect_true(all(cdr1_end >= cdr1_start + 2L))
    expect_true(all(cdr1_end %% 3L == 0L))

    fwr2_start <- intdata[ , "fwr2_start"]
    fwr2_end   <- intdata[ , "fwr2_end"]
    expect_true(all(fwr2_start == cdr1_end + 1L))
    expect_true(all(fwr2_end >= fwr2_start + 2L))
    expect_true(all(fwr2_end %% 3L == 0L))

    cdr2_start <- intdata[ , "cdr2_start"]
    cdr2_end   <- intdata[ , "cdr2_end"]
    expect_true(all(cdr2_start == fwr2_end + 1L))
    expect_true(all(cdr2_end >= cdr2_start + 2L))
    if (organism == "mouse" && domain_system == "kabat") {
        ## In mouse.ndm.kabat, fwr3_end is not a multiple of 3 for
        ## allele "3609.1.84".
	expect_identical(allele_name[cdr2_end %% 3L != 0L], "3609.1.84")
    } else {
        expect_true(all(cdr2_end %% 3L == 0L))
    }

    fwr3_start <- intdata[ , "fwr3_start"]
    fwr3_end   <- intdata[ , "fwr3_end"]
    if (organism == "mouse" && domain_system == "imgt") {
        ## In mouse.ndm.imgt, cdr2_end = 201 and fwr3_start = 292 for
        ## allele "J558.1.85".
        expect_identical(allele_name[fwr3_start != cdr2_end + 1L], "J558.1.85")
    } else {
        expect_true(all(fwr3_start == cdr2_end + 1L))
    }
    expect_true(all(fwr3_end >= fwr3_start + 2L))
    if (organism == "human" && domain_system == "imgt") {
        ## In human.ndm.imgt, fwr3_end is not a multiple of 3 for
        ## allele "IGHV5-a*02".
        expect_identical(allele_name[fwr3_end %% 3L != 0L], "IGHV5-a*02")
    } else {
        expect_true(all(fwr3_end %% 3L == 0L))
    }

    expect_true(all(intdata[ , "coding_frame_start"] == 0L))

    intdata2 <- load_intdata(organism, for.aa=TRUE, domain_system=domain_system,
                             which="original")
    expect_true(is.data.frame(intdata2))
    expect_identical(colnames(intdata2), .IGBLAST_INTDATA_COLNAMES)
}

test_that("load_intdata()", {
    organisms <- list_igblast_organisms()
    for (organism in organisms)
        for (domain_system in c("imgt", "kabat"))
            .check_intdata(organism, domain_system)
})

