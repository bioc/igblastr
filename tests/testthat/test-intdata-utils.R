
test_that("load_intdata()", {
    organisms <- list_igblast_organisms()
    for (organism in organisms)
        for (domain_system in c("imgt", "kabat")) {
            intdata <- load_intdata(organism, domain_system=domain_system,
                                    which="original")
            ## 13 rows are repeated in original human.ndm.imgt
            ## 24 rows are repeated in original mouse.ndm.imgt
            ## --> we set 'allow.dup.entries' to TRUE for these two datasets.
            allow.dup.entries <- (organism %in% c("human", "mouse")) &&
                                 (domain_system == "imgt")
            ok <- check_V_ndm_data(intdata, allow.dup.entries=allow.dup.entries)
            if (organism == "human" && domain_system == "imgt") {
                ## We expect all *_end columns to contain multiples of 3.
                ## However, in original human.ndm.imgt, fwr3_end is not a
                ## multiple of 3 for allele "IGHV5-a*02".
                expect_identical(intdata[!ok, "allele_name"], "IGHV5-a*02")
            } else if (organism == "mouse" && domain_system == "imgt") {
                ## We expect fwr3_start == cdr2_end + 1. However, in
                ## original mouse.ndm.imgt, cdr2_end = 201 and fwr3_start = 292
                ## for allele "J558.1.85".
                expect_identical(intdata[!ok, "allele_name"], "J558.1.85")
            } else if (organism == "mouse" && domain_system == "kabat") {
                ## We expect all *_end columns to contain multiples of 3.
                ## However, in mouse.ndm.kabat, fwr3_end is not a multiple
                ## of 3 for allele "3609.1.84".
                expect_identical(intdata[!ok, "allele_name"], "3609.1.84")
            } else {
                expect_true(all(ok))
            }

            ## Very loose check of the internal *.pdm.* files.
            ## Note that, assuming these files are used by igblastp, which
            ## we don't support yet, we don't really care about these files
            ## for now.
            intdata2 <- load_intdata(organism, for.aa=TRUE,
                                     domain_system=domain_system,
                                     which="original")
            expect_true(is.data.frame(intdata2))
            expected_colnames <- names(igblastr:::.IGBLAST_INTDATA_COL2CLASS)
            expect_identical(colnames(intdata2), expected_colnames)
        }
})

