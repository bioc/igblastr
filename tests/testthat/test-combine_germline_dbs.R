.check_combined_germline_db <- function(db_name, db_name1, db_name2,
                                                 suffix1, suffix2)
{
    ## Check combined counts.
    all_counts <- list_germline_dbs(long.listing=TRUE)
    counts  <- all_counts[[db_name]]
    counts1 <- all_counts[[db_name1]]
    counts2 <- all_counts[[db_name2]]
    expect_identical(counts, counts1 + counts2)

    ## Check combined allele names.
    allele_names  <- names(load_germline_db(db_name))
    allele_names1 <- names(load_germline_db(db_name1))
    allele_names2 <- names(load_germline_db(db_name2))
    expect_identical(length(allele_names),  sum(counts))
    allele_names1 <- add_suffix(allele_names1, suffix1)
    allele_names2 <- add_suffix(allele_names2, suffix2)
    expect_true(setequal(allele_names, c(allele_names1, allele_names2)))

    ## Check combined intdata.
    intdata  <- load_intdata(db_name)
    intdata1 <- load_intdata(db_name1)
    intdata2 <- load_intdata(db_name2)
    intdata1$allele_name <- add_suffix(intdata1[ , "allele_name"], suffix1)
    intdata2$allele_name <- add_suffix(intdata2[ , "allele_name"], suffix2)
    expected <- rbind(intdata1, intdata2)
    expect_true(igblastr:::have_same_rows(intdata, expected, "allele_name"))

    ## Check combined auxdata.
    auxdata  <- load_auxdata(db_name)
    auxdata1 <- load_auxdata(db_name1)
    auxdata2 <- load_auxdata(db_name2)
    auxdata1$allele_name <- add_suffix(auxdata1[ , "allele_name"], suffix1)
    auxdata2$allele_name <- add_suffix(auxdata2[ , "allele_name"], suffix2)
    expected <- rbind(auxdata1, auxdata2)
    expect_true(igblastr:::have_same_rows(auxdata, expected, "allele_name"))
}

test_that("combine_germline_dbs()", {
    db_name1 <- "_OGRDB.human.IGH+IGK+IGL.202410"
    db_name2 <- "_OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202410"

    ## Combine 'db_name1' and 'db_name2'.
    db_name12 <- "comb.OGRDB.human+mouse.IGH+IGK+IGL"
    combine_germline_dbs(db_name12, db_name1, db_name2, "_Hs", "_Mm")
    .check_combined_germline_db(db_name12, db_name1, db_name2, "_Hs", "_Mm")

    ## Combine 'db_name2' and 'db_name1'.
    db_name21 <- "comb.OGRDB.mouse+human.IGH+IGK+IGL"
    combine_germline_dbs(db_name21, db_name2, db_name1, "_Mm", "_Hs")
    .check_combined_germline_db(db_name21, db_name2, db_name1, "_Mm", "_Hs")

    ## Use combined germline dbs 'db_name12' and 'db_name21' on a small
    ## set of human heavy chain BCR sequences.
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "heavy_sequences.fasta")
    use_germline_db(db_name12)
    AIRR_df12 <- igblastn(query, num_alignments_V=1,
                                 num_alignments_D=1,
                                 num_alignments_J=1)
    use_germline_db(db_name21)
    AIRR_df21 <- igblastn(query, num_alignments_V=1,
                                 num_alignments_D=1,
                                 num_alignments_J=1)

    ## The order in which we combined the germline dbs doesn't affect
    ## igblastn()'s results.
    expect_identical(AIRR_df12, AIRR_df21)

    expect_true(all(AIRR_df12$locus %in% "IGH"))
    expect_true(all(igblastr:::has_suffix(AIRR_df21$v_call, "_Hs")))
    ## Looks like 'AIRR_df21$d_call' contains some mouse D allele names!
    #expect_true(all(igblastr:::has_suffix(AIRR_df21$d_call, "_Hs")))
    expect_true(all(igblastr:::has_suffix(AIRR_df21$j_call, "_Hs")))

    ## Use combined germline dbs 'db_name12' and 'db_name21' on a small
    ## set of human light chain BCR sequences.
    query <- system.file(package="igblastr", "extdata",
                         "BCR", "light_sequences.fasta")
    use_germline_db(db_name12)
    AIRR_df12 <- igblastn(query, num_alignments_V=1,
                                 num_alignments_D=1,
                                 num_alignments_J=1)
    use_germline_db(db_name21)
    AIRR_df21 <- igblastn(query, num_alignments_V=1,
                                 num_alignments_D=1,
                                 num_alignments_J=1)

    ## The order in which we combined the germline dbs doesn't affect
    ## igblastn()'s results.
    expect_identical(AIRR_df12, AIRR_df21)

    expect_true(all(AIRR_df12$locus %in% c("IGK", "IGL")))
    expect_true(all(igblastr:::has_suffix(AIRR_df21$v_call, "_Hs")))
    expect_true(all(igblastr:::has_suffix(AIRR_df21$j_call, "_Hs")))

    rm_germline_db(db_name12)
    rm_germline_db(db_name21)
})

