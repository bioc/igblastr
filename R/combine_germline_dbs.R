### =========================================================================
### combine_germline_dbs()
### -------------------------------------------------------------------------
###


.combine_locus_alleles <- function(locus, alleles1, alleles2, suffix1, suffix2)
{
    stopifnot(isSingleNonWhiteString(locus),
              is(alleles1, "DNAStringSet"), is(alleles2, "DNAStringSet"),
              isSingleNonWhiteString(suffix1), isSingleNonWhiteString(suffix2))
    locus_alleles1 <- alleles1[substr(names(alleles1), 1L, 3L) == locus]
    locus_alleles2 <- alleles2[substr(names(alleles2), 1L, 3L) == locus]
    ## Add user-supplied suffixes to disambiguate allele names.
    names(locus_alleles1) <- add_suffix(names(locus_alleles1), suffix1)
    names(locus_alleles2) <- add_suffix(names(locus_alleles2), suffix2)
    c(locus_alleles1, locus_alleles2)
}

.combine_region_dbs <- function(destdir, db_name1, db_name2, suffix1, suffix2,
                                region_type, verbose=FALSE)
{
    alleles1 <- load_germline_db(db_name1, region_types=region_type)
    alleles2 <- load_germline_db(db_name2, region_types=region_type)
    loci1 <- substr(names(alleles1), 1, 3)
    loci2 <- substr(names(alleles2), 1, 3)
    stopifnot(all(loci1 %in% IG_LOCI), all(loci2 %in% IG_LOCI))

    locus2alleles <- lapply(
        setNames(IG_LOCI, IG_LOCI),
        .combine_locus_alleles, alleles1, alleles2, suffix1, suffix2
    )

    locus2alleles <- locus2alleles[lengths(locus2alleles) != 0L]
    stopifnot(length(locus2alleles) != 0L)
    names(locus2alleles) <- paste0(names(locus2alleles), region_type, ".fasta")

    fasta_files <- vapply(locus2alleles,
        function(alleles) {
            fasta_file <- tempfile(fileext=".fasta")
            writeXStringSet(alleles, fasta_file)
            fasta_file
        },
        character(1))

    create_region_db(fasta_files, destdir, region_type=region_type,
                     verbose=verbose)
    lapply(fasta_files, nuke_file)
}

.combine_db_intdata <- function(destdir, db_name1, db_name2, suffix1, suffix2,
                                verbose=FALSE)
{
    if (verbose)
        message(wmsg("Adding the intdata to the combined germline db"), " ... ",
                appendLF=FALSE)
    intdata1 <- load_intdata(db_name1)
    intdata2 <- load_intdata(db_name2)
    intdata1$allele_name <- add_suffix(intdata1[ , "allele_name"], suffix1)
    intdata2$allele_name <- add_suffix(intdata2[ , "allele_name"], suffix2)
    intdata <- rbind(intdata1, intdata2)
    write_ndm_data_to_db(intdata, destdir, check.and.reorder=TRUE)
    if (verbose)
        message("ok.\n")
}

.combine_db_auxdata <- function(destdir, db_name1, db_name2, suffix1, suffix2,
                                verbose=FALSE)
{
    if (verbose)
        message(wmsg("Adding the auxdata to the combined germline db"), " ... ",
                appendLF=FALSE)
    auxdata1 <- load_auxdata(db_name1)
    auxdata2 <- load_auxdata(db_name2)
    auxdata1$allele_name <- add_suffix(auxdata1[ , "allele_name"], suffix1)
    auxdata2$allele_name <- add_suffix(auxdata2[ , "allele_name"], suffix2)
    auxdata <- rbind(auxdata1, auxdata2)
    write_auxdata_to_db(auxdata, destdir, check.and.reorder=TRUE)
    if (verbose)
        message("ok.\n")
}

.do_combine_germline_dbs <- function(destdir,
                                     db_name1, db_name2, suffix1, suffix2,
                                     verbose=FALSE)
{
    if (!isSingleNonWhiteString(db_name1))
        stop(wmsg("'db_name1' must be a single (non-empty) string"))
    if (!isSingleNonWhiteString(db_name2))
        stop(wmsg("'db_name2' must be a single (non-empty) string"))
    if (!isSingleString(suffix1))
        stop(wmsg("'suffix1' must be a single string"))
    if (has_whitespace(suffix1))
        stop(wmsg("'suffix1' cannot contain whitespace characters"))
    if (!isSingleString(suffix2))
        stop(wmsg("'suffix2' must be a single string"))
    if (has_whitespace(suffix2))
        stop(wmsg("'suffix2' cannot contain whitespace characters"))
    for (region_type in VDJ_REGION_TYPES)
        .combine_region_dbs(destdir, db_name1, db_name2, suffix1, suffix2,
                            region_type, verbose=verbose)
    .combine_db_intdata(destdir, db_name1, db_name2, suffix1, suffix2,
                        verbose=verbose)
    .combine_db_auxdata(destdir, db_name1, db_name2, suffix1, suffix2,
                        verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### combine_germline_dbs()
###

### Combine two existing germline dbs, typically from two different
### organisms, into one hybrid germline db.
### Note that allele names are usually not unique across species.
### Arguments 'suffix1' and 'suffix2' can be used to disambiguate them.
combine_germline_dbs <- function(db_name,
                                 db_name1, db_name2, suffix1, suffix2,
                                 overwrite=FALSE, verbose=FALSE)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    if (has_whitespace(db_name))
        stop(wmsg("'db_name' cannot contain whitespace characters"))
    if (!has_prefix(db_name, "comb"))
        stop(wmsg("'db_name' must start with \"comb\""))

    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    db_path <- get_germline_db_path(db_name)
    if (dir.exists(db_path) && !overwrite)
        stop_on_existing_cached_germline_db(db_name)

    ## We first create the combined germline db in a temporary folder, and,
    ## only if successful, we replace 'db_path' with the temporary folder.
    ## Otherwise we destroy the temporary folder and raise an error. This
    ## achieves atomicity and avoids loosing the content of the existing
    ## 'db_path' in case something goes wrong.
    tmp_destdir <- tempfile("germline_db_")
    dir.create(tmp_destdir)
    on.exit(nuke_file(tmp_destdir))
    .do_combine_germline_dbs(tmp_destdir, db_name1, db_name2, suffix1, suffix2,
                             verbose=verbose)
    rename_file(tmp_destdir, db_path, replace=TRUE)

    invisible(db_name)
}

