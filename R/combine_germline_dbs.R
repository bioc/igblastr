### =========================================================================
### combine_germline_dbs() & combine_c_region_dbs()
### -------------------------------------------------------------------------
###


.check_comb_db_name <- function(db_name)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    if (has_whitespace(db_name))
        stop(wmsg("'db_name' cannot contain whitespace characters"))
    if (!has_prefix(db_name, "comb"))
        stop(wmsg("'db_name' must start with \"comb\""))
}

.check_suffixes <- function(suffix1, suffix2)
{
    if (!isSingleString(suffix1))
        stop(wmsg("'suffix1' must be a single string"))
    if (has_whitespace(suffix1))
        stop(wmsg("'suffix1' cannot contain whitespace characters"))
    if (!isSingleString(suffix2))
        stop(wmsg("'suffix2' must be a single string"))
    if (has_whitespace(suffix2))
        stop(wmsg("'suffix2' cannot contain whitespace characters"))
    if (suffix1 == suffix2)
        stop(wmsg("'suffix1' and 'suffix2' must be different"))
}

.stop_on_colliding_allele_names <- function()
{
    msg1 <- "Allele name collisions!"
    msg2 <- c("Even after adding 'suffix1' to the allele names of the first ",
              "db, and 'suffix2' to the allele names of the second db, some ",
	      "allele names are still colliding. This can be avoided by using ",
              "more specific suffixes.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2))
}

.combine_locus_alleles <- function(locus, alleles1, alleles2, suffix1, suffix2)
{
    stopifnot(isSingleNonWhiteString(locus),
              is(alleles1, "DNAStringSet"), is(alleles2, "DNAStringSet"),
              isSingleString(suffix1), isSingleString(suffix2))
    locus_alleles1 <- alleles1[substr(names(alleles1), 1L, 3L) == locus]
    locus_alleles2 <- alleles2[substr(names(alleles2), 1L, 3L) == locus]
    ## Add user-supplied suffixes to disambiguate allele names.
    names(locus_alleles1) <- add_suffix(names(locus_alleles1), suffix1)
    names(locus_alleles2) <- add_suffix(names(locus_alleles2), suffix2)
    locus_alleles <- c(locus_alleles1, locus_alleles2)
    if (anyDuplicated(names(locus_alleles)))
        .stop_on_colliding_allele_names()
    locus_alleles
}

### We only support combining region dbs for the IG loci at the moment.
.combine_region_dbs <- function(destdir, db_path1, db_path2, suffix1, suffix2,
                                region_type, verbose=FALSE)
{
    alleles1 <- readDNAStringSet(get_db_fasta_file(db_path1, region_type))
    alleles2 <- readDNAStringSet(get_db_fasta_file(db_path2, region_type))
    loci1 <- substr(names(alleles1), 1, 3)
    loci2 <- substr(names(alleles2), 1, 3)
    stopifnot(all(loci1 %in% IG_LOCI), all(loci2 %in% IG_LOCI))

    locus2alleles <- lapply(
        setNames(IG_LOCI, IG_LOCI),
        .combine_locus_alleles, alleles1, alleles2, suffix1, suffix2
    )

    locus2alleles <- locus2alleles[lengths(locus2alleles) != 0L]
    stopifnot(length(locus2alleles) != 0L)

    loci <- names(locus2alleles)
    fasta_files <- vapply(loci,
        function(locus) {
            alleles <- locus2alleles[[locus]]
            fasta_file <- tempfile(fileext=".fasta")
            writeXStringSet(alleles, fasta_file)
            if (verbose) {
                msg <- c("Combined and saved ", length(alleles), " ",
                         region_type, " alleles for locus ", locus, " to ",
                         "temp file ", basename(fasta_file))
                message(wmsg(msg))
            }
            fasta_file
        },
        character(1))
    if (verbose)
        message("")
    on.exit(unlink(fasta_files))

    names(fasta_files) <- paste0(loci, region_type, ".fasta")
    create_region_db(fasta_files, destdir, region_type=region_type,
                     verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### combine_germline_dbs()
###

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
                                     db_path1, db_path2, suffix1, suffix2,
                                     verbose=FALSE)
{
    for (region_type in VDJ_REGION_TYPES)
        .combine_region_dbs(destdir, db_path1, db_path2, suffix1, suffix2,
                            region_type, verbose=verbose)
    db_name1 <- basename(db_path1)
    db_name2 <- basename(db_path2)
    .combine_db_intdata(destdir, db_name1, db_name2, suffix1, suffix2,
                        verbose=verbose)
    .combine_db_auxdata(destdir, db_name1, db_name2, suffix1, suffix2,
                        verbose=verbose)
}

### Combine two existing germline dbs, typically (but not necessarily) from
### two different organisms, into one hybrid germline db.
### Note that allele names are usually not unique across species.
### Arguments 'suffix1' and 'suffix2' can be used to disambiguate them.
combine_germline_dbs <- function(db_name,
                                 db_name1, db_name2, suffix1, suffix2,
                                 overwrite=FALSE, verbose=FALSE)
{
    .check_comb_db_name(db_name)
    check_germline_db_name(db_name1, what="'db_name1'")
    check_germline_db_name(db_name2, what="'db_name2'")
    .check_suffixes(suffix1, suffix2)

    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    db_path <- get_germline_db_path(db_name)
    if (dir.exists(db_path) && !overwrite)
        stop_on_existing_germline_db(db_name)

    db_path1 <- get_germline_db_path(db_name1)
    db_path2 <- get_germline_db_path(db_name2)

    ## We first create the combined germline db in a temporary folder, and,
    ## only if successful, we replace 'db_path' with the temporary folder.
    ## Otherwise we destroy the temporary folder and raise an error. This
    ## achieves atomicity and avoids loosing the content of the existing
    ## 'db_path' in case something goes wrong.
    tmp_destdir <- tempfile("germline_db_")
    dir.create(tmp_destdir)
    on.exit(nuke_file(tmp_destdir))
    .do_combine_germline_dbs(tmp_destdir, db_path1, db_path2, suffix1, suffix2,
                             verbose=verbose)
    rename_file(tmp_destdir, db_path, replace=TRUE)

    invisible(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### combine_c_region_dbs()
###

### Combine two existing C-region dbs, typically (but not necessarily) from
### two different organisms, into one hybrid C-region db.
### Note that allele names are usually not unique across species.
### Arguments 'suffix1' and 'suffix2' can be used to disambiguate them.
combine_c_region_dbs <- function(db_name,
                                 db_name1, db_name2, suffix1, suffix2,
                                 overwrite=FALSE, verbose=FALSE)
{
    .check_comb_db_name(db_name)
    check_c_region_db_name(db_name1, what="'db_name1'")
    check_c_region_db_name(db_name2, what="'db_name2'")
    .check_suffixes(suffix1, suffix2)

    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    db_path <- get_c_region_db_path(db_name)
    if (dir.exists(db_path) && !overwrite)
        stop_on_existing_c_region_db(db_name)

    db_path1 <- get_c_region_db_path(db_name1)
    db_path2 <- get_c_region_db_path(db_name2)

    ## We first create the combined C-region db in a temporary folder, and,
    ## only if successful, we replace 'db_path' with the temporary folder.
    ## Otherwise we destroy the temporary folder and raise an error. This
    ## achieves atomicity and avoids loosing the content of the existing
    ## 'db_path' in case something goes wrong.
    tmp_destdir <- tempfile("c_region_db_")
    dir.create(tmp_destdir)
    on.exit(nuke_file(tmp_destdir))
    .combine_region_dbs(tmp_destdir, db_path1, db_path2, suffix1, suffix2,
                        "C", verbose=verbose)
    rename_file(tmp_destdir, db_path, replace=TRUE)

    invisible(db_name)
}

