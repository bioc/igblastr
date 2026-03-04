### =========================================================================
### download_OGRDB_germline_sequences()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normalize_OGRDB_germline_sets()
### infer_loci_from_OGRDB_set_names()
###

### Not exported!
normalize_OGRDB_germline_sets <- function(germline_sets)
{
    if (!is.numeric(germline_sets) || length(germline_sets) == 0L)
        stop(wmsg("'germline_sets' must be a non-empty integer vector"))
    set_names <- names(germline_sets)
    if (is.null(set_names))
        stop(wmsg("'germline_sets' must have names"))
    if (anyNA(set_names))
        stop(wmsg("the names on 'germline_sets' cannot contain NAs"))
    set_names <- trimws2(set_names)
    if (any(nchar(set_names) == 0L))
        stop(wmsg("the names on 'germline_sets' cannot be empty"))
    if (anyDuplicated(set_names))
        stop(wmsg("the names on 'germline_sets' cannot contain duplicates"))
    if (!is.integer(germline_sets))
        germline_sets <- setNames(as.integer(germline_sets), set_names)
    if (anyNA(germline_sets))
        stop(wmsg("'germline_sets' cannot contain NAs"))
    germline_sets
}

### Not exported!
### Tries to infer the locus embedded in the name of each germline set.
infer_loci_from_OGRDB_set_names <- function(set_names)
{
    stopifnot(is.character(set_names))
    pattern1 <- "IG[HKL]"
    bad_ix <- grep(pattern1, set_names, invert=TRUE)
    if (length(bad_ix) != 0L) {
        in1string <- paste0(set_names[bad_ix], collapse=", ")
        stop(wmsg("cannot guess locus for: ", in1string))
    }
    pattern2 <- paste0("^.*(", pattern1, ").*$")
    sub(pattern2, "\\1", set_names)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_OGRDB_germline_sequences()
###

### Produces one FASTA file per "IMGT group".
### Returns the names of the FASTA files that got created.
.split_OGRDB_fasta_file <- function(fasta_file, set_name, set_locus,
                                    destdir=".")
{
    stopifnot(isSingleNonWhiteString(fasta_file),
              isSingleNonWhiteString(set_name),
              isSingleNonWhiteString(set_locus),
              isSingleNonWhiteString(destdir), dir.exists(destdir))
    alleles <- readDNAStringSet(fasta_file)
    allele_names <- names(alleles)

    allele_loci <- substr(allele_names, 1L, 3L)
    bad_idx <- which(allele_loci != set_locus)
    ## Hopefully this will never happen.
    if (length(bad_idx) != 0L)
        stop(wmsg("based on their names, ", length(bad_idx), " alleles ",
                  "(out of ", length(alleles), ") in OGRDB germline set ",
                  set_name, " don't belong to the ", set_locus, " locus ",
                  "(e.g. ", allele_names[[bad_idx[[1L]]]], ")"))

    allele_regions <- substr(allele_names, 4L, 4L)
    bad_idx <- which(!(allele_regions %in% VDJC_REGION_TYPES))
    ## Hopefully this will never happen either.
    if (length(bad_idx) != 0L) {
        in1string <- paste(VDJC_REGION_TYPES, collapse="/")
        stop(wmsg("based on their names, ", length(bad_idx), " alleles ",
                  "(out of ", length(alleles), ") in OGRDB germline set ",
                  set_name, " don't belong to any of the ", in1string, " ",
                  "regions (e.g. ", allele_names[[bad_idx[[1L]]]], ")"))
    }

    ans <- character(0)
    for (region_type in VDJC_REGION_TYPES) {
        selected_alleles <- alleles[allele_regions == region_type]
        if (length(selected_alleles) == 0L)
            next
        group <- paste0(set_locus, region_type)
        filename <- paste0(group, ".fasta")
        destfile <- file.path(destdir, filename)
        if (file.exists(destfile))
            stop(wmsg("At least 2 germline sets in 'germline_sets' contain ",
                      "alleles for the same \"IMGT group\" (", group, "). ",
                      "This is not supported."))
        writeXStringSet(selected_alleles, destfile)
        ans <- c(ans, filename)
    }
    ans
}

### Returns the list of FASTA files that were produced in an invisible
### character vector that carries the names of the corresponding germline
### sets.
download_OGRDB_germline_sequences <- function(organism="Homo sapiens",
                                              germline_sets, gapped=TRUE,
                                              source_set=FALSE,
                                              destdir=".", overwrite=FALSE,
                                              recache=FALSE, ...)
{
    organism <- normalize_OGRDB_organism(organism)
    germline_sets <- normalize_OGRDB_germline_sets(germline_sets)
    set_names <- names(germline_sets)
    set_loci <- infer_loci_from_OGRDB_set_names(set_names)
    if (!isTRUEorFALSE(gapped))
        stop(wmsg("'gapped' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(source_set))
        stop(wmsg("'source_set' must be TRUE or FALSE"))
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir)) {
        if (file.exists(destdir))
            stop(wmsg(destdir, ": not a directory"))
        stop(wmsg(destdir, ": no such directory"))
    }
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))

    format <- if (gapped) "gapped" else "ungapped"
    tmp_destdir <- tempfile("OGRDB_germline_sequences_")
    dir.create(tmp_destdir)
    on.exit(nuke_file(tmp_destdir))

    filenames_list <- lapply(setNames(seq_along(germline_sets), set_names),
        function(i) {
            set_name <- set_names[[i]]
            set_version <- germline_sets[[i]]
            set_locus <- set_loci[[i]]
            ## Download OGRDB germline set to local store if it's not
            ## already there.
            local_file <- download_OGRDB_germline_set_to_OGRDB_store(organism,
                                         set_name, set_version,
                                         format=format, source_set=source_set,
                                         recache=recache, ...)
            .split_OGRDB_fasta_file(local_file, set_name, set_locus,
                                    destdir=tmp_destdir)
        }
    )

    filenames <- unlist(filenames_list, use.names=FALSE)
    names(filenames) <- rep.int(names(filenames_list), lengths(filenames_list))
    stopifnot(!anyDuplicated(filenames))  # should never happen

    ## Copy files from 'tmp_destdir' to 'destdir'.
    fasta_files <- file.path(tmp_destdir, filenames)
    copy_files_to_dir(fasta_files, destdir, overwrite=overwrite)

    invisible(filenames)
}

