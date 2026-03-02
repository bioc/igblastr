### =========================================================================
### download_IMGT_germline_sequences()
### -------------------------------------------------------------------------


.get_path_to_IMGT_store <- function(release=NULL)
{
    IMGT_store <- igblastr_cache(IMGT_STORE)
    if (!is.null(release)) {
        stopifnot(isSingleNonWhiteString(release))
        IMGT_store <- file.path(IMGT_store, release)
    }
    IMGT_store
}

.stop_on_missing_release <- function()
{
    all_releases <- list_IMGT_releases()
    stop(wmsg("Argument 'release' is required and must be set ",
              "to a valid IMGT/V-QUEST release."),
         "\n  ",
         wmsg("Latest IMGT/V-QUEST release is \"", all_releases[[1L]],
              "\" (recommended). Use list_IMGT_releases() to list ",
              "all releases."))
}

.validate_IMGT_release <- function(release)
{
    if (!isSingleNonWhiteString(release))
        stop(wmsg("'release' must be a single (non-empty) string"))
    ## First we try offline validation by checking the IMGT local store.
    if (dir.exists(.get_path_to_IMGT_store(release)))
        return(release)
    ## Off-line validation above failed so we try online validation.
    all_releases <- list_IMGT_releases()
    if (!(release %in% all_releases)) {
        stop(wmsg("\"", release, "\" is not a valid IMGT/V-QUEST release."),
             "\n  ",
             wmsg("Latest IMGT/V-QUEST release is \"", all_releases[[1L]],
                  "\" (recommended). Use list_IMGT_releases() to list ",
                  "all releases."))
    }
    release
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### IMGT_is_up()
###

IMGT_is_up <- function()
    websiteIsUp(IMGT_URL, connecttimeout=get_IMGT_connecttimeout())


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_IMGT_releases()
###

### Returns the IMGT/V-QUEST releases from newest to oldest (latest first).
list_IMGT_releases <- function(recache=FALSE)
{
    latest_release <- get_latest_IMGT_release(recache=recache)
    all_zips <- list_archived_IMGT_zips(recache=recache)
    archived_releases <- sub("^[^0-9]*([-0-9]+).*$", "\\1", all_zips)
    c(latest_release, sort(archived_releases, decreasing=TRUE))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_IMGT_organisms()
###

list_IMGT_organisms <- function(release)
{
    if (missing(release))
        .stop_on_missing_release()
    release <- .validate_IMGT_release(release)

    ## Download IMGT/V-QUEST release to local store if it's not there already.
    IMGT_store <- .get_path_to_IMGT_store(release)
    if (!dir.exists(IMGT_store))
        download_and_unzip_IMGT_release(release, IMGT_store)
    list_organisms_in_IMGT_store(IMGT_store)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_IMGT_release_to_IMGT_store()
###

### Not exported!
### Returns path to the local store.
download_IMGT_release_to_IMGT_store <- function(release, recache=FALSE, ...)
{
    ## Check arguments.
    if (missing(release))
        .stop_on_missing_release()
    release <- .validate_IMGT_release(release)
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))

    ## Download IMGT/V-QUEST release to local store if it's not already there.
    IMGT_store <- .get_path_to_IMGT_store(release)
    if (!dir.exists(IMGT_store) || recache)
        download_and_unzip_IMGT_release(release, IMGT_store, ...)

    IMGT_store
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### find_organism_fasta_store_in_IMGT_store()
###

### Not exported!
find_organism_fasta_store_in_IMGT_store <- function(IMGT_store, organism,
                                                    loci_prefix)
{
    stopifnot(isSingleNonWhiteString(loci_prefix))
    organism_path <- find_organism_in_IMGT_store(organism, IMGT_store)
    fasta_store <- file.path(organism_path, loci_prefix)
    if (!dir.exists(fasta_store)) {
        organism <- basename(organism_path)
        stop(wmsg("cannot find ", loci_prefix, " germline sequences ",
                  "for ", organism, " in IMGT release ", basename(IMGT_store)))
    }
    fasta_store
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_IMGT_germline_sequences()
###

### The downloaded files are cached in the "IMGT local store".
### Returns the names of the downloaded FASTA files in an invisible character
### vector.
download_IMGT_germline_sequences <- function(release, organism="Homo sapiens",
                                             tcr.db=FALSE,
                                             destdir=".", overwrite=FALSE,
                                             recache=FALSE, ...)
{
    organism <- normalize_IMGT_organism(organism)
    if (!isTRUEorFALSE(tcr.db))
        stop(wmsg("'tcr.db' must be TRUE or FALSE"))
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir)) {
        if (file.exists(destdir))
            stop(wmsg(destdir, ": not a directory"))
        stop(wmsg(destdir, ": no such directory"))
    }
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))

    ## Download IMGT/V-QUEST release to local store if it's not already there.
    IMGT_store <- download_IMGT_release_to_IMGT_store(release,
                                                      recache=recache, ...)

    ## Get path to local FASTA store for IMGT germline sequences.
    loci_prefix <- if (tcr.db) "TR" else "IG"
    fasta_store <- find_organism_fasta_store_in_IMGT_store(IMGT_store,
                                                           organism,
                                                           loci_prefix)

    ## Copy files from local FASTA store to 'destdir'.
    fasta_files <- list_fasta_files(fasta_store)
    copy_files_to_dir(fasta_files, destdir, overwrite=overwrite)

    invisible(basename(fasta_files))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### validate_redit_imgt_file_on_IMGT_release()
###
### On Sep 9, 2025, running validate_redit_imgt_file_on_IMGT_release() on
### releases 202343-3, 202405-2, 202518-3 and 202530-1, produced the
### following results:
### - 202330-1:  0 failures
### - 202343-3:  5 failures (5/5 TR files for Mus_musculus_C57BL6J)
### - 202405-2: 16 failures (6/7 IG + 10/10 TR files for Mus_musculus_C57BL6J)
### - 202518-3: 16 failures (6/7 IG + 10/10 TR files for Mus_musculus_C57BL6J)
### - 202530-1:  0 failures
### The failures on various Mus_musculus_C57BL6J files are expected and due
### to Perl script edit_imgt_file.pl not working properly on these files.
### See R/edit_imgt_file.R for more information.
###

### Used in unit tests. Requires Perl!
### Returns number of failures.
validate_redit_imgt_file_on_IMGT_release <- function(release, ...)
{
    ## Download IMGT/V-QUEST release to local store if it's not already there.
    IMGT_store <- download_IMGT_release_to_IMGT_store(release, ...)
    validate_redit_imgt_file(IMGT_store, recursive=TRUE)
}

