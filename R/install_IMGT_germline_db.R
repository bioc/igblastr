### =========================================================================
### install_IMGT_germline_db()
### -------------------------------------------------------------------------


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

.path_to_IMGT_local_store <- function(release=NULL)
{
    local_store <- igblastr_cache(IMGT_LOCAL_STORE)
    if (!is.null(release)) {
        stopifnot(isSingleNonWhiteString(release))
        local_store <- file.path(local_store, release)
    }
    local_store
}

.validate_IMGT_release <- function(release)
{
    if (!isSingleNonWhiteString(release))
        stop(wmsg("'release' must be a single (non-empty) string"))
    ## First we try offline validation by checking the IMGT local store.
    if (dir.exists(.path_to_IMGT_local_store(release)))
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

list_IMGT_organisms <- function(release)
{
    if (missing(release))
        .stop_on_missing_release()
    release <- .validate_IMGT_release(release)

    ## Download IMGT/V-QUEST release to local store if it's not there already.
    local_store <- .path_to_IMGT_local_store(release)
    if (!dir.exists(local_store))
        download_and_unzip_IMGT_release(release, local_store)
    list_organisms_in_IMGT_local_store(local_store)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .form_IMGT_germline_db_name()
###

.form_IMGT_germline_db_name <- function(fasta_store)
{
    stopifnot(isSingleNonWhiteString(fasta_store), dir.exists(fasta_store),
              basename(fasta_store) %in% c("IG", "TR"))
    organism_path <- dirname(fasta_store)
    organism <- basename(organism_path)
    refdir <- dirname(organism_path)
    stopifnot(basename(refdir) == VQUEST_REFERENCE_DIRECTORY)
    local_store <- dirname(refdir)
    release <- basename(local_store)
    tcr.db <- basename(fasta_store) == "TR"
    loci <- get_loci_from_input_germline_fasta_set(fasta_store, tcr.db=tcr.db)
    sprintf("IMGT-%s.%s.%s", release, organism, paste(loci, collapse="+"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_IMGT_germline_db()
###

install_IMGT_germline_db <- function(release, organism="Homo sapiens",
                                     tcr.db=FALSE, force=FALSE, ...)
{
    ## Check arguments.
    if (missing(release))
        .stop_on_missing_release()
    release <- .validate_IMGT_release(release)
    organism <- normalize_IMGT_organism(organism)
    if (!isTRUEorFALSE(tcr.db))
        stop(wmsg("'tcr.db' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))

    ## Download IMGT/V-QUEST release to local store if it's not there already.
    local_store <- .path_to_IMGT_local_store(release)
    if (!dir.exists(local_store))
        download_and_unzip_IMGT_release(release, local_store, ...)

    ## Compute 'fasta_store'.
    organism_path <- find_organism_in_IMGT_local_store(organism, local_store)
    organism <- basename(organism_path)
    fasta_store <- file.path(organism_path, if (tcr.db) "TR" else "IG")
    if (!dir.exists(fasta_store))
        stop(wmsg("cannot find ", basename(fasta_store), " germline ",
                  "sequences for ", organism, " in IMGT release ", release))

    ## Compute 'db_name'.
    db_name <- .form_IMGT_germline_db_name(fasta_store)

    ## Create IMGT germline db.
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    db_path <- file.path(germline_dbs_home, db_name)
    create_germline_db(fasta_store, db_path, tcr.db=tcr.db, force=force)

    ## Success!
    message("Germline db ", db_name, " successfully installed.")
    message("Call use_germline_db(\"", db_name, "\") to select it")
    message("as the germline database to use with igblastn().")

    invisible(db_name)
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
    release <- .validate_IMGT_release(release)
    ## Download IMGT/V-QUEST release to local store if it's not there already.
    local_store <- .path_to_IMGT_local_store(release)
    if (!dir.exists(local_store))
        download_and_unzip_IMGT_release(release, local_store, ...)
    validate_redit_imgt_file(local_store, recursive=TRUE)
}

