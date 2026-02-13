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
### .path_to_IMGT_germline_fasta_store()
###

.path_to_IMGT_germline_fasta_store <- function(local_store, organism,
                                               loci_prefix)
{
    stopifnot(isSingleNonWhiteString(loci_prefix))
    organism_path <- find_organism_in_IMGT_local_store(organism, local_store)
    organism <- basename(organism_path)
    fasta_store <- file.path(organism_path, loci_prefix)
    if (!dir.exists(fasta_store))
        stop(wmsg("cannot find ", loci_prefix, " germline sequences ",
                  "for ", organism, " in IMGT release ", basename(local_store)))
    fasta_store
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_effective_loci()
###

### Returns the intersection between 'wanted_loci' and 'found_loci'
### but with lots of sanity checks, bells and whistles.
.get_effective_loci <- function(wanted_loci, found_loci)
{
    stopifnot(is.character(wanted_loci), is.character(found_loci))
    keep_idx <- which(wanted_loci %in% found_loci)
    if (length(keep_idx) == 0L) {
        what <- if (length(wanted_loci) == 1L) "locus" else "loci"
        in1string <- paste0(wanted_loci, collapse=", ")
        stop(wmsg("no FASTA files found for ", what, " ", in1string))
    }
    missing_loci <- wanted_loci[-keep_idx]
    ## Like 'intersect(found_loci, wanted_loci)' but the returned
    ## intersection is guaranteed to be ordered like in 'found_loci'.
    loci <- found_loci[found_loci %in% wanted_loci]
    if (length(missing_loci) != 0L) {
        what1 <- if (length(missing_loci) == 1L) "locus" else "loci"
        in1string1 <- paste0(missing_loci, collapse=", ")
        what2 <- if (length(loci) == 1L) "locus" else "loci"
        in1string2 <- paste0(loci, collapse=", ")
        warning(wmsg("No FASTA files found for ", what1, " ", in1string1, " ",
                     "--> Installing germline db for ", what2, " ", in1string2,
                     "."),
                immediate.=TRUE)
    }
    loci
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_IMGT_germline_db()
###

.form_IMGT_germline_db_name <- function(fasta_store, loci)
{
    stopifnot(isSingleNonWhiteString(fasta_store), dir.exists(fasta_store))
    stop_if_malformed_loci_vector(loci)
    organism_path <- dirname(fasta_store)
    organism <- basename(organism_path)
    refdir <- dirname(organism_path)
    stopifnot(basename(refdir) == VQUEST_REFERENCE_DIRECTORY)
    local_store <- dirname(refdir)
    release <- basename(local_store)
    sprintf("IMGT-%s.%s.%s", release, organism, paste(loci, collapse="+"))
}

install_IMGT_germline_db <- function(release, organism="Homo sapiens",
                                     tcr.db=FALSE, loci="auto",
                                     without.intdata=FALSE,
                                     force=FALSE, verbose=FALSE, ...)
{
    ## Check arguments.
    if (missing(release))
        .stop_on_missing_release()
    release <- .validate_IMGT_release(release)
    organism <- normalize_IMGT_organism(organism)
    loci <- normalize_user_supplied_loci(loci, tcr.db=tcr.db,
                                         stop.if.missing.regions=TRUE)
    loci_prefix <- extract_loci_prefix(loci)
    if (!isTRUEorFALSE(without.intdata))
        stop(wmsg("'without.intdata' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    ## Download IMGT/V-QUEST release to local store if it's not there already.
    local_store <- .path_to_IMGT_local_store(release)
    if (!dir.exists(local_store))
        download_and_unzip_IMGT_release(release, local_store, ...)

    ## Get path to local FASTA store for IMGT germline sequences.
    fasta_store <- .path_to_IMGT_germline_fasta_store(local_store, organism,
                                                      loci_prefix)

    ## Keep loci for which IMGT actually provides FASTA files.
    found_loci <- list_loci_in_germline_fasta_dir(fasta_store, loci_prefix)
    loci <- .get_effective_loci(loci, found_loci)

    ## Compute 'db_name'.
    db_name <- .form_IMGT_germline_db_name(fasta_store, loci)

    ## Create IMGT germline db.
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    db_path <- file.path(germline_dbs_home, db_name)
    create_germline_db(fasta_store, loci, db_path, force=force, verbose=verbose)

    ## Add internal data to germline db (only for IG dbs for now).
    if (!(without.intdata || tcr.db)) {
        if (verbose)
            message(wmsg("Computing and adding the intdata to ", db_name),
                    " ... ", appendLF=FALSE)
        add_V_ndm_data_to_IMGT_germline_db(db_path)
        if (verbose)
            message("ok\n")
    }

    ## Success!
    message("Germline db ", db_name, " successfully installed.")
    message("Call use_germline_db(\"", db_name, "\") to select it")
    message("as the germline db to use with igblastn().")

    invisible(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_IMGT_c_region_db()
###

.form_IMGT_c_region_db_name <- function(organism, loci, version)
{
    stopifnot(isSingleNonWhiteString(organism))
    sprintf("IMGT.%s.%s.%s", organism, paste(loci, collapse="+"), version)
}

install_IMGT_c_region_db <- function(organism, loci, force=FALSE, verbose=FALSE)
{
    organism <- find_organism_shortname(normalize_IMGT_organism(organism))
    loci <- normalize_user_supplied_loci(loci)
    loci_prefix <- extract_loci_prefix(loci)
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    ## Get path to local FASTA store for IMGT C-region sequences.
    fasta_store <- path_to_IMGT_c_region_fasta_store(organism, loci_prefix)
    if (!dir.exists(fasta_store))
        stop(wmsg("we're not aware of any ", loci_prefix, " C-region ",
                  "sequences available at IMGT for ", organism, ", sorry"))

    ## Keep loci for which IMGT actually provides FASTA files.
    found_loci <- list_loci_in_c_region_fasta_dir(fasta_store, loci_prefix)
    loci <- .get_effective_loci(loci, found_loci)

    ## Compute 'db_name'.
    version <- read_version_file(fasta_store)
    db_name <- .form_IMGT_c_region_db_name(organism, loci, version)

    ## Create IMGT C-region db.
    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    db_path <- file.path(c_region_dbs_home, db_name)
    create_c_region_db(fasta_store, loci, db_path, force=force, verbose=verbose)

    ## Success!
    message("C-region db ", db_name, " successfully installed.")
    message("Call use_c_region_db(\"", db_name, "\") to select it as the")
    message("C-region db to use with igblastn().")

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

