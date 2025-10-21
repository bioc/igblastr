### =========================================================================
### list_germline_dbs() and related
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_germline_dbs_home()
###

### Not exported!
### Returns path to GERMLINE_DBS cache compartment (see R/cache-utils.R for
### details about igblastr's cache organization).
### When 'init.path=TRUE':
### - if the path to return exists then no further action is performed;
### - if the path to return does NOT exist then it's created and populated
###   with the built-in germline dbs.
### This means that the returned path is only guaranteed to exist
### when 'init.path' is set to TRUE.
get_germline_dbs_home <- function(init.path=FALSE)
{
    stopifnot(isTRUEorFALSE(init.path))
    germline_dbs_home <- igblastr_cache(GERMLINE_DBS)
    if (init.path) {
        if (dir.exists(germline_dbs_home)) {
            create_missing_builtin_germline_dbs(germline_dbs_home)
        } else {
            create_all_builtin_germline_dbs(germline_dbs_home)
        }
    }
    germline_dbs_home
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_germline_dbs()
###

.OLD_BUILTIN_AIRR_HUMAN_DB <- "_AIRR.human.IGH+IGK+IGL.202501"

.warn_if_old_builtin_AIRR_human_db_exists <- function()
{
    germline_dbs_home <- get_germline_dbs_home()
    old_db <- .OLD_BUILTIN_AIRR_HUMAN_DB
    db_path <- file.path(germline_dbs_home, old_db)
    if (dir.exists(db_path)) {
        new_dbs <- paste0("_AIRR.human.IGH+IGK+IGL.",
                          c("202309", "202309.src", "202410", "202410.src"))
        new_db <- new_dbs[[4L]]
        msg1 <- c("In igblastr 0.99.23, the following built-in germline dbs ",
                  "were added: ", paste(new_dbs, collapse=", "), ".")
        msg2 <- c("Note that ", new_db, " is exactly the same as ",
                  old_db, ", only the name of the db is different.")
        msg3 <- c("The new name is the result of a revisited naming ",
                  "scheme for the built-in AIRR germline dbs for human. ",
                  "See the Value section in '?germline_dbs_home' for ",
                  "more information.")
        msg4 <- c("From now on, please make sure to always use ",
                  "\"", new_db, "\" instead of \"", old_db, "\" in your code.")
        msg5 <- c("To get rid of this warning, remove germline db ",
                  old_db, " with 'rm_germline_db(\"", old_db, "\")'")
        warning(wmsg(msg1), "\n\n  ", wmsg(msg2), "\n  ",
                wmsg(msg3), "\n\n  ", wmsg(msg4), "\n\n  ", wmsg(msg5))
    }
}

### 'long.listing' is ignored when 'names.only' is TRUE.
### Returns a germline_dbs_df object (data.frame extension) by default.
list_germline_dbs <- function(builtin.only=FALSE,
                              names.only=FALSE, long.listing=FALSE)
{
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    ans <- list_dbs(germline_dbs_home, what="germline",
                    builtin.only=builtin.only,
                    names.only=names.only, long.listing=long.listing)
    if (is.data.frame(ans))
        class(ans) <- c("germline_dbs_df", class(ans))
    if (!names.only)
        .warn_if_old_builtin_AIRR_human_db_exists()
    ans
}

print.germline_dbs_df <- function(x, ...)
{
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    print_dbs_df(x, germline_dbs_home, what="germline")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### check_germline_db_name()
###

.stop_on_no_installed_germline_db_yet <- function()
{
    msg <- c("You don't have any installed germline database yet. ",
             "Use any of the install_*_germline_db() function (e.g. ",
             "install_IMGT_germline_db()) to install at least one.")
    stop(wmsg(msg))
}

.stop_on_invalid_germline_db_name <- function(db_name)
{
    msg1 <- c("\"", db_name, "\" is not the name of a cached germline db.")
    msg2 <- c("Use list_germline_dbs() to list the germline dbs ",
              "currently installed in the cache (see '?list_germline_dbs').")
    msg3 <- c("Note that you can use any of the install_*_germline_db() ",
              "function (e.g. install_IMGT_germline_db()) to install ",
              "additional germline dbs in the cache.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

### Not exported!
check_germline_db_name <- function(db_name)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    all_db_names <- list_germline_dbs(names.only=TRUE)
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    if (!(db_name %in% all_db_names))
        .stop_on_invalid_germline_db_name(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_germline_db_path()
###

### Not exported!
### Note that the returned path is NOT guaranteed to exist.
make_germline_db_path <- function(db_name)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    stopifnot(db_name != "USING")
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    file.path(germline_dbs_home, db_name)             # NOT guaranteed to exist
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### use_germline_db()
###

.stop_on_no_selected_germline_db_yet <- function()
{
    msg <- c("You haven't selected the germline database to use ",
             "with igblastn() yet. Please select one with ",
             "use_germline_db(\"<db_name>\"). ",
             "See '?use_germline_db' for more information.")
    stop(wmsg(msg))
}

.get_germline_db_in_use <- function(verbose=FALSE)
{
    all_db_names <- list_germline_dbs(names.only=TRUE)
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    db_path <- get_db_in_use(germline_dbs_home, what="germline")
    if (db_path == "")
        .stop_on_no_selected_germline_db_yet()
    make_blastdbs(db_path, verbose=verbose)
    basename(db_path)
}

use_germline_db <- function(db_name=NULL, verbose=FALSE)
{
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))
    if (is.null(db_name))
        return(.get_germline_db_in_use(verbose=verbose))

    check_germline_db_name(db_name)
    if (db_name == .OLD_BUILTIN_AIRR_HUMAN_DB)
        .warn_if_old_builtin_AIRR_human_db_exists()

    db_path <- make_germline_db_path(db_name)
    make_blastdbs(db_path, verbose=verbose)

    ## Returns 'db_name' invisibly.
    set_db_in_use("germline", db_name, verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_germline_db()
###

.normarg_region_types <- function(region_types=NULL)
{
    if (is.null(region_types))
        return(VDJ_REGION_TYPES)
    if (!is.character(region_types) || anyNA(region_types))
        stop(wmsg("'region_types' must be NULL or ",
                  "a character vector with no NAs"))
    region_types <- toupper(region_types)
    if (length(region_types) == 1L) {
        region_types <- safeExplode(region_types)
    } else if (any(nchar(region_types) != 1L)) {
        stop(wmsg("'region_types' must have single-letter elements"))
    }
    if (!all(region_types %in% VDJ_REGION_TYPES))
        stop(wmsg("'region_types' can only contain letters V, D, or J"))
    region_types
}

### Returns the V, D, and/or J regions in a DNAStringSet object.
load_germline_db <- function(db_name, region_types=NULL)
{
    check_germline_db_name(db_name)
    db_path <- make_germline_db_path(db_name)
    region_types <- .normarg_region_types(region_types)
    fasta_files <- vapply(region_types,
        function(region_type) get_db_fasta_file(db_path, region_type),
        character(1), USE.NAMES=FALSE)
    readDNAStringSet(fasta_files)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_germline_blastdbs()
###

### Not used at the moment and not exported!
clean_germline_blastdbs <- function()
{
    germline_dbs_home <- get_germline_dbs_home()  # NOT guaranteed to exist
    if (dir.exists(germline_dbs_home)) {
        all_db_names <- list_germline_dbs(names.only=TRUE)
        for (db_name in all_db_names) {
            db_path <- make_germline_db_path(db_name)
            clean_blastdbs(db_path)
        }
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rm_germline_db()
###

rm_germline_db <- function(db_name)
{
    check_germline_db_name(db_name)
    if (has_prefix(db_name, "_") && db_name != .OLD_BUILTIN_AIRR_HUMAN_DB)
        stop(wmsg("cannot remove a built-in germline db"))

    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    db_in_use_path <- get_db_in_use(germline_dbs_home, what="germline")
    if (db_in_use_path != "" && basename(db_in_use_path) == db_name)
        set_db_in_use("germline", "")  # cancel current selection

    db_path <- make_germline_db_path(db_name)
    nuke_file(db_path)
}

