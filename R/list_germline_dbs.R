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
            reset_germline_dbs()
        }
    }
    germline_dbs_home
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Handle old built-in human AIRR db graciously
###

OLD_BUILTIN_AIRR_HUMAN_DB <- "_AIRR.human.IGH+IGK+IGL.202501"

.warn_about_old_builtin_AIRR_human_db <- function()
{
    old_db <- OLD_BUILTIN_AIRR_HUMAN_DB
    new_dbs <- paste0("_AIRR.human.IGH+IGK+IGL.",
                      c("202309", "202309.src", "202410", "202410.src"))
    new_db <- new_dbs[[4L]]
    msg1 <- c("In igblastr 0.99.23, the following built-in germline dbs ",
              "were added: ", paste(new_dbs, collapse=", "), ".")
    msg2 <- c("Note that ", new_db, " is exactly the same as ",
              old_db, ", only the name of the db is different.")
    msg3 <- c("The new name is the result of a revisited naming ",
              "scheme for the built-in AIRR germline dbs for human. ",
              "See the Value section in '?list_germline_dbs' for ",
              "more information.")
    msg4 <- c("From now on, please make sure to always use ",
              "\"", new_db, "\" instead of \"", old_db, "\" in your code.")
    msg5 <- c("To get rid of this warning, remove germline db ",
              old_db, " with 'rm_germline_db(\"", old_db, "\")'")
    warning(wmsg(msg1), "\n\n  ", wmsg(msg2), "\n  ",
            wmsg(msg3), "\n\n  ", wmsg(msg4), "\n\n  ", wmsg(msg5))
}

warn_if_old_builtin_AIRR_human_db_exists <- function()
{
    germline_dbs_home <- get_germline_dbs_home()
    old_db <- OLD_BUILTIN_AIRR_HUMAN_DB
    db_path <- file.path(germline_dbs_home, old_db)
    if (dir.exists(db_path))
        .warn_about_old_builtin_AIRR_human_db()
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_germline_dbs()
###

### 'long.listing' is ignored when 'names.only' is TRUE.
### Returns a germline_dbs_df object (data.frame extension) by default.
list_germline_dbs <- function(builtin.only=FALSE, with.intdata.only=FALSE,
                              names.only=FALSE, long.listing=FALSE)
{
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    ans <- list_dbs(germline_dbs_home, what="germline",
                    builtin.only=builtin.only,
                    with.intdata.only=with.intdata.only,
                    names.only=names.only, long.listing=long.listing)
    if (is.data.frame(ans))
        class(ans) <- c("germline_dbs_df", class(ans))
    if (!names.only)
        warn_if_old_builtin_AIRR_human_db_exists()
    ans
}

print.germline_dbs_df <- function(x, ...)
{
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    print_dbs_df(x, germline_dbs_home, what="germline")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_germline_db_names()
###

.stop_on_no_installed_germline_db_yet <- function()
{
    msg <- c("You don't have any installed germline database yet. ",
             "Use any of the install_*_germline_db() function (e.g. ",
             "install_IMGT_germline_db()) to install at least one.")
    stop(wmsg(msg))
}

### Not exported!
list_germline_db_names <- function()
{
    all_db_names <- list_germline_dbs(names.only=TRUE)
    ## Should never happen because of the built-in AIRR dbs but we keep this
    ## check anyways just in case we get rid of the built-in AIRR dbs in the
    ## future.
    if (length(all_db_names) == 0L)
        .stop_on_no_installed_germline_db_yet()
    all_db_names
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### check_germline_db_name()
###

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
    all_db_names <- list_germline_db_names()
    if (!(db_name %in% all_db_names))
        .stop_on_invalid_germline_db_name(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_germline_db_path()
###

### Not exported!
### Note that the returned path is NOT guaranteed to exist.
get_germline_db_path <- function(db_name)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    stopifnot(db_name != "USING")
    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    file.path(germline_dbs_home, db_name)             # NOT guaranteed to exist
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rm_germline_db()
###

rm_germline_db <- function(db_name)
{
    check_germline_db_name(db_name)
    if (has_prefix(db_name, "_") && db_name != OLD_BUILTIN_AIRR_HUMAN_DB)
        stop(wmsg("cannot remove a built-in germline db"))

    germline_dbs_home <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    db_in_use_path <- get_db_in_use(germline_dbs_home, what="germline")
    if (db_in_use_path != "" && basename(db_in_use_path) == db_name)
        set_db_in_use("germline", "")  # cancel current selection

    db_path <- get_germline_db_path(db_name)
    nuke_file(db_path)
}

