### =========================================================================
### list_c_region_dbs() and related
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_c_region_dbs_home()
###

### The first built-in IMGT TR C-region dbs were added in Sep 2025 in
### igblastr 0.99.15 for human and mouse.
.has_builtin_IMGT_TR_c_region_dbs <- function(c_region_dbs_home)
{
    stopifnot(isSingleNonWhiteString(c_region_dbs_home),
              dir.exists(c_region_dbs_home))
    db_path <- list.files(c_region_dbs_home,
                          pattern="_IMGT\\.human\\.TR",
                          full.names=TRUE)
    length(db_path) == 1L && dir.exists(db_path)
}

### Returns whether the built-in C-region dbs need to be (re)created.
### Creating them is needed if folder 'c_region_dbs_home' does not exist.
### Recreating them is needed if folder 'c_region_dbs_home' exists but is
### out-of-sync with the content of igblastr/inst/extdata/constant_regions/.
.need_to_create_builtin_c_region_dbs <- function(c_region_dbs_home)
{
    stopifnot(isSingleNonWhiteString(c_region_dbs_home))
    if (!dir.exists(c_region_dbs_home))
        return(TRUE)
    ## In igblastr <= 0.99.12, the list of built-in C-region dbs is expected
    ## to be:
    ##     _IMGT.human.IGH+IGK+IGL.202412
    ##     _IMGT.mouse.IGH.202412
    ##     _IMGT.rabbit.IGH.202412
    ## In igblastr 0.99.13 (Aug 2025), we added:
    ##     _IMGT.rat.IGH.202508
    ## In igblastr 0.99.15 (Sep 2025), we replaced _IMGT.mouse.IGH.202412
    ## with _IMGT.mouse.IGH.202509 and added:
    ##     _IMGT.rhesus_monkey.IGH.202509
    ## In igblastr 0.99.16 (Sep 2025), we added:
    ##     _IMGT.human.TRA+TRB+TRG+TRD.202509
    ##     _IMGT.mouse.TRA+TRB+TRG+TRD.202509
    ## So we only check for the presence of the built-in IMGT C-region db
    ## for rhesus_monkey to decide whether the built-in C-region dbs need to
    ## be recreated or not.
    !.has_builtin_IMGT_TR_c_region_dbs(c_region_dbs_home)
}

### Returns path to C_REGION_DBS cache compartment (see R/cache-utils.R for
### details about igblastr's cache organization).
### When 'init.path=TRUE':
### - if the path to return exists then no further action is performed;
### - if the path to return does NOT exist then it's created and populated
###   with the built-in C-region dbs.
### This means that the returned path is only guaranteed to exist
### when 'init.path' is set to TRUE.
get_c_region_dbs_home <- function(init.path=FALSE)
{
    stopifnot(isTRUEorFALSE(init.path))
    c_region_dbs_home <- igblastr_cache(C_REGION_DBS)
    if (init.path && .need_to_create_builtin_c_region_dbs(c_region_dbs_home))
        create_all_builtin_c_region_dbs(c_region_dbs_home)
    c_region_dbs_home
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_c_region_dbs()
###

### 'long.listing' is ignored when 'names.only' is TRUE.
### Returns a c_region_dbs_df object (data.frame extension) by default.
list_c_region_dbs <- function(builtin.only=FALSE,
                              names.only=FALSE, long.listing=FALSE)
{
    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    ans <- list_dbs(c_region_dbs_home, what="C-region",
                    builtin.only=builtin.only,
                    names.only=names.only, long.listing=long.listing)
    if (is.data.frame(ans))
        class(ans) <- c("c_region_dbs_df", class(ans))
    ans
}

print.c_region_dbs_df <- function(x, ...)
{
    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    print_dbs_df(x, c_region_dbs_home, what="C-region")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .check_c_region_db_name()
###

.stop_on_invalid_c_region_db_name <- function(db_name)
{
    msg1 <- c("\"", db_name, "\" is not the name of a cached C-region db.")
    msg2 <- c("Use list_c_region_dbs() to list the C-region dbs ",
              "currently installed in the cache (see '?list_c_region_dbs').")
    stop(wmsg(msg1), "\n  ", wmsg(msg2))
}

.check_c_region_db_name <- function(db_name)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    all_db_names <- list_c_region_dbs(names.only=TRUE)
    if (!(db_name %in% all_db_names))
        .stop_on_invalid_c_region_db_name(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### c_region_db_path()
###

### Not exported!
### Note that the returned path is NOT guaranteed to exist.
c_region_db_path <- function(db_name)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    stopifnot(db_name != "USING")
    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    file.path(c_region_dbs_home, db_name)             # NOT guaranteed to exist
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### use_c_region_db()
###

### Returns "" if no db is currently in use.
.get_c_region_db_in_use <- function(verbose=FALSE)
{
    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    db_path <- get_db_in_use(c_region_dbs_home, what="C-region")
    if (db_path == "")
        return(db_path)
    make_blastdbs(db_path, verbose=verbose)
    basename(db_path)
}

.how_to_suppress_use_c_region_db_msg <- function(db_name)
{
    msg1 <- "To suppress this message, use:"
    msg2 <- c("suppressMessages(use_c_region_db(\"", db_name, "\"))")
    c(wmsg(msg1), "\n    ", wmsg(msg2))
}

.note_on_selecting_IMGT_c_region_db <- function(db_name)
{
    is_imgt_db <- grepl("^_?IMGT\\.", db_name)
    if (!is_imgt_db)
        return()
    message("  ", wmsg(IMGT_TERMS_OF_USE), "\n\n  ",
            .how_to_suppress_use_c_region_db_msg(db_name))
}

.select_c_region_db <- function(db_name, verbose=FALSE)
{
    if (!isSingleString(db_name))
        stop(wmsg("'db_name' must be a single string"))

    if (db_name != "") {
        .check_c_region_db_name(db_name)
        .note_on_selecting_IMGT_c_region_db(db_name)
        db_path <- c_region_db_path(db_name)
        make_blastdbs(db_path, verbose=verbose)
    }

    ## Returns 'db_name' invisibly.
    set_db_in_use("C-region", db_name, verbose=verbose)
}

### Passing 'db_name=""' will cancel the current selection.
use_c_region_db <- function(db_name=NULL, verbose=FALSE)
{
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))
    if (is.null(db_name))
        return(.get_c_region_db_in_use(verbose=verbose))
    .select_c_region_db(db_name, verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_c_region_db()
###

### Returns the C regions in a DNAStringSet object.
load_c_region_db <- function(db_name)
{
    .check_c_region_db_name(db_name)
    db_path <- c_region_db_path(db_name)
    db_fasta_file <- get_db_fasta_file(db_path, "C")
    readDNAStringSet(db_fasta_file)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_c_region_blastdbs()
###

### Not used at the moment and not exported!
clean_c_region_blastdbs <- function()
{
    c_region_dbs_home <- get_c_region_dbs_home()  # NOT guaranteed to exist
    if (dir.exists(c_region_dbs_home)) {
        all_db_names <- list_c_region_dbs(names.only=TRUE)
        for (db_name in all_db_names) {
            db_path <- c_region_db_path(db_name)
            clean_blastdbs(db_path)
        }
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rm_c_region_db()
###

rm_c_region_db <- function(db_name)
{
    .check_c_region_db_name(db_name)
    if (has_prefix(db_name, "_"))
        stop(wmsg("cannot remove a built-in C-region db"))

    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    db_in_use_path <- get_db_in_use(c_region_dbs_home, what="C-region")
    if (db_in_use_path != "" && basename(db_in_use_path) == db_name)
        set_db_in_use("C-region", "")  # cancel current selection

    db_path <- c_region_db_path(db_name)
    nuke_file(db_path)
}

