### =========================================================================
### list_c_region_dbs() and related
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_builtin_c_region_dbs()
###

### Do NOT call this in .onLoad()! It relies on create_c_region_db()
### which requires that Perl and a valid IgBLAST installation (for
### the 'edit_imgt_file.pl' script) are already available on the machine.
### However, none of these things are guaranteed to be available at load-time,
### especially if it's the first time that the package gets loaded on the user
### machine (e.g. right after installing the package from source).
.create_builtin_c_region_dbs <- function(destdir)
{
    stopifnot(isSingleNonWhiteString(destdir), !dir.exists(destdir))

    ## We first create the dbs in a temporary folder, and, only if successful,
    ## rename the temporary folder to 'destdir'. Otherwise we destroy the
    ## temporary folder and raise an error. This achieves atomicity.
    tmp_destdir <- tempfile("builtin_c_region_dbs_")
    dir.create(tmp_destdir, recursive=TRUE)
    on.exit(nuke_file(tmp_destdir))

    ## Create IMGT C-region dbs.
    IMGT_c_region_dir <- system.file(package="igblastr",
                                     "extdata", "constant_regions", "IMGT",
                                     mustWork=TRUE)
    organism_paths <- list.dirs(IMGT_c_region_dir, recursive=FALSE)
    for (organism_path in organism_paths) {
        db_name <- form_IMGT_c_region_db_name(organism_path)
        db_path <- file.path(tmp_destdir, db_name)
        create_c_region_db(organism_path, db_path)
    }

    ## Any other built-in C-region dbs to create?

    ## Everyting went fine so we can rename 'tmp_destdir' to 'destdir'.
    rename_file(tmp_destdir, destdir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_c_region_dbs_path()
###

### Returns path to C_REGION_DBS cache compartment (see R/cache-utils.R for
### details about igblastr's cache organization).
### When 'init.path=TRUE':
### - if the path to return exists then no further action is performed;
### - if the path to return does NOT exist then it's created and populated
###   with the built-in C-region dbs.
### This means that the returned path is only guaranteed to exist
### when 'init.path' is set to TRUE.
.get_c_region_dbs_path <- function(init.path=FALSE)
{
    stopifnot(isTRUEorFALSE(init.path))
    c_region_dbs_path <- igblastr_cache(C_REGION_DBS)
    if (!dir.exists(c_region_dbs_path) && init.path)
        .create_builtin_c_region_dbs(c_region_dbs_path)
    c_region_dbs_path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_c_region_dbs()
###

### 'long.listing' is ignored when 'names.only' is TRUE.
### Returns a c_region_dbs_df object (data.frame extension) by default.
list_c_region_dbs <- function(builtin.only=FALSE,
                              names.only=FALSE, long.listing=FALSE)
{
    c_region_dbs_path <- .get_c_region_dbs_path(TRUE)  # guaranteed to exist
    ans <- list_dbs(c_region_dbs_path, what="C-region",
                    builtin.only=builtin.only,
                    names.only=names.only, long.listing=long.listing)
    if (is.data.frame(ans))
        class(ans) <- c("c_region_dbs_df", class(ans))
    ans
}

print.c_region_dbs_df <- function(x, ...)
{
    c_region_dbs_path <- .get_c_region_dbs_path(TRUE)  # guaranteed to exist
    print_dbs_df(x, c_region_dbs_path, what="C-region")
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
### make_c_region_db_path()
###

### Not exported!
### Note that the returned path is NOT guaranteed to exist.
make_c_region_db_path <- function(db_name)
{
    if (!isSingleNonWhiteString(db_name))
        stop(wmsg("'db_name' must be a single (non-empty) string"))
    stopifnot(db_name != "USING")
    c_region_dbs_path <- .get_c_region_dbs_path(TRUE)  # guaranteed to exist
    file.path(c_region_dbs_path, db_name)              # NOT guaranteed to exist
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### use_c_region_db()
###

### Returns "" if no db is currently in use.
.get_c_region_db_in_use <- function(verbose=FALSE)
{
    c_region_dbs_path <- .get_c_region_dbs_path(TRUE)  # guaranteed to exist
    db_path <- get_db_in_use(c_region_dbs_path, what="C-region")
    if (db_path == "")
        return(db_path)
    make_blastdbs(db_path, verbose=verbose)
    basename(db_path)
}

### Passing 'db_name=""' will cancel the current selection.
use_c_region_db <- function(db_name=NULL, verbose=FALSE)
{
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))
    if (is.null(db_name))
        return(.get_c_region_db_in_use(verbose=verbose))

    if (!isSingleString(db_name))
        stop(wmsg("'db_name' must be a single string"))

    if (db_name != "") {
        .check_c_region_db_name(db_name)
        db_path <- make_c_region_db_path(db_name)
        make_blastdbs(db_path, verbose=verbose)
    }
    ## Returns 'db_name' invisibly.
    set_db_in_use("C-region", db_name, verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_c_region_db()
###

### Returns the C regions in a DNAStringSet object.
load_c_region_db <- function(db_name)
{
    .check_c_region_db_name(db_name)
    db_path <- make_c_region_db_path(db_name)
    fasta_file <- get_db_fasta_file(db_path, "C")
    readDNAStringSet(fasta_file)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_c_region_blastdbs()
###

### Not used at the moment and not exported!
clean_c_region_blastdbs <- function()
{
    c_region_dbs_path <- .get_c_region_dbs_path()  # NOT guaranteed to exist
    if (dir.exists(c_region_dbs_path)) {
        all_db_names <- list_c_region_dbs(names.only=TRUE)
        for (db_name in all_db_names) {
            db_path <- make_c_region_db_path(db_name)
            clean_blastdbs(db_path)
        }
    }
}

