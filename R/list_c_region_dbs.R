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
    stopifnot(isSingleNonWhiteString(destdir))

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
    for (organism in names(LATIN_NAMES)) {
        fastadir <- file.path(IMGT_c_region_dir, organism, "14.1")
        db_name <- form_builtin_IMGT_c_region_db_name(fastadir)
        db_path <- file.path(tmp_destdir, db_name)
        create_c_region_db(fastadir, db_path)
    }

    ## Any other built-in C-region dbs to create?

    ## Everyting went fine so we can rename 'tmp_destdir' to 'destdir'.
    rename_file(tmp_destdir, destdir, replace=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_c_region_dbs_path()
###

### The built-in IMGT C-region db for rhesus monkey was added in igblastr
### 0.99.15 (Sep 2025).
.has_builtin_IMGT_c_region_db_for_rhesus_monkey <- function(c_region_dbs_path)
{
    stopifnot(isSingleNonWhiteString(c_region_dbs_path),
              dir.exists(c_region_dbs_path))
    db_path <- list.files(c_region_dbs_path, pattern="_IMGT\\.rhesus_monkey\\.",
                          full.names=TRUE)
    length(db_path) == 1L && dir.exists(db_path)
}

### Returns whether the built-in C-region dbs need to be (re)created.
### Creating them is needed if folder 'c_region_dbs_path' does not exist.
### Recreating them is needed if folder 'c_region_dbs_path' exists but is
### out-of-sync with the content of igblastr/inst/extdata/constant_regions/.
.need_to_create_builtin_c_region_dbs <- function(c_region_dbs_path)
{
    stopifnot(isSingleNonWhiteString(c_region_dbs_path))
    if (!dir.exists(c_region_dbs_path))
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
    ## So we only check for the presence of the built-in IMGT C-region db
    ## for rhesus_monkey to decide whether the built-in C-region dbs need to
    ## be recreated or not.
    !.has_builtin_IMGT_c_region_db_for_rhesus_monkey(c_region_dbs_path)
}

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
    if (.need_to_create_builtin_c_region_dbs(c_region_dbs_path) && init.path)
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

