### =========================================================================
### reset_c_region_dbs()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_builtin_IMGT_c_region_db()
###

.form_builtin_IMGT_c_region_db_name <- function(organism, loci, version)
{
    stopifnot(isSingleNonWhiteString(organism))
    sprintf("_IMGT.%s.%s.%s", organism, paste(loci, collapse="+"), version)
}

.create_builtin_IMGT_c_region_db <- function(fasta_store, loci_prefix,
                                             organism, destdir,
                                             verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(verbose))
    loci <- list_loci_in_c_region_fasta_dir(fasta_store, loci_prefix)
    version <- read_version_file(fasta_store)
    db_name <- .form_builtin_IMGT_c_region_db_name(organism, loci, version)
    db_path <- file.path(destdir, db_name)
    create_c_region_db(fasta_store, loci, db_path, verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### path_to_IMGT_c_region_fasta_store()
###

### Not exported!
path_to_IMGT_c_region_fasta_store <- function(organism, loci_prefix)
{
    stopifnot(isSingleNonWhiteString(organism),
              organism %in% names(LATIN_NAMES),
              isSingleNonWhiteString(loci_prefix))
    IMGT_c_region_store <- system.file(package="igblastr",
                                       "extdata", "constant_regions", "IMGT",
                                       mustWork=TRUE)
    fasta_store <- file.path(IMGT_c_region_store, organism, loci_prefix)
    if (loci_prefix == "IG")
        fasta_store <- file.path(fasta_store, "14.1")
    fasta_store
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reset_c_region_dbs()
###

### If 'destdir' exists, it gets destroyed and replaced with a freshly
### populated directory. In other words, its final content doesn't depend
### on whether it already exists or not.
.create_all_builtin_c_region_dbs <- function(destdir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    ## We first create the dbs in a temporary folder, and, only if successful,
    ## rename the temporary folder to 'destdir'. Otherwise we destroy the
    ## temporary folder and raise an error. This achieves atomicity.
    tmp_destdir <- tempfile("builtin_c_region_dbs_")
    dir.create(tmp_destdir)
    on.exit(nuke_file(tmp_destdir))

    ## Create IMGT C-region dbs.
    for (organism in names(LATIN_NAMES)) {
        fasta_store <- path_to_IMGT_c_region_fasta_store(organism, "IG")
        .create_builtin_IMGT_c_region_db(fasta_store, "IG", organism,
                                         tmp_destdir, verbose=verbose)
        fasta_store <- path_to_IMGT_c_region_fasta_store(organism, "TR")
        if (dir.exists(fasta_store))
            .create_builtin_IMGT_c_region_db(fasta_store, "TR", organism,
                                             tmp_destdir, verbose=verbose)
    }

    ## Any other built-in C-region dbs to create?

    ## Everything went fine so we can rename 'tmp_destdir' to 'destdir'.
    rename_file(tmp_destdir, destdir, replace=TRUE)
}

### Will nuke any user-installed db!
reset_c_region_dbs <- function(verbose=FALSE)
{
    set_db_in_use("C-region", "")  # cancel current selection
    c_region_dbs_home <- igblastr_cache(C_REGION_DBS)
    .create_all_builtin_c_region_dbs(c_region_dbs_home, verbose=verbose)
}

