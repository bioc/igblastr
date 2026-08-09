### =========================================================================
### reset_c_region_dbs()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .preinstall_IMGT_c_region_db()
###

.form_preinstalled_IMGT_c_region_db_name <- function(organism, loci, version)
{
    stopifnot(isSingleNonWhiteString(organism))
    sprintf("_IMGT.%s.%s.%s", organism, paste(loci, collapse="+"), version)
}

.preinstall_IMGT_c_region_db <- function(fasta_store, loci_prefix,
                                         organism, destdir,
                                         verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(verbose))
    loci <- list_loci_in_c_region_fasta_dir(fasta_store, loci_prefix)
    version <- read_version_file(fasta_store)
    db_name <- .form_preinstalled_IMGT_c_region_db_name(organism, loci, version)
    db_path <- file.path(destdir, db_name)
    ## Disambiguation needed for "IG" dbs for human, mouse, rabbit, rat, and
    ## rhesus_monkey.
    disambiguate <- (loci_prefix == "IG") &&
        (organism %in% c("human", "mouse", "rabbit", "rat", "rhesus_monkey"))
    create_c_region_db(fasta_store, loci, db_path,
                       disambiguate.allele.names=disambiguate,
                       verbose=verbose)
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

### Has its own unit tests in tests/testthat/test-reset_c_region_dbs.R!
### If 'install_dir' exists, it gets destroyed and replaced with a freshly
### populated directory. In other words, its final content doesn't depend
### on whether it already exists or not.
.preinstall_c_region_dbs <- function(install_dir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(install_dir))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    ## We first create the dbs in a temporary folder, and, only if successful,
    ## rename the temporary folder to 'install_dir'. Otherwise we destroy the
    ## temporary folder and raise an error. This achieves atomicity.
    tmp_install_dir <- tempfile("preinstalled_c_region_dbs_")
    dir.create(tmp_install_dir)
    on.exit(nuke_file(tmp_install_dir))

    ## Create IMGT C-region dbs.
    for (organism in names(LATIN_NAMES)) {
        fasta_store <- path_to_IMGT_c_region_fasta_store(organism, "IG")
        .preinstall_IMGT_c_region_db(fasta_store, "IG", organism,
                                     tmp_install_dir, verbose=verbose)
        fasta_store <- path_to_IMGT_c_region_fasta_store(organism, "TR")
        if (dir.exists(fasta_store))
            .preinstall_IMGT_c_region_db(fasta_store, "TR", organism,
                                         tmp_install_dir, verbose=verbose)
    }

    ## Any other C-region db to preinstall?

    ## Everything went fine so we can rename 'tmp_install_dir' to 'install_dir'.
    rename_file(tmp_install_dir, install_dir, replace=TRUE)
}

### Will nuke any user-installed db!
reset_c_region_dbs <- function(verbose=FALSE)
{
    set_db_in_use("C-region", "")  # cancel current selection
    c_region_dbs_home <- igblastr_cache(C_REGION_DBS)
    .preinstall_c_region_dbs(c_region_dbs_home, verbose=verbose)
}

