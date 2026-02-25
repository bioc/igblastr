### =========================================================================
### reset_germline_dbs()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .install_missing_builtin_AIRR_human_germline_dbs()
###

.form_builtin_AIRR_human_germline_db_name <- function(fasta_store)
{
    stopifnot(isSingleNonWhiteString(fasta_store))
    loci_in1string <- paste(IG_LOCI, collapse="+")
    version <- basename(dirname(fasta_store))
    db_name <- sprintf("_AIRR.human.%s.%s", loci_in1string, version)
    flavor <- basename(fasta_store)
    stopifnot(flavor %in% c("ref", "src"))
    if (flavor == "src")
        db_name <- paste0(db_name, ".src")
    db_name
}

### Install db only if missing.
.install_builtin_AIRR_human_germline_db <-
    function(install_dir, fasta_store, verbose=FALSE)
{
    db_name <- .form_builtin_AIRR_human_germline_db_name(fasta_store)
    install_germline_db(install_dir, db_name, fasta_store, IG_LOCI,
                        if.exists="no-op", verbose=verbose)
    if (basename(fasta_store) == "ref") {
        db_path <- file.path(install_dir, db_name)
        add_V_ndm_data_to_germline_db(db_path, fasta_store,
                                      domain_system="imgt")
    }
}

.install_missing_builtin_AIRR_human_germline_dbs <-
    function(install_dir, human_dir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(human_dir), dir.exists(human_dir),
              isTRUEorFALSE(verbose))
    subdirs <- list.dirs(human_dir, full.names=FALSE, recursive=FALSE)
    subdirs <- setdiff(subdirs, "diffs")
    fasta_stores <- file.path(human_dir, rep(subdirs, each=2L), c("ref", "src"))
    for (fasta_store in fasta_stores)
        .install_builtin_AIRR_human_germline_db(install_dir, fasta_store,
                                                verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .install_missing_builtin_AIRR_mouse_germline_dbs()
###

.form_builtin_AIRR_mouse_germline_db_name <- function(fasta_store)
{
    stopifnot(isSingleNonWhiteString(fasta_store))
    strain <- basename(fasta_store)
    loci_in1string <- paste(IG_LOCI, collapse="+")
    version <- read_version_file(fasta_store)
    sprintf("_AIRR.mouse.%s.%s.%s", strain, loci_in1string, version)
}

### Install db only if missing.
.install_builtin_AIRR_mouse_germline_db <-
    function(install_dir, fasta_store, verbose=FALSE)
{
    db_name <- .form_builtin_AIRR_mouse_germline_db_name(fasta_store)
    install_germline_db(install_dir, db_name, fasta_store, IG_LOCI,
                        if.exists="no-op", verbose=verbose)
}

.install_missing_builtin_AIRR_mouse_germline_dbs <-
    function(install_dir, mouse_dir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(mouse_dir), dir.exists(mouse_dir),
              isTRUEorFALSE(verbose))
    fasta_stores <- list.dirs(mouse_dir, recursive=FALSE)
    for (fasta_store in fasta_stores)
        .install_builtin_AIRR_mouse_germline_db(install_dir, fasta_store,
                                                verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .install_missing_builtin_AIRR_rhesus_monkey_germline_dbs()
###

.form_builtin_AIRR_rhesus_monkey_germline_db_name <- function(fasta_store)
{
    stopifnot(isSingleNonWhiteString(fasta_store))
    loci_in1string <- paste(IG_LOCI, collapse="+")
    version <- basename(fasta_store)
    sprintf("_AIRR.rhesus_monkey.%s.%s", loci_in1string, version)
}

### Install db only if missing.
.install_builtin_AIRR_rhesus_monkey_germline_db <-
    function(install_dir, fasta_store, verbose=FALSE)
{
    db_name <- .form_builtin_AIRR_rhesus_monkey_germline_db_name(fasta_store)
    install_germline_db(install_dir, db_name, fasta_store, IG_LOCI,
                        if.exists="no-op", verbose=verbose)
    db_path <- file.path(install_dir, db_name)
    add_V_ndm_data_to_germline_db(db_path, fasta_store, domain_system="imgt")
}

.install_missing_builtin_AIRR_rhesus_monkey_germline_dbs <-
    function(install_dir, rhesus_monkey_dir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(rhesus_monkey_dir),
              dir.exists(rhesus_monkey_dir), isTRUEorFALSE(verbose))
    subdirs <- list.dirs(rhesus_monkey_dir, full.names=FALSE, recursive=FALSE)
    subdirs <- setdiff(subdirs, "diffs")
    fasta_stores <- file.path(rhesus_monkey_dir, subdirs)
    for (fasta_store in fasta_stores)
        .install_builtin_AIRR_rhesus_monkey_germline_db(install_dir, fasta_store,
                                                        verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_missing_builtin_germline_dbs()
###

### Not exported!
### 'install_dir' must exist.
install_missing_builtin_germline_dbs <- function(install_dir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(install_dir), dir.exists(install_dir),
              isTRUEorFALSE(verbose))

    AIRR_germline_seq_dir <- system.file(package="igblastr",
                                 "extdata", "germline_sets", "AIRR",
                                 mustWork=TRUE)

    human_dir <- file.path(AIRR_germline_seq_dir, "human")
    .install_missing_builtin_AIRR_human_germline_dbs(install_dir, human_dir,
                                                     verbose=verbose)

    mouse_dir <- file.path(AIRR_germline_seq_dir, "mouse")
    .install_missing_builtin_AIRR_mouse_germline_dbs(install_dir, mouse_dir,
                                                     verbose=verbose)

    rhesus_monkey_dir <- file.path(AIRR_germline_seq_dir, "rhesus_monkey")
    .install_missing_builtin_AIRR_rhesus_monkey_germline_dbs(install_dir,
                                                             rhesus_monkey_dir,
                                                             verbose=verbose)

    ## Any other built-in germline dbs to install?
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reset_germline_dbs()
###

### Has its own tests in tests/testthat/test-reset_germline_dbs.R!
### If 'install_dir' exists, it gets destroyed and replaced with a freshly
### populated directory. In other words, its final content doesn't depend
### on whether it already exists or not. Note that the destroy/replace
### operation is atomic.
.install_all_builtin_germline_dbs <- function(install_dir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(install_dir))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    ## We first create the dbs in a temporary folder, and, only if successful,
    ## rename the temporary folder to 'install_dir'. Otherwise we destroy the
    ## temporary folder and raise an error. This achieves atomicity.
    tmp_install_dir <- tempfile("builtin_germline_dbs_")
    dir.create(tmp_install_dir)
    on.exit(nuke_file(tmp_install_dir))

    install_missing_builtin_germline_dbs(tmp_install_dir, verbose=verbose)

    ## Everything went fine so we can rename 'tmp_install_dir' to 'install_dir'.
    rename_file(tmp_install_dir, install_dir, replace=TRUE)
}

### Will nuke any user-installed db!
reset_germline_dbs <- function(verbose=FALSE)
{
    set_db_in_use("germline", "")  # cancel current selection
    germline_dbs_home <- igblastr_cache(GERMLINE_DBS)
    .install_all_builtin_germline_dbs(germline_dbs_home, verbose=verbose)
}

