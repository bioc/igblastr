### =========================================================================
### reset_germline_dbs()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_missing_builtin_AIRR_human_germline_dbs()
###

.form_builtin_AIRR_human_germline_db_name <- function(fasta_store)
{
    loci_in1string <- paste(IG_LOCI, collapse="+")
    version <- basename(dirname(fasta_store))
    db_name <- sprintf("_AIRR.human.%s.%s", loci_in1string, version)
    flavor <- basename(fasta_store)
    stopifnot(flavor %in% c("ref", "src"))
    if (flavor == "src")
        db_name <- paste0(db_name, ".src")
    db_name
}

.create_builtin_AIRR_human_germline_db <-
    function(fasta_store, destdir, only.if.missing=FALSE, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(only.if.missing), isTRUEorFALSE(verbose))
    db_name <- .form_builtin_AIRR_human_germline_db_name(fasta_store)
    db_path <- file.path(destdir, db_name)
    if (!(dir.exists(db_path) && only.if.missing))
        create_germline_db(fasta_store, IG_LOCI, db_path, verbose=verbose)
    if (basename(fasta_store) == "ref")
        add_V_ndm_data_to_germline_db(db_path, fasta_store,
                                      domain_system="imgt")
}

.create_missing_builtin_AIRR_human_germline_dbs <-
    function(human_dir, destdir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(human_dir), dir.exists(human_dir),
              isTRUEorFALSE(verbose))
    subdirs <- list.dirs(human_dir, full.names=FALSE, recursive=FALSE)
    subdirs <- setdiff(subdirs, "diffs")
    fasta_stores <- file.path(human_dir, rep(subdirs, each=2L), c("ref", "src"))
    for (fasta_store in fasta_stores)
        .create_builtin_AIRR_human_germline_db(fasta_store, destdir,
                                               only.if.missing=TRUE,
                                               verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_missing_builtin_AIRR_mouse_germline_dbs()
###

.form_builtin_AIRR_mouse_germline_db_name <- function(fasta_store)
{
    strain <- basename(fasta_store)
    loci_in1string <- paste(IG_LOCI, collapse="+")
    version <- read_version_file(fasta_store)
    sprintf("_AIRR.mouse.%s.%s.%s", strain, loci_in1string, version)
}

.create_builtin_AIRR_mouse_germline_db <-
    function(fasta_store, destdir, only.if.missing=FALSE, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(only.if.missing), isTRUEorFALSE(verbose))
    db_name <- .form_builtin_AIRR_mouse_germline_db_name(fasta_store)
    db_path <- file.path(destdir, db_name)
    if (!(dir.exists(db_path) && only.if.missing))
        create_germline_db(fasta_store, IG_LOCI, db_path, verbose=verbose)
}

.create_missing_builtin_AIRR_mouse_germline_dbs <-
    function(mouse_dir, destdir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(mouse_dir), dir.exists(mouse_dir),
              isTRUEorFALSE(verbose))
    fasta_stores <- list.dirs(mouse_dir, recursive=FALSE)
    for (fasta_store in fasta_stores)
        .create_builtin_AIRR_mouse_germline_db(fasta_store, destdir,
                                               only.if.missing=TRUE,
                                               verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_missing_builtin_AIRR_rhesus_monkey_germline_dbs()
###

.form_builtin_AIRR_rhesus_monkey_germline_db_name <- function(fasta_store)
{
    loci_in1string <- paste(IG_LOCI, collapse="+")
    version <- basename(fasta_store)
    sprintf("_AIRR.rhesus_monkey.%s.%s", loci_in1string, version)
}

.create_builtin_AIRR_rhesus_monkey_germline_db <-
    function(fasta_store, destdir, only.if.missing=FALSE, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(only.if.missing), isTRUEorFALSE(verbose))
    db_name <- .form_builtin_AIRR_rhesus_monkey_germline_db_name(fasta_store)
    db_path <- file.path(destdir, db_name)
    if (!(dir.exists(db_path) && only.if.missing))
        create_germline_db(fasta_store, IG_LOCI, db_path, verbose=verbose)
    add_V_ndm_data_to_germline_db(db_path, fasta_store, domain_system="imgt")
}

.create_missing_builtin_AIRR_rhesus_monkey_germline_dbs <-
    function(rhesus_monkey_dir, destdir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(rhesus_monkey_dir),
              dir.exists(rhesus_monkey_dir), isTRUEorFALSE(verbose))
    subdirs <- list.dirs(rhesus_monkey_dir, full.names=FALSE, recursive=FALSE)
    subdirs <- setdiff(subdirs, "diffs")
    fasta_stores <- file.path(rhesus_monkey_dir, subdirs)
    for (fasta_store in fasta_stores)
        .create_builtin_AIRR_rhesus_monkey_germline_db(fasta_store, destdir,
                                                       only.if.missing=TRUE,
                                                       verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_missing_builtin_germline_dbs()
###

### Not exported!
### 'destdir' must exist.
create_missing_builtin_germline_dbs <- function(destdir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(verbose))

    AIRR_germline_seq_dir <- system.file(package="igblastr",
                                 "extdata", "germline_sets", "AIRR",
                                 mustWork=TRUE)

    human_dir <- file.path(AIRR_germline_seq_dir, "human")
    .create_missing_builtin_AIRR_human_germline_dbs(human_dir, destdir,
                                                    verbose=verbose)

    mouse_dir <- file.path(AIRR_germline_seq_dir, "mouse")
    .create_missing_builtin_AIRR_mouse_germline_dbs(mouse_dir, destdir,
                                                    verbose=verbose)

    rhesus_monkey_dir <- file.path(AIRR_germline_seq_dir, "rhesus_monkey")
    .create_missing_builtin_AIRR_rhesus_monkey_germline_dbs(rhesus_monkey_dir,
                                                            destdir,
                                                            verbose=verbose)

    ## Any other built-in germline dbs to create?
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reset_germline_dbs()
###

### If 'destdir' exists, it gets destroyed and replaced with a freshly
### populated directory. In other words, its final content doesn't depend
### on whether it already exists or not.
.create_all_builtin_germline_dbs <- function(destdir, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    ## We first create the dbs in a temporary folder, and, only if successful,
    ## rename the temporary folder to 'destdir'. Otherwise we destroy the
    ## temporary folder and raise an error. This achieves atomicity.
    tmp_destdir <- tempfile("builtin_germline_dbs_")
    dir.create(tmp_destdir)
    on.exit(nuke_file(tmp_destdir))

    create_missing_builtin_germline_dbs(tmp_destdir, verbose=verbose)

    ## Everything went fine so we can rename 'tmp_destdir' to 'destdir'.
    rename_file(tmp_destdir, destdir, replace=TRUE)
}

### Will nuke any user-installed db!
reset_germline_dbs <- function(verbose=FALSE)
{
    set_db_in_use("germline", "")  # cancel current selection
    germline_dbs_home <- igblastr_cache(GERMLINE_DBS)
    .create_all_builtin_germline_dbs(germline_dbs_home, verbose=verbose)
}

