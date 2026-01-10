### =========================================================================
### Creation of built-in germline and C-region dbs
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


read_version_file <- function(dirpath)
{
    stopifnot(isSingleNonWhiteString(dirpath))
    version_path <- file.path(dirpath, "version")
    if (!file.exists(version_path))
        stop(wmsg("missing 'version' file in ", dirpath, "/"))
    version <- readLines(version_path)
    if (length(version) != 1L)
        stop(wmsg("file '", version_path, "' should contain exactly one line"))
    version <- trimws2(version)
    if (version == "")
        stop(wmsg("file '", version_path, "' contains only white spaces"))
    version
}


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

.create_builtin_AIRR_human_germline_db <- function(fasta_store, destdir,
                                                   only.if.missing=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(only.if.missing))
    db_name <- .form_builtin_AIRR_human_germline_db_name(fasta_store)
    db_path <- file.path(destdir, db_name)
    if (!(only.if.missing && dir.exists(db_path)))
        create_germline_db(fasta_store, IG_LOCI, db_path)
}

.create_missing_builtin_AIRR_human_germline_dbs <- function(human_dir, destdir)
{
    stopifnot(isSingleNonWhiteString(human_dir), dir.exists(human_dir))
    subdirs <- list.dirs(human_dir, full.names=FALSE, recursive=FALSE)
    subdirs <- setdiff(subdirs, "diffs")
    fasta_stores <- file.path(human_dir, rep(subdirs, each=2L), c("ref", "src"))
    for (fasta_store in fasta_stores)
        .create_builtin_AIRR_human_germline_db(fasta_store, destdir,
                                               only.if.missing=TRUE)
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

.create_builtin_AIRR_mouse_germline_db <- function(fasta_store, destdir,
                                                   only.if.missing=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(only.if.missing))
    db_name <- .form_builtin_AIRR_mouse_germline_db_name(fasta_store)
    db_path <- file.path(destdir, db_name)
    if (!(only.if.missing && dir.exists(db_path)))
        create_germline_db(fasta_store, IG_LOCI, db_path)
}

.create_missing_builtin_AIRR_mouse_germline_dbs <- function(mouse_dir, destdir)
{
    stopifnot(isSingleNonWhiteString(mouse_dir), dir.exists(mouse_dir))
    fasta_stores <- list.dirs(mouse_dir, recursive=FALSE)
    for (fasta_store in fasta_stores)
        .create_builtin_AIRR_mouse_germline_db(fasta_store, destdir,
                                               only.if.missing=TRUE)
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
    function(fasta_store, destdir, only.if.missing=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(only.if.missing))
    db_name <- .form_builtin_AIRR_rhesus_monkey_germline_db_name(fasta_store)
    db_path <- file.path(destdir, db_name)
    if (!(only.if.missing && dir.exists(db_path)))
        create_germline_db(fasta_store, IG_LOCI, db_path)
}

.create_missing_builtin_AIRR_rhesus_monkey_germline_dbs <-
    function(rhesus_monkey_dir, destdir)
{
    stopifnot(isSingleNonWhiteString(rhesus_monkey_dir),
              dir.exists(rhesus_monkey_dir))
    subdirs <- list.dirs(rhesus_monkey_dir, full.names=FALSE, recursive=FALSE)
    subdirs <- setdiff(subdirs, "diffs")
    fasta_stores <- file.path(rhesus_monkey_dir, subdirs)
    for (fasta_store in fasta_stores)
        .create_builtin_AIRR_rhesus_monkey_germline_db(fasta_store, destdir,
                                                       only.if.missing=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_missing_builtin_germline_dbs()
###

### 'destdir' must exist.
create_missing_builtin_germline_dbs <- function(destdir)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir))

    AIRR_germline_seq_dir <- system.file(package="igblastr",
                                 "extdata", "germline_sequences", "AIRR",
                                 mustWork=TRUE)

    human_dir <- file.path(AIRR_germline_seq_dir, "human")
    .create_missing_builtin_AIRR_human_germline_dbs(human_dir, destdir)

    mouse_dir <- file.path(AIRR_germline_seq_dir, "mouse")
    .create_missing_builtin_AIRR_mouse_germline_dbs(mouse_dir, destdir)

    rhesus_monkey_dir <- file.path(AIRR_germline_seq_dir, "rhesus_monkey")
    .create_missing_builtin_AIRR_rhesus_monkey_germline_dbs(rhesus_monkey_dir,
                                                            destdir)

    ## Any other built-in germline dbs to create?
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_all_builtin_germline_dbs()
###

### 'destdir' must NOT exist.
create_all_builtin_germline_dbs <- function(destdir)
{
    stopifnot(isSingleNonWhiteString(destdir), !dir.exists(destdir))

    ## We first create the dbs in a temporary folder, and, only if successful,
    ## rename the temporary folder to 'destdir'. Otherwise we destroy the
    ## temporary folder and raise an error. This achieves atomicity.
    tmp_destdir <- tempfile("builtin_germline_dbs_")
    dir.create(tmp_destdir)
    on.exit(nuke_file(tmp_destdir))

    create_missing_builtin_germline_dbs(tmp_destdir)

    ## Everything went fine so we can rename 'tmp_destdir' to 'destdir'.
    rename_file(tmp_destdir, destdir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_builtin_IMGT_c_region_db()
###

.form_builtin_IMGT_c_region_db_name <- function(organism, loci, version)
{
    stopifnot(isSingleNonWhiteString(organism))
    sprintf("_IMGT.%s.%s.%s", organism, paste(loci, collapse="+"), version)
}

.create_builtin_IMGT_c_region_db <- function(fasta_store, loci_prefix,
                                             organism, destdir)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir))
    loci <- list_loci_in_c_region_fasta_dir(fasta_store, loci_prefix)
    version <- read_version_file(fasta_store)
    db_name <- .form_builtin_IMGT_c_region_db_name(organism, loci, version)
    db_path <- file.path(destdir, db_name)
    create_c_region_db(fasta_store, loci, db_path)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_all_builtin_c_region_dbs()
###

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

create_all_builtin_c_region_dbs <- function(destdir)
{
    stopifnot(isSingleNonWhiteString(destdir))

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
                                         tmp_destdir)
        fasta_store <- path_to_IMGT_c_region_fasta_store(organism, "TR")
        if (dir.exists(fasta_store))
            .create_builtin_IMGT_c_region_db(fasta_store, "TR", organism,
                                             tmp_destdir)
    }

    ## Any other built-in C-region dbs to create?

    ## Everything went fine so we can rename 'tmp_destdir' to 'destdir'.
    rename_file(tmp_destdir, destdir, replace=TRUE)
}

