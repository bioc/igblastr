### =========================================================================
### Creation of built-in germline and C-region dbs
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_missing_builtin_AIRR_human_germline_dbs()
###

.form_builtin_AIRR_human_germline_db_name <- function(fasta_dir)
{
    loci <- get_loci_from_input_germline_fasta_set(fasta_dir)
    flavor <- basename(fasta_dir)
    stopifnot(flavor %in% c("ref", "src"))
    version <- basename(dirname(fasta_dir))
    db_name <- sprintf("_AIRR.human.%s.%s", paste(loci, collapse="+"), version)
    if (flavor == "src")
        db_name <- paste0(db_name, ".src")
    db_name
}

.create_builtin_AIRR_human_germline_db <- function(fasta_dir, destdir,
                                                   only.if.missing=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(only.if.missing))
    db_name <- .form_builtin_AIRR_human_germline_db_name(fasta_dir)
    db_path <- file.path(destdir, db_name)
    if (!(only.if.missing && dir.exists(db_path)))
        create_germline_db(fasta_dir, db_path)
}

.create_missing_builtin_AIRR_human_germline_dbs <- function(human_dir, destdir)
{
    stopifnot(isSingleNonWhiteString(human_dir), dir.exists(human_dir))
    human_subdirs <- list.dirs(human_dir, recursive=FALSE)
    fasta_dirs <- file.path(rep(human_subdirs, each=2L), c("ref", "src"))
    for (fasta_dir in fasta_dirs) {
        .create_builtin_AIRR_human_germline_db(fasta_dir, destdir,
                                               only.if.missing=TRUE)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_missing_builtin_AIRR_mouse_germline_dbs()
###

.form_builtin_AIRR_mouse_germline_db_name <- function(fasta_dir)
{
    loci <- get_loci_from_input_germline_fasta_set(fasta_dir)
    strain <- basename(fasta_dir)
    version <- read_version_file(fasta_dir)
    sprintf("_AIRR.mouse.%s.%s.%s", strain, paste(loci, collapse="+"), version)
}

.create_builtin_AIRR_mouse_germline_db <- function(fasta_dir, destdir,
                                                   only.if.missing=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(only.if.missing))
    db_name <- .form_builtin_AIRR_mouse_germline_db_name(fasta_dir)
    db_path <- file.path(destdir, db_name)
    if (!(only.if.missing && dir.exists(db_path)))
        create_germline_db(fasta_dir, db_path)
}

.create_missing_builtin_AIRR_mouse_germline_dbs <- function(mouse_dir, destdir)
{
    fasta_dirs <- list.dirs(mouse_dir, recursive=FALSE)
    for (fasta_dir in fasta_dirs) {
        .create_builtin_AIRR_mouse_germline_db(fasta_dir, destdir,
                                               only.if.missing=TRUE)
    }
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
    dir.create(tmp_destdir, recursive=TRUE)
    on.exit(nuke_file(tmp_destdir))

    create_missing_builtin_germline_dbs(tmp_destdir)

    ## Everyting went fine so we can rename 'tmp_destdir' to 'destdir'.
    rename_file(tmp_destdir, destdir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .create_builtin_IMGT_c_region_db()
###

.form_builtin_IMGT_c_region_db_name <- function(organism, fasta_dir,
                                                tcr.db=FALSE)
{
    stopifnot(isSingleNonWhiteString(organism))
    loci <- get_loci_from_input_c_region_fasta_set(fasta_dir, tcr.db=tcr.db)
    version <- read_version_file(fasta_dir)
    sprintf("_IMGT.%s.%s.%s", organism, paste(loci, collapse="+"), version)
}

.create_builtin_IMGT_c_region_db <- function(fasta_dir, organism, destdir,
                                             tcr.db=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir))
    db_name <- .form_builtin_IMGT_c_region_db_name(organism, fasta_dir,
                                                   tcr.db=tcr.db)
    db_path <- file.path(destdir, db_name)
    create_c_region_db(fasta_dir, db_path, tcr.db=tcr.db)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_all_builtin_c_region_dbs()
###

create_all_builtin_c_region_dbs <- function(destdir)
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
        organism_dir <- file.path(IMGT_c_region_dir, organism)
        fasta_dir <- file.path(organism_dir, "IG", "14.1")
        .create_builtin_IMGT_c_region_db(fasta_dir, organism, tmp_destdir)
        fasta_dir <- file.path(organism_dir, "TR")
        if (dir.exists(fasta_dir))
            .create_builtin_IMGT_c_region_db(fasta_dir, organism, tmp_destdir,
                                             tcr.db=TRUE)
    }

    ## Any other built-in C-region dbs to create?

    ## Everyting went fine so we can rename 'tmp_destdir' to 'destdir'.
    rename_file(tmp_destdir, destdir, replace=TRUE)
}

