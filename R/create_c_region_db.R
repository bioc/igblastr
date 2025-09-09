### =========================================================================
### create_c_region_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


.list_C_fasta_files <- function(fasta_dir, tcr.db=FALSE)
{
    stopifnot(isSingleNonWhiteString(fasta_dir), dir.exists(fasta_dir))
    if (!isTRUEorFALSE(tcr.db))
        stop(wmsg("'tcr.db' must be TRUE or FALSE"))
    prefix <- if (tcr.db) "TR" else "IG"
    pattern <- paste0("^", prefix, ".C\\.fasta$")
    fasta_files <- list.files(fasta_dir, pattern=pattern)
    stopifnot(length(fasta_files) != 0L)
    fasta_files
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_loci_from_input_c_region_fasta_set()
###

### Returns a character vector of loci in canonical order.
.get_loci_from_c_region_fasta_set <- function(fasta_files, tcr.db=FALSE)
{
    stopifnot(is.character(fasta_files), isTRUEorFALSE(tcr.db))
    loci <- unique(sub("C\\.fasta$", "", fasta_files))
    valid_loci <- if (tcr.db) TR_LOCI else IG_LOCI
    stopifnot(all(loci %in% valid_loci))
    valid_loci[valid_loci %in% loci]  # return loci in canonical order
}

get_loci_from_input_c_region_fasta_set <-
    function(fasta_dir, tcr.db=FALSE)
{
    fasta_files <- .list_C_fasta_files(fasta_dir, tcr.db=tcr.db)
    .get_loci_from_c_region_fasta_set(fasta_files, tcr.db=tcr.db)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_c_region_db()
###

.stop_on_existing_c_region_db <- function(destdir)
{
    db_name <- basename(destdir)
    msg1 <- c("C-region db ", db_name, " is already installed.")
    msg2 <- c("Use list_c_region_dbs() to list the C-region databases ",
              "already installed on your machine (see '?list_c_region_dbs').")
    msg3 <- c("Use 'force=TRUE' to reinstall.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

### Creates a C-region db (constant regions) from a collection of FASTA
### files (typically obtained from IMGT) for a given organism.
### Note that 'destdir' will typically be the path to a subdir of the
### C_REGION_DBS cache compartment (see R/cache-utils.R for details about
### igblastr's cache organization). This subdir or any of its parent
### directories don't need to exist yet.
create_c_region_db <- function(fasta_dir, destdir, tcr.db=FALSE, force=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir))
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    if (dir.exists(destdir) && !force)
        .stop_on_existing_c_region_db(destdir)

    fasta_files <- .list_C_fasta_files(fasta_dir, tcr.db=tcr.db)
    fasta_files <- file.path(fasta_dir, fasta_files)

    ## We first create the db in a temporary folder, and, only if successful,
    ## replace 'destdir' with the temporary folder. Otherwise we destroy the
    ## temporary folder and raise an error. This achieves atomicity and avoids
    ## loosing the content of the existing 'destdir' in case something goes
    ## wrong.
    tmp_destdir <- tempfile("c_region_db_")
    dir.create(tmp_destdir, recursive=TRUE)
    on.exit(nuke_file(tmp_destdir))
    create_region_db(fasta_files, tmp_destdir, region_type="C")
    rename_file(tmp_destdir, destdir, replace=TRUE)
}

