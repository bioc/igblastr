### =========================================================================
### create_germline_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


.list_VDJ_fasta_files <- function(fasta_dir, tcr.db=FALSE)
{
    stopifnot(isSingleNonWhiteString(fasta_dir), dir.exists(fasta_dir))
    if (!isTRUEorFALSE(tcr.db))
        stop(wmsg("'tcr.db' must be TRUE or FALSE"))
    prefix <- if (tcr.db) "TR" else "IG"
    pattern <- paste0("^", prefix, ".[VDJ]\\.fasta$")
    fasta_files <- list.files(fasta_dir, pattern=pattern)
    stopifnot(length(fasta_files) != 0L)
    fasta_files
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_loci_from_input_germline_fasta_set()
###

### Returns a character vector of loci in canonical order.
.get_loci_from_germline_fasta_set <- function(fasta_files, tcr.db=FALSE)
{
    stopifnot(is.character(fasta_files), isTRUEorFALSE(tcr.db))
    loci <- unique(sub("[VDJ]\\.fasta$", "", fasta_files))
    valid_loci <- if (tcr.db) TR_LOCI else IG_LOCI
    stopifnot(all(loci %in% valid_loci))
    valid_loci[valid_loci %in% loci]  # return loci in canonical order
}

.check_fasta_set <- function(fasta_files, loci)
{
    stopifnot(is.character(fasta_files))
    loci2regiontypes <- map_loci_to_region_types(loci)
    for (locus in loci) {
        pattern <- paste0("^", locus)
        current_files <- grep(pattern, fasta_files, value=TRUE)
        expected_files <- paste0(locus, loci2regiontypes[[locus]], ".fasta")
        missing_files <- setdiff(expected_files, current_files)
        n <- length(missing_files)
        if (n != 0L) {
            verb <- if (n == 1L) " is" else "s are"
            in1string <- paste(missing_files, collapse=", ")
            warning(wmsg("the following file", verb, " missing ",
                         "for locus ", locus, ": ", in1string))
        }
        unexpected_files <- setdiff(current_files, expected_files)
        n <- length(unexpected_files)
        if (n != 0L) {
            verb <- if (n == 1L) " is" else "s are"
            in1string <- paste(unexpected_files, collapse=", ")
            warning(wmsg("the following file", verb, " usually not expected ",
                         "for locus ", locus, ": ", in1string))
        }
    }
}

get_loci_from_input_germline_fasta_set <-
    function(fasta_dir, tcr.db=FALSE, check.fasta.set=FALSE)
{
    stopifnot(isTRUEorFALSE(check.fasta.set))
    fasta_files <- .list_VDJ_fasta_files(fasta_dir, tcr.db=tcr.db)
    loci <- .get_loci_from_germline_fasta_set(fasta_files, tcr.db=tcr.db)
    if (check.fasta.set)
        .check_fasta_set(fasta_files, loci)
    loci
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_germline_db()
###

.stop_on_existing_germline_db <- function(destdir)
{
    db_name <- basename(destdir)
    msg1 <- c("Germline db ", db_name, " is already installed.")
    msg2 <- c("Use list_germline_dbs() to list the germline databases ",
              "already installed on your machine (see '?list_germline_dbs').")
    msg3 <- c("Use 'force=TRUE' to reinstall.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

.list_fasta_files <- function(fasta_dir, region_type=VDJ_REGION_TYPES)
{
    region_type <- match.arg(region_type)
    pattern <- paste0(region_type, "\\.fasta$")
    files <- list.files(fasta_dir, pattern=pattern)
    if (length(files) == 0L)
        stop(wmsg("Anomaly: no ", region_type, " files found in ", fasta_dir))
    file.path(fasta_dir, files)
}

### Create the three "region dbs": one V-, one D-, and one J-region db.
.create_VDJ_region_dbs <- function(fasta_dir, destdir, tcr.db,
                                   edit_imgt_file_Perl_script)
{
    for (region_type in VDJ_REGION_TYPES) {
        fasta_files <- .list_fasta_files(fasta_dir, region_type)
        create_region_db(fasta_files, destdir, region_type=region_type,
                         edit_imgt_file_Perl_script=edit_imgt_file_Perl_script)
    }
}

### Perl required!
###
### A "germline db" is made of three "region dbs": one V-, one D-, and one
### J-region db. Calls create_region_db() to create each "region db".
### Note that 'destdir' will typically be the path to a subdir of the
### GERMLINE_DBS cache compartment (see R/cache-utils.R for details about
### igblastr's cache organization). This subdir or any of its parent
### directories don't need to exist yet.
create_germline_db <- function(fasta_dir, destdir, tcr.db=FALSE, force=FALSE)
{
    if (!isTRUEorFALSE(tcr.db))
        stop(wmsg("'tcr.db' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    stopifnot(isSingleNonWhiteString(destdir))
    if (dir.exists(destdir) && !force)
        .stop_on_existing_germline_db(destdir)
    ## get_edit_imgt_file_Perl_script() also checks that Perl script
    ## edit_imgt_file.pl is available and that Perl is functioning.
    edit_imgt_file_Perl_script <- get_edit_imgt_file_Perl_script()

    ## We ignore the returned loci. Only purpose is to check the set
    ## of input germline FASTA files.
    get_loci_from_input_germline_fasta_set(fasta_dir, tcr.db=tcr.db,
                                           check.fasta.set=TRUE)

    ## We first create the three region dbs in a temporary folder, and, only
    ## if successful, we replace 'destdir' with the temporary folder. Otherwise
    ## we destroy the temporary folder and raise an error. This achieves
    ## atomicity and avoids loosing the content of the existing 'destdir' in
    ## case something goes wrong.
    tmp_destdir <- tempfile("germline_db_")
    dir.create(tmp_destdir, recursive=TRUE)
    on.exit(nuke_file(tmp_destdir))
    .create_VDJ_region_dbs(fasta_dir, tmp_destdir, tcr.db,
                           edit_imgt_file_Perl_script)
    rename_file(tmp_destdir, destdir, replace=TRUE)
}

