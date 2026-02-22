### =========================================================================
### create_germline_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_loci_in_germline_fasta_dir()
###

.list_VDJ_fasta_files <- function(fasta_dir, loci_prefix)
{
    stopifnot(isSingleNonWhiteString(fasta_dir), dir.exists(fasta_dir),
              isSingleNonWhiteString(loci_prefix))
    pattern <- paste0("^", loci_prefix, ".[VDJ]\\.fasta$")
    fasta_files <- list_fasta_files(fasta_dir, pattern, full.names=FALSE)
    stopifnot(length(fasta_files) != 0L)
    fasta_files
}

### Returns a character vector of loci in canonical order.
.get_loci_from_germline_fasta_set <- function(fasta_files, loci_prefix)
{
    stopifnot(is.character(fasta_files),
              isSingleString(loci_prefix), loci_prefix %in% c("IG", "TR"))
    loci <- unique(sub("[VDJ]\\.fasta$", "", fasta_files))
    valid_loci <- if (loci_prefix == "IG") IG_LOCI else TR_LOCI
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

list_loci_in_germline_fasta_dir <-
    function(fasta_dir, loci_prefix, check.fasta.set=FALSE)
{
    stopifnot(isTRUEorFALSE(check.fasta.set))
    fasta_files <- .list_VDJ_fasta_files(fasta_dir, loci_prefix)
    loci <- .get_loci_from_germline_fasta_set(fasta_files, loci_prefix)
    if (check.fasta.set)
        .check_fasta_set(fasta_files, loci)
    loci
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .collect_fasta_files()
###

.list_fasta_files_for_region_type <- function(fasta_dir,
                                              region_type=VDJ_REGION_TYPES)
{
    region_type <- match.arg(region_type)
    pattern <- paste0(region_type, "\\.fasta$")
    fasta_files <- list_fasta_files(fasta_dir, pattern, full.names=FALSE)
    if (length(fasta_files) == 0L)
        stop(wmsg("Anomaly: no ", region_type, " files found in ", fasta_dir))
    fasta_files
}

.collect_fasta_files <- function(fasta_dir, region_type, loci)
{
    wanted_loci <- get_region_type_loci(region_type, loci)
    ## 'loci' should have gone thru .check_loci_for_missing_regions()
    ## so this is not supposed to happen. However it also went thru
    ## .get_effective_loci() which could have removed some loci from
    ## the original user selection. So yes, it's actually still possible
    ## that 'wanted_loci' will be empty!
    if (length(wanted_loci) == 0L)
        stop(wmsg("no fasta files found for region ", region_type, " for ",
                  "the selected loci"))
    wanted_files <- paste0(wanted_loci, region_type, ".fasta")
    found_files <- .list_fasta_files_for_region_type(fasta_dir, region_type)
    fasta_files <- intersect(wanted_files, found_files)
    if (length(fasta_files) == 0L)
        stop(wmsg("no fasta files found for region ", region_type, " for ",
                  "the selected loci"))
    missing_files <- setdiff(wanted_files, found_files)
    n <- length(missing_files)
    if (n != 0L) {
        verb <- if (n == 1L) " is" else "s are"
        in1string <- paste(missing_files, collapse=", ")
        warning(wmsg("the following file", verb, " missing ",
                     "for ", region_type, ": ", in1string), immediate.=TRUE)
    }
    file.path(fasta_dir, fasta_files)
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

### Create the three "region dbs": one V-, one D-, and one J-region db.
.create_VDJ_region_dbs <- function(fasta_dir, loci, destdir,
                                   gapped=FALSE, with.intdata=FALSE,
                                   verbose=FALSE)
{
    if (!isSingleNonWhiteString(fasta_dir))
        stop(wmsg("'fasta_dir' must be a single (non-empty) string"))
    if (!dir.exists(fasta_dir))
        stop(wmsg("directory ", fasta_dir, " not found"))
    stopifnot(isTRUEorFALSE(gapped))
    stopifnot(isTRUEorFALSE(with.intdata))
    stopifnot(isTRUEorFALSE(verbose))
    V_fasta_files <- .collect_fasta_files(fasta_dir, "V", loci)
    create_region_db(V_fasta_files, destdir, region_type="V",
                     gapped=gapped, with.intdata=with.intdata, verbose=verbose)
    D_fasta_files <- .collect_fasta_files(fasta_dir, "D", loci)
    create_region_db(D_fasta_files, destdir, region_type="D", verbose=verbose)
    J_fasta_files <- .collect_fasta_files(fasta_dir, "J", loci)
    create_region_db(J_fasta_files, destdir, region_type="J", verbose=verbose)
}

### A "germline db" is made of three "region dbs": one V-, one D-, and one
### J-region db. Calls create_region_db() to create each "region db".
### Note that 'destdir' will typically be the path to a subdir of the
### GERMLINE_DBS cache compartment (see R/cache-utils.R for details about
### igblastr's cache organization). This subdir or any of its parent
### directories don't need to exist yet.
### See create_region_db() in R/create_region_db.R for the roles of
### the 'gapped' and 'with.intdata' arguments.
create_germline_db <- function(fasta_dir, loci, destdir,
                               gapped=FALSE, with.intdata=FALSE,
                               force=FALSE, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir))
    if (!isTRUEorFALSE(gapped))
        stop(wmsg("'gapped' must be TRUE or FALSE"))
    check_with.intdata(with.intdata, gapped)
    if (!isTRUEorFALSE(force))
        stop(wmsg("'force' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    if (verbose) {
        what <- c("CREATING ", basename(destdir))
        message("\n====== START ", wmsg(what), " ======\n")
    }

    if (dir.exists(destdir) && !force)
        .stop_on_existing_germline_db(destdir)

    ## We first create the three region dbs in a temporary folder, and, only
    ## if successful, we replace 'destdir' with the temporary folder. Otherwise
    ## we destroy the temporary folder and raise an error. This achieves
    ## atomicity and avoids loosing the content of the existing 'destdir' in
    ## case something goes wrong.
    tmp_destdir <- tempfile("germline_db_")
    dir.create(tmp_destdir)
    on.exit(nuke_file(tmp_destdir))
    .create_VDJ_region_dbs(fasta_dir, loci, tmp_destdir,
                           gapped=gapped, with.intdata=with.intdata,
                           verbose=verbose)
    rename_file(tmp_destdir, destdir, replace=TRUE)

    if (verbose)
        message("====== DONE ", wmsg(what), " ======\n")
}

