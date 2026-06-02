### =========================================================================
### install_IMGT_germline_db()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_effective_loci()
###

### Returns the intersection between 'wanted_loci' and 'found_loci'
### but with lots of sanity checks, bells and whistles.
.get_effective_loci <- function(wanted_loci, found_loci)
{
    stopifnot(is.character(wanted_loci), is.character(found_loci))
    keep_idx <- which(wanted_loci %in% found_loci)
    if (length(keep_idx) == 0L) {
        what <- if (length(wanted_loci) == 1L) "locus" else "loci"
        in1string <- paste0(wanted_loci, collapse=", ")
        stop(wmsg("no FASTA files found for ", what, " ", in1string))
    }
    missing_loci <- wanted_loci[-keep_idx]
    ## Like 'intersect(found_loci, wanted_loci)' but the returned
    ## intersection is guaranteed to be ordered like in 'found_loci'.
    loci <- found_loci[found_loci %in% wanted_loci]
    if (length(missing_loci) != 0L) {
        what1 <- if (length(missing_loci) == 1L) "locus" else "loci"
        in1string1 <- paste0(missing_loci, collapse=", ")
        what2 <- if (length(loci) == 1L) "locus" else "loci"
        in1string2 <- paste0(loci, collapse=", ")
        warning(wmsg("No FASTA files found for ", what1, " ", in1string1, " ",
                     "--> Installing germline db for ", what2, " ", in1string2,
                     "."),
                immediate.=TRUE)
    }
    loci
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_IMGT_germline_db()
###

.form_IMGT_germline_db_name <- function(fasta_store, loci)
{
    stopifnot(isSingleNonWhiteString(fasta_store), dir.exists(fasta_store))
    check_selected_loci(loci)
    organism_path <- dirname(fasta_store)
    organism <- basename(organism_path)
    refdir <- dirname(organism_path)
    stopifnot(basename(refdir) == VQUEST_REFERENCE_DIRECTORY)
    IMGT_store <- dirname(refdir)
    release <- basename(IMGT_store)
    sprintf("IMGT-%s.%s.%s", release, organism, paste(loci, collapse="+"))
}

### Why does IMGT human J allele IGLJ2A*01 have a codon start set to 1?
### This is unexpected because:
### - there's no FGXG motif in this coding frame (the allele sequence
###   translates to CGIRRRDQADRP in this coding frame);
### - IGLJ2A*01 sequence is exactly the same as IGLJ2*01 and IGLJ3*01
###   sequences, which both have a codon start set to 2;
### - IGLJ2*01 and IGLJ3*01 actually contain the FGXG motif in the coding
###   frame that starts at position 2 (the 12 codons in the sequence
###   translate to VVFGGGTKLTVL).
### TODO: Ask the IMGT folks about this.
### In the mean time, we exclude it. This shouldn't make any difference
### from an IgBLAST operations point of view because its sequence is
### the same as IGLJ2*01 and IGLJ3*01 which we keep (note for that matter
### that we could also exclude IGLJ3*01 and it shouldn't make any difference
### either).
.EXCLUDED_IMGT_HUMAN_J_ALLELES <- "IGLJ2A*01"

### Loads IgBLAST auxdata for the organism associated with the germline db.
### Will be used to "complete" the computed auxiliary data and verify
### agreement with it.
.load_igblast_auxdata <- function(db_name)
{
    igblast_organism <- infer_igblast_organism_from_db_name(db_name)
    if (is.na(igblast_organism))
        return(NULL)
    auxdata <- load_auxdata(igblast_organism)
    ## We know that the human auxdata shipped with IgBLAST is incorrect
    ## for IGHJ6*02 and IGHJ6*03 -- and so is discordant with what
    ## install_germline_db() will compute -- so we correct these rows.
    ## Note that if the user updated their "live" IgBLAST data with
    ## update_live_igdata(), then these rows should already be correct
    ## in 'auxdata'.
    if (igblast_organism == "human") {
        fixme <- auxdata[ , "allele_name"] %in% c("IGHJ6*02", "IGHJ6*03")
        auxdata[fixme, "extra_bps"] <- 1L  # replace 0L with 1L
    }
    ## We also remove rows for which the "coding_frame_start" / "cdr3_end"
    ## combination doesn't make sense. This is the case for example for
    ## human allele TRAJ31*01 and for mouse alleles TRAJ21*02, TRAJ24*02,
    ## and TRDJ2*02 (as of May 31, 2026).
    coding_frame_start <- auxdata[ , "coding_frame_start"]
    cdr3_end <- auxdata[ , "cdr3_end"]
    drop_idx <- which(((cdr3_end + 1L - coding_frame_start) %% 3L) != 0L)
    if (length(drop_idx) != 0L)
        auxdata <- auxdata[-drop_idx, ]
    auxdata
}

.check_concordance_with_igblast_intdata <- function(db_name)
{
    igblast_organism <- infer_igblast_organism_from_db_name(db_name)
    if (is.na(igblast_organism))
        return()
    computed_intdata <- load_intdata(db_name)
    igblast_intdata <- load_intdata(igblast_organism)
    disc_rowpairs <- find_discordant_intdata(computed_intdata, igblast_intdata)
    if (nrow(disc_rowpairs) == 0L)
        return()
    msg1 <- c("The computed \"internal data\" that we included in the ",
              "germline db has some disagreements with the \"internal data\" ",
              "for ", igblast_organism, " shipped with IgBLAST.")
    msg2 <- c("To display the disagreements, run:")
    msg3 <- c("  show_intdata_disagreements(\"", db_name, "\")")
    warning(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", msg3)
}

install_IMGT_germline_db <- function(release, organism="Homo sapiens",
                                     tcr.db=FALSE, loci="auto",
                                     without.intdata=FALSE,
                                     without.auxdata=FALSE,
                                     overwrite=FALSE, verbose=FALSE, ...)
{
    ## Check arguments.
    organism <- normalize_IMGT_organism(organism)
    loci <- normalize_user_supplied_loci(loci, tcr.db=tcr.db,
                                         stop.if.missing.regions=TRUE)
    loci_prefix <- extract_selected_loci_prefix(loci)
    if (!isTRUEorFALSE(without.intdata))
        stop(wmsg("'without.intdata' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(without.auxdata))
        stop(wmsg("'without.auxdata' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    ## Download IMGT/V-QUEST release to local store if it's not already there.
    IMGT_store <- download_IMGT_release_to_IMGT_store(release, ...)

    ## Get path to local FASTA store for IMGT germline sequences.
    fasta_store <- find_organism_fasta_store_in_IMGT_store(IMGT_store,
                                                           organism,
                                                           loci_prefix)

    ## Keep loci for which IMGT actually provides FASTA files.
    found_loci <- list_loci_in_germline_fasta_dir(fasta_store, loci_prefix)
    loci <- .get_effective_loci(loci, found_loci)

    ## Compute 'db_name'.
    db_name <- .form_IMGT_germline_db_name(fasta_store, loci)

    ## Do we need to exclude any J allele known to be problematic?
    if (organism == "Homo_sapiens" && loci_prefix == "IG") {
        excluded_J_alleles <- .EXCLUDED_IMGT_HUMAN_J_ALLELES
    } else {
        excluded_J_alleles <- character(0)
    }

    ## Create and install germline db.
    install_dir <- get_germline_dbs_home(TRUE)  # guaranteed to exist
    with.auxdata <- !without.auxdata && (loci_prefix == "IG" ||
                                         organism == "Homo_sapiens")
    igblast_auxdata <- if (with.auxdata) .load_igblast_auxdata(db_name)
                       else NULL
    if.exists <- if (overwrite) "overwrite" else "error"
    install_germline_db(install_dir, db_name, fasta_store, loci,
                        gapped=TRUE, with.intdata=!without.intdata,
                        excluded_J_alleles=excluded_J_alleles,
                        with.auxdata=with.auxdata, imgt.fasta=TRUE,
                        ref_auxdata=igblast_auxdata,
                        if.exists=if.exists, verbose=verbose,
                        cheer.if.success=TRUE)

    ## Check concordance of computed intdata with IgBLAST intdata.
    if (!without.intdata)
        .check_concordance_with_igblast_intdata(db_name)

    invisible(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_IMGT_c_region_db()
###

.form_IMGT_c_region_db_name <- function(organism, loci, version)
{
    stopifnot(isSingleNonWhiteString(organism))
    sprintf("IMGT.%s.%s.%s", organism, paste(loci, collapse="+"), version)
}

install_IMGT_c_region_db <- function(organism, loci,
                                     disambiguate.allele.names=FALSE,
                                     overwrite=FALSE, verbose=FALSE)
{
    organism <- normalize_IMGT_organism(organism)
    igblast_organism <- lookup_igblast_organism(organism)
    loci <- normalize_user_supplied_loci(loci)
    loci_prefix <- extract_selected_loci_prefix(loci)
    if (!isTRUEorFALSE(disambiguate.allele.names))
        stop(wmsg("'disambiguate.allele.names' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    ## Get path to local FASTA store for IMGT C-region sequences.
    fasta_store <-
        path_to_IMGT_c_region_fasta_store(igblast_organism, loci_prefix)
    if (!dir.exists(fasta_store))
        stop(wmsg("we're not aware of any ", loci_prefix, " C-region ",
                  "sequences available at IMGT for ", igblast_organism,
                  ", sorry"))

    ## Keep loci for which IMGT actually provides FASTA files.
    found_loci <- list_loci_in_c_region_fasta_dir(fasta_store, loci_prefix)
    loci <- .get_effective_loci(loci, found_loci)

    ## Compute 'db_name'.
    version <- read_version_file(fasta_store)
    db_name <- .form_IMGT_c_region_db_name(igblast_organism, loci, version)

    ## Create IMGT C-region db.
    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    db_path <- file.path(c_region_dbs_home, db_name)
    create_c_region_db(fasta_store, loci, db_path,
                       disambiguate.allele.names=disambiguate.allele.names,
                       overwrite=overwrite, verbose=verbose)

    ## Success!
    message("C-region db ", db_name, " successfully installed.")
    message("Call use_c_region_db(\"", db_name, "\") to select it as the")
    message("C-region db to use with igblastn().")

    invisible(db_name)
}

