### =========================================================================
### create_region_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


.check_input_fasta_files <- function(fasta_files)
{
    if (!is.character(fasta_files) || length(fasta_files) == 0L)
        stop(wmsg("'fasta_files' must be a non-empty character vector"))
    if (anyNA(fasta_files) || anyDuplicated(fasta_files))
        stop(wmsg("'fasta_files' cannot contain NAs or duplicates"))
}

.check_destdir <- function(destdir)
{
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
}

.check_gapped <- function(gapped, region_type)
{
    if (!isTRUEorFALSE(gapped))
        stop(wmsg("'gapped' must be TRUE or FALSE"))
    if (gapped && !identical(region_type, "V"))
        stop(wmsg("'gapped=TRUE' can only be used when 'region_type' is \"V\""))
}

.get_final_fasta_path <- function(destdir, region_type)
{
    .check_destdir(destdir)
    file.path(destdir, paste0(region_type, ".fasta"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .load_original_V_alleles_with_loci()
###

### The names of the supplied files must be of the form <locus>V.fasta
### where locus is a 3-letter string. Note that we don't care if <locus>
### is a valid locus or not (see IG_LOCI and TR_LOCI in R/loci-utils.R for
### the 7 valid loci), as long as it's a 3-letter string, so we don't bother
### to check.
.infer_loci_from_filenames <- function(fasta_files)
{
    .check_input_fasta_files(fasta_files)
    filenames <- basename(fasta_files)
    suffix <- substr(filenames, 4L, nchar(filenames))
    if (!all(suffix == "V.fasta")) {
        in1string <- paste(filenames, collapse=", ")
        stop(wmsg("The names of the FASTA files must be of the form ",
                  "<locus>V.fasta, where <locus> is a 3-letter string. ",
                  "Got: ", in1string))
    }
    substr(filenames, 1L, 3L)
}

.load_original_V_alleles_with_loci <- function(fasta_files)
{
    loci <- .infer_loci_from_filenames(fasta_files)
    dna <- lapply(seq_along(fasta_files),
        function(i) {
            dna <- readDNAStringSet(fasta_files[[i]])
            mcols(dna)$locus <- loci[[i]]
            dna
        }
    )
    do.call(c, dna)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .clean_and_merge_fasta_files()
###
### See procedure described at
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### for how to create a germline or C-region db from the FASTA files
### available at IMGT.
### This is a 3-step procedure: (1) combine, (2) clean/edit, (3) compile.
### The .clean_and_merge_fasta_files() function only performs steps (1)
### and (2). Compilation (with makeblastdb) will happen at a latter time.
### Note that .clean_and_merge_fasta_files() calls clean_allele_set() to
### perform step (2). See clean_allele_set() in R/clean_allele_set.R for
### more information.

### See create_region_db() below in this file for the roles of
### the 'gapped' and 'with.intdata' arguments.
### Returns a DNAStringSet object containing the "clean" alleles
### that went into the db.
.clean_and_merge_fasta_files <- function(fasta_files, destdir,
                                         region_type=VDJC_REGION_TYPES,
                                         gapped=FALSE, with.intdata=FALSE,
                                         verbose=FALSE)
{
    .check_input_fasta_files(fasta_files)
    region_type <- match.arg(region_type)
    final_fasta <- .get_final_fasta_path(destdir, region_type)
    .check_gapped(gapped, region_type)
    check_with.intdata(with.intdata, gapped)
    stopifnot(isTRUEorFALSE(verbose))

    if (verbose) {
        in1string <- paste(basename(fasta_files), collapse=", ")
        msg <- c("Input file(s): ", in1string)
        message("  o ", wmsg(msg, margin=4L), "\n")
    }

    ## (1) Combine.
    if (with.intdata) {
        ## We use .load_original_V_alleles_with_loci() to get the locus info
        ## (returned as a metadata column on 'dna'). clean_allele_set()
        ## requires and uses this metadata column when 'with.intdata' is TRUE.
        dna <- .load_original_V_alleles_with_loci(fasta_files)
    } else {
        dna <- readDNAStringSet(fasta_files)
    }

    ## (2) Clean/edit.
    dna <- clean_allele_set(dna, gapped=gapped, with.intdata=with.intdata,
                            verbose=verbose)

    writeXStringSet(dna, final_fasta)
    dna
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .add_intdata_to_db()
###

.add_intdata_to_db <- function(dna, destdir)
{
    stopifnot(is(dna, "DNAStringSet"))
    allele_names <- names(dna)
    stopifnot(!is.null(allele_names))
    dna_mcols <- mcols(dna, use.names=FALSE)
    stopifnot(is(dna_mcols, "DataFrame"),
              identical(allele_names, dna_mcols[ , "allele_name"]))

    .check_destdir(destdir)
    internal_data_path <- file.path(destdir, "internal_data")
    stopifnot(!dir.exists(internal_data_path))

    V_ndm_data <- as.data.frame(dna_mcols[ , names(IGBLAST_INTDATA_COL2CLASS)])
    check_V_ndm_data_col2class(V_ndm_data)
    dir.create(internal_data_path)
    destfile <- file.path(internal_data_path, "V.ndm.imgt")
    write_V_ndm_data(V_ndm_data, destfile)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_region_db()
###

.region_db_already_exists <- function(destdir, region_type)
{
    original_fasta_dir <- get_db_original_fasta_dir(destdir, region_type)
    final_fasta <- .get_final_fasta_path(destdir, region_type)
    file.exists(original_fasta_dir) || file.exists(final_fasta)
}

.stop_on_existing_region_db <- function(destdir, region_type)
{
    msg1 <- c("There already seems to be a ", region_type, "-region ",
              "database in ", destdir)
    msg2 <- c("Use 'overwrite=TRUE' to overwrite or choose another ",
              "destination directory.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2))
}

.nuke_existing_region_db <- function(destdir, region_type)
{
    pattern <- paste0("^", region_type, "\\.fasta$")
    clean_blastdbs(destdir, pattern)
    original_fasta_dir <- get_db_original_fasta_dir(destdir, region_type)
    final_fasta <- .get_final_fasta_path(destdir, region_type)
    nuke_file(original_fasta_dir)
    nuke_file(final_fasta)
}

### Creates a "region db" (V-, D-, J-, or C-region) from a collection of
### FASTA files (typically obtained from IMGT or AIRR-community/OGRDB) for
### a given organism.
### See .clean_and_merge_fasta_files() above in this file for the workhorse
### behind create_region_db().
### 'destdir' must be the path to a writable directory that already exists!
### Set 'gapped' to TRUE if the supplied sequences are gapped V allele
### sequences (note that only V allele sequences are allowed to have gaps).
### Set 'with.intdata' to TRUE if the supplied sequences are gapped V allele
### sequences **and** the associated internal data should be computed and
### added to the db.
### The following subdirectory and files will be added to 'destdir':
### - V_original_fasta/: subdirectory containing the input FASTA files
###       corresponding to the V regions, one FASTA file per region;
### - V.fasta: FASTA file produced by calling .clean_and_merge_fasta_files()
###       on the files in V_original_fasta/, with allele names disambiguated
###       if needed;
### - internal_data/V.ndm.imgt: internal data associated with the V alleles,
###       when 'with.intdata' is set to TRUE.
create_region_db <- function(fasta_files, destdir,
                             region_type=VDJC_REGION_TYPES,
                             gapped=FALSE, with.intdata=FALSE,
                             overwrite=FALSE, verbose=FALSE)
{
    .check_input_fasta_files(fasta_files)
    .check_destdir(destdir)
    region_type <- match.arg(region_type)
    .check_gapped(gapped, region_type)
    check_with.intdata(with.intdata, gapped)
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    if (verbose)
        message("Creating ", region_type, "-region db ...\n")

    if (.region_db_already_exists(destdir, region_type)) {
        if (!overwrite)
            .stop_on_existing_region_db(destdir, region_type)
        .nuke_existing_region_db(destdir, region_type)
    }

    ## Create "original fasta" subdir and copy fasta files to it.
    original_fasta_dir <- get_db_original_fasta_dir(destdir, region_type)
    stopifnot(dir.create(original_fasta_dir))
    destfiles <- names(fasta_files)
    if (is.null(destfiles)) {
        stopifnot(all(file.copy(fasta_files, original_fasta_dir)))
    } else {
        destfiles <- file.path(original_fasta_dir, destfiles)
        stopifnot(all(file.copy(fasta_files, destfiles)))
    }

    ## Clean and merge the original fasta files.
    original_fasta_files <- list_fasta_files(original_fasta_dir)
    dna <- .clean_and_merge_fasta_files(original_fasta_files, destdir,
                                        region_type=region_type,
                                        gapped=gapped,
                                        with.intdata=with.intdata,
                                        verbose=verbose)
    if (verbose)
        message("... done. (Final number of ", region_type, " alleles ",
                "in db = ", length(dna), ")\n")

    if (with.intdata) {
        if (verbose)
            message(wmsg("Adding the intdata to the db"), " ... ",
                    appendLF=FALSE)
        .add_intdata_to_db(dna, destdir)
        if (verbose)
            message("ok.\n")
    }

    dna
}

