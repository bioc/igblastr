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

.get_final_fasta_path <- function(destdir, region_type)
{
    file.path(destdir, paste0(region_type, ".fasta"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .delete_repeated_fasta_records()
###

### In-place replacement!
### "Repeated" alleles are alleles with identical names **and** identical
### sequences. Note that:
### - At this point, the gaps in the allele sequences should have been removed.
### - We keep alleles with identical sequences but different names.
### - We also keep alleles with identical names but different sequences.
###   HOWEVER, we will disambiguate their names later with
###   .disambiguate_allele_names() (see below).
.delete_repeated_fasta_records <- function(fasta_file, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(fasta_file), isTRUEorFALSE(verbose))
    dna <- readDNAStringSet(fasta_file)
    what <- if (verbose) "FASTA record(s)" else ""
    from <- if (verbose) basename(fasta_file) else ""
    dna2 <- drop_repeated_sequences(dna, what=what, from=from)
    if (length(dna2) != length(dna))
        writeXStringSet(dna2, fasta_file)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .disambiguate_allele_names()
###

### Similar to base::make.unique() but mangles with suffixes made of
### lowercase letters.
.make_pool_of_suffixes <- function(min_pool_size)
{
    max_pool_size <- (length(letters)**8 - 1) / (length(letters) - 1) - 1
    if (min_pool_size > max_pool_size)
        stop(wmsg("too many duplicate seq ids"))
    ans <- character(0)
    for (i in 1:7) {
        ans <- c(ans, mkAllStrings(letters, i))
        if (length(ans) >= min_pool_size)
            return(ans)
    }
    ## Should never happen because we checked for this condition earlier (see
    ## above).
    stop(wmsg("too many duplicate seq ids"))
}

.make_unique_seqids <- function(seqids)
{
    stopifnot(is.character(seqids))
    if (length(seqids) <= 1L)
        return(seqids)
    oo <- order(seqids)
    seqids2 <- seqids[oo]
    ir <- IRanges(1L, runLength(Rle(seqids2)))
    pool_of_suffixes <- .make_pool_of_suffixes(max(width(ir)))
    suffixes <- extractList(pool_of_suffixes, ir)  # CharacterList
    suffixes[lengths(suffixes) == 1L] <- ""
    seqids2 <- paste0(seqids2, unlist(suffixes, use.names=FALSE))
    ans <- seqids2[S4Vectors:::reverseIntegerInjection(oo, length(oo))]
    setNames(ans, names(seqids))
}

### In-place replacement!
### Does not touch the sequences, only the allele names.
### Returns the number of alleles.
.disambiguate_allele_names <- function(fasta_file, verbose=FALSE)
{
    stopifnot(isSingleNonWhiteString(fasta_file), isTRUEorFALSE(verbose))
    fasta_lines <- readLines(fasta_file)
    header_idx <- grep("^>", fasta_lines)
    allele_names <- trimws2(sub("^>", "", fasta_lines[header_idx]))
    if (anyDuplicated(allele_names)) {
        new_allele_names <- .make_unique_seqids(allele_names)
        fasta_lines[header_idx] <- paste0(">", new_allele_names)
        writeLines(fasta_lines, fasta_file)
        if (verbose) {
            idx <- which(allele_names != new_allele_names)
            in1string <- paste0(allele_names[idx], "->", new_allele_names[idx],
                                collapse=", ")
            msg <- c("Renamed the ", length(idx), " following ",
                     "allele(s) in ", basename(fasta_file), ": ", in1string)
            message("  o ", wmsg(msg, margin=4L), "\n")
        }
    }
    length(allele_names)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .combine_and_edit_fasta_files()
###

### The workhorse behind create_region_db().
###
### See procedure described at
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### for how to create a germline or C-region db from the FASTA files
### available at IMGT.
### This is a 3-step procedure: (1) combine, (2) edit, (3) compile.
### The .combine_and_edit_fasta_files() function below only implements
### steps (1) and (2). Compilation (with makeblastdb) will happen at a
### latter time.
### IMPORTANT NOTE: .combine_and_edit_fasta_files() performs an **enhanced**
### edit step that consists of the 3 following sub-steps:
###   2a. basic edit as performed by the original edit_imgt_file.pl script;
###   2b. drop repeated alleles;
###   2c. disambiguate allele names.
### Finally note that the same procedure can be applied as-is to the FASTA
### files provided by AIRR-community/OGRDB.
### Returns the final number of alleles.
.combine_and_edit_fasta_files <- function(fasta_files, destdir,
                                          region_type=VDJC_REGION_TYPES,
                                          verbose=FALSE)
{
    .check_input_fasta_files(fasta_files)
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    region_type <- match.arg(region_type)
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    if (verbose) {
        in1string <- paste(basename(fasta_files), collapse=", ")
        msg <- c("Input file(s): ", in1string)
        message("  o ", wmsg(msg, margin=4L), "\n")
    }

    ## (1) Combine FASTA files.
    combined_fasta <- file.path(destdir, paste0(".", region_type, ".fasta"))
    concatenate_files(fasta_files, combined_fasta)

    ## (2a) Edit combined FASTA file. In igblastr 0.99.17, we switched
    ##      from edit_imgt_file() to redit_imgt_file() to perform this step.
    ##      This allowed us to no longer depend on Perl.
    final_fasta <- .get_final_fasta_path(destdir, region_type)
    #errfile <- file.path(destdir,
    #                     paste0(region_type, "_edit_imgt_file_errors.txt"))
    #edit_imgt_file(combined_fasta, final_fasta, errfile, check.output=TRUE)
    redit_imgt_file(combined_fasta, final_fasta)
    unlink(combined_fasta, force=TRUE)

    ## (2b) Drop repeated alleles.
    .delete_repeated_fasta_records(final_fasta, verbose=verbose)

    ## (2c) Mangle allele names to make them unique if they're not.
    ##      Because we did (2b), remaining repeated allele names are
    ##      guaranteed to be associated with distinct DNA sequences.
    .disambiguate_allele_names(final_fasta, verbose=verbose)
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
### See .combine_and_edit_fasta_files() above in this file for the workhorse
### behind create_region_db().
### 'destdir' must be the path to a writable directory that already exists!
### The following subdirectory and files will be added to 'destdir':
### - V_original_fasta/: subdirectory containing the input FASTA files
###       corresponding to the V regions, one FASTA file per region;
### - V.fasta: the combined and edited FASTA file produced by calling
###       .combine_and_edit_fasta_files() on the files in V_original_fasta/,
###       with allele names disambiguated if needed.
create_region_db <- function(fasta_files, destdir,
                             region_type=VDJC_REGION_TYPES,
                             overwrite=FALSE, verbose=FALSE)
{
    .check_input_fasta_files(fasta_files)
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    region_type <- match.arg(region_type)
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

    ## Combine and edit the original fasta files.
    original_fasta_files <- list_fasta_files(original_fasta_dir)
    num_alleles <- .combine_and_edit_fasta_files(original_fasta_files, destdir,
                                                 region_type=region_type,
                                                 verbose=verbose)
    if (verbose)
        message("... done (number of alleles in db: ", num_alleles, ").\n")
}

