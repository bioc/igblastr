### =========================================================================
### create_region_db()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


.get_original_fasta_dir <- function(destdir, region_type)
{
    file.path(destdir, paste0(region_type, "_original_fasta"))
}

.get_final_fasta_path <- function(destdir, region_type)
{
    file.path(destdir, paste0(region_type, ".fasta"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .disambiguate_fasta_seqids()
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
    unlist(suffixes, use.names=FALSE)
    seqids2 <- paste0(seqids2, unlist(suffixes, use.names=FALSE))
    ans <- seqids2[S4Vectors:::reverseIntegerInjection(oo, length(oo))]
    setNames(ans, names(seqids))
}

### In-place replacement!
.disambiguate_fasta_seqids <- function(filepath)
{
    stopifnot(isSingleNonWhiteString(filepath))
    fasta_lines <- readLines(filepath)
    header_idx <- grep("^>", fasta_lines)
    header_lines <- fasta_lines[header_idx]
    if (anyDuplicated(header_lines)) {
        fasta_lines[header_idx] <- .make_unique_seqids(header_lines)
        writeLines(fasta_lines, filepath)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .combine_and_edit_fasta_files()
###

### The workhorse behind create_region_db().
###
### See procedure described at
###   https://ncbi.github.io/igblast/cook/How-to-set-up.html
### for how to create a germline or C-region db from the FASTA files
### available at IMGT. Note that the same procedure can be applied to
### the FASTA files available at AIRR-community/OGRDB.
### This is a 3-step procedure: (1) combine, (2) edit, (3) compile.
### The .combine_and_edit_fasta_files() function below implements
### steps (1) and (2).
### Compilation (with makeblastdb) will happen at a latter time.
.combine_and_edit_fasta_files <- function(fasta_files, destdir,
                                          region_type=c(VDJ_REGION_TYPES, "C"))
{
    if (!is.character(fasta_files) || anyNA(fasta_files))
        stop(wmsg("'fasta_files' must be a character vector with no NAs"))
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    region_type <- match.arg(region_type)

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

    ## (2b) Mangle seq ids to make them unique if they're not.
    .disambiguate_fasta_seqids(final_fasta)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_region_db()
###

.region_db_already_exists <- function(destdir, region_type)
{
    original_fasta_dir <- .get_original_fasta_dir(destdir, region_type)
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
    original_fasta_dir <- .get_original_fasta_dir(destdir, region_type)
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
###   - V_original_fasta/: subdirectory containing the input FASTA files
###         corresponding to the V regions, one FASTA file per region;
###   - V.fasta: the combined and edited FASTA file produced by calling
###         .combine_and_edit_fasta_files() on the files in V_original_fasta/,
###         with allele names disambiguated if needed.
create_region_db <- function(fasta_files, destdir,
                             region_type=c(VDJ_REGION_TYPES, "C"),
                             overwrite=FALSE)
{
    if (!is.character(fasta_files) || anyNA(fasta_files))
        stop(wmsg("'fasta_files' must be a character vector with no NAs"))
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir))
        stop(wmsg("'destdir' must be the path to an existing directory"))
    region_type <- match.arg(region_type)
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))

    if (.region_db_already_exists(destdir, region_type)) {
        if (!overwrite)
            .stop_on_existing_region_db(destdir, region_type)
        .nuke_existing_region_db(destdir, region_type)
    }

    ## Create "original fasta" subdir and copy fasta files to it.
    original_fasta_dir <- .get_original_fasta_dir(destdir, region_type)
    stopifnot(dir.create(original_fasta_dir))
    destfiles <- names(fasta_files)
    if (is.null(destfiles)) {
        stopifnot(all(file.copy(fasta_files, original_fasta_dir)))
    } else {
        destfiles <- file.path(original_fasta_dir, destfiles)
        stopifnot(all(file.copy(fasta_files, destfiles)))
    }

    ## Combine and edit the original fasta files.
    original_files <- list.files(original_fasta_dir, full.names=TRUE)
    .combine_and_edit_fasta_files(original_files, destdir,
                                  region_type=region_type)
}

