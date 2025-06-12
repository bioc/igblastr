### =========================================================================
### augment_germline_db_[VDJ]()
### -------------------------------------------------------------------------


### 'novel_alleles' must be the path to a FASTA file (possibly
### gz-compressed), or a named DNAStringSet object.
### Returns **absolute** path to the **uncompressed** FASTA file
### containing the germline sequences.
.normarg_novel_alleles <- function(novel_alleles)
{
    if (isSingleNonWhiteString(novel_alleles))
        return(fasta_files_as_one_uncompressed_file(novel_alleles,
                                                    "novel_alleles"))
    if (is(novel_alleles, "DNAStringSet")) {
        if (is.null(names(novel_alleles)))
            stop(wmsg("DNAStringSet object 'novel_alleles' ",
                      "must have names"))
        seqlens <- setNames(width(novel_alleles), names(novel_alleles))
        check_seqlens(seqlens, "query")
        path <- tempfile("novel_alleles_", fileext=".fasta")
        writeXStringSet(novel_alleles, path)
        attr(path, "safe_to_remove") <- TRUE
        return(path)
    }
    stop(wmsg("'novel_alleles' must be a single (non-empty) string ",
              "that is the path to a FASTA file, or a named DNAStringSet ",
              "object"))
}

.extract_allele_names <- function(novel_alleles)
{
    stopifnot(isSingleNonWhiteString(novel_alleles))
    allele_names <- trimws2(names(fasta.seqlengths(novel_alleles)))
    if (!all(nzchar(allele_names)))
        stop(wmsg("all the sequences in 'novel_alleles' must have a name"))
    allele_names <- trimws2(sub("\\|.*$", "", allele_names))
    if (!all(nzchar(allele_names)))
        stop(wmsg("all the sequences in 'novel_alleles' must have a name"))
    if (anyDuplicated(allele_names))
        stop(wmsg("some allele names in 'novel_alleles' are duplicated"))
    allele_names
}

.check_novel_allele_names <- function(novel_alleles, db_fasta, db_name)
{
    stopifnot(isSingleNonWhiteString(novel_alleles),
              isSingleNonWhiteString(db_fasta))
    new_allele_names <- .extract_allele_names(novel_alleles)
    cur_allele_names <- names(fasta.seqlengths(db_fasta))
    clashing_names <- intersect(new_allele_names, cur_allele_names)
    if (length(clashing_names) != 0L) {
        in1string <- paste(clashing_names, collapse=", ")
        stop(wmsg("the following allele names in 'novel_alleles' are ",
                  "already present in germline db ", db_name, ": ",
                  in1string))
    }
}

### 'novel_alleles' must be the path to a FASTA file or a named
### DNAStringSet object.
.augment_region_db <- function(db_name, region_type=VDJ_REGION_TYPES,
                               novel_alleles, destdir=".", overwrite=FALSE)
{
    check_germline_db_name(db_name)
    db_path <- make_germline_db_path(db_name)
    region_type <- match.arg(region_type)
    db_fasta <- get_db_fasta_file(db_path, region_type=region_type)
    novel_fasta <- .normarg_novel_alleles(novel_alleles)
    if (isTRUE(attr(novel_fasta, "safe_to_remove")))
        on.exit(unlink(novel_fasta))
    .check_novel_allele_names(novel_fasta, db_fasta, db_name)
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!dir.exists(destdir))
        dir.create(destdir)

    fasta_files <- c(db_fasta, novel_fasta)
    clean_name1 <- paste0(db_name, ":", region_type, ".fasta")
    if (isSingleNonWhiteString(novel_alleles)) {
        ## 'novel_alleles' is the path to a FASTA file (possibly
        ## gz-compressed).
        clean_name2 <- sub("\\.gz$", "", basename(novel_alleles))
    } else {
        ## 'novel_alleles' is a named DNAStringSet object.
        clean_name2 <- basename(novel_fasta)
    }
    names(fasta_files) <- c(clean_name1, clean_name2)
    create_region_db(fasta_files, destdir, region_type=region_type,
                     overwrite=overwrite)
    pattern <- paste0("^", region_type, "\\.fasta$")
    make_blastdbs(destdir, pattern=pattern, force=TRUE)

    message("New augmented ", region_type, " germline db successfully ",
            "created in ", destdir)
    message("To use it with igblastn(), do something like:")
    message("")
    message("    igblastn(..., germline_db_", region_type,
            "=\"", destdir, "\")")
    message("")
    message("See '?augment_germline_db_", region_type, "' for ",
            "more information.")
}

### 'novel_alleles' must be the path to a FASTA file or a named
### DNAStringSet object.
augment_germline_db_V <- function(db_name, novel_alleles, destdir=".",
                                  overwrite=FALSE)
{
    .augment_region_db(db_name, region_type="V", novel_alleles,
                       destdir=destdir, overwrite=overwrite)
}

augment_germline_db_D <- function(db_name, novel_alleles, destdir=".",
                                  overwrite=FALSE)
{
    .augment_region_db(db_name, region_type="D", novel_alleles,
                       destdir=destdir, overwrite=overwrite)
}

augment_germline_db_J <- function(db_name, novel_alleles, destdir=".",
                                  overwrite=FALSE)
{
    .augment_region_db(db_name, region_type="J", novel_alleles,
                       destdir=destdir, overwrite=overwrite)
}

