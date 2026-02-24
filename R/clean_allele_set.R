### =========================================================================
### clean_allele_set()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


checkarg_with.intdata <- function(with.intdata, gapped)
{
    if (!isTRUEorFALSE(with.intdata))
        stop(wmsg("'with.intdata' must be TRUE or FALSE"))
    if (with.intdata && !isTRUE(gapped))
        stop(wmsg("'with.intdata=TRUE' can only be used when 'gapped' is TRUE"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Detection/handling of "repeated" vector elements
###
### By "repeated" vector elements here we mean elements with identical
### values **and** names. Note that, even though vector-like object 'x' will
### always be a DNAStringSet object in the context of our use cases, we've
### kept the implementation of these low-level utilities agnostic of the
### exact nature of 'x'.

.match_first_repeated <- function(x)
{
    x_names <- names(x)
    stopifnot(!is.null(x_names))
    m1 <- match(x, x)
    m2 <- match(x_names, x_names)
    selfmatchIntegerPairs(m1, m2)
}

### Similar to duplicated() but also takes into account the names on
### vector-like object 'x' to decide whether an element is "repeated" or not.
### Works on any vector-like object that has names and supports match().
.is_repeated <- function(x)
{
    .match_first_repeated(x) != seq_along(x)
}

### Alternate implementation.
### Note that:
### - This is a lot less generic than the implementation above because it
###   relies on splitAsList() working on 'x' and on duplicated() working on
###   the List derivative returned by splitAsList(). Works if 'x' is a
###   DNAStringSet object though, which is all we care about in the context
###   of this file.
### - Not as efficient as .is_repeated() above.
.is_repeated2 <- function(x)
{
    x_names <- names(x)
    stopifnot(!is.null(x_names))
    split_by_name <- splitAsList(x, x_names)
    unsplit(duplicated(split_by_name), x_names)
}

### MCOI: A character vector naming the Metadata Columns Of Interest.
.repeated_vector_elts_have_identical_metadata <- function(x, MCOI=NULL)
{
    stopifnot(is(x, "Vector"))
    x_mcols <- mcols(x, use.names=FALSE)
    stopifnot(is(x_mcols, "DataFrame"))
    if (!is.null(MCOI)) {
        stopifnot(is.character(MCOI))
        x_mcols <- x_mcols[ , MCOI, drop=FALSE]
    }
    m <- .match_first_repeated(x)
    identical(x_mcols[m, , drop=FALSE], x_mcols)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .annotate_V_alleles()
###

.check_locus <- function(locus, what)
{
    if (!is.character(locus))
        stop(wmsg(what, " must be of type character"))
    if (anyNA(locus))
        stop(wmsg(what, " cannot contain NAs"))
    if (any(nchar(locus) != 3L))
        stop(wmsg(what, " must contain 3-letter strings"))
}

### 'dna' must have a "locus" metadata column.
### Returns a DNAStringSet object that:
### - is parallel to input DNAStringSet object 'dna';
### - carries the names and metadata columns of input object 'dna';
### - also carries the annotations obtained with compute_imgt_intdata()
###   in additional metadata columns.
.annotate_V_alleles <- function(dna)
{
    stopifnot(is(dna, "DNAStringSet"))
    dna_mcols <- mcols(dna, use.names=FALSE)
    if (is.null(dna_mcols) || is.null(locus <- dna_mcols$locus))
        stop(wmsg("'dna' must have a \"locus\" metadata column ",
                  "when 'gapped' is TRUE"))
    .check_locus(locus, "the \"locus\" metadata column on 'dna'")

    ## Annotate the V alleles based on their gaps.
    intdata <- compute_imgt_intdata(dna)
    stopifnot(identical(names(dna), intdata[ , "allele_name"]))

    ## Add "chain_type" column.
    intdata$chain_type <- paste0("V", substr(locus, 3L, 3L))

    mcols(dna) <- cbind(dna_mcols, intdata)
    dna
}

### Check that "repeated" V alleles (i.e. V alleles with identical
### **ungapped** DNA sequences **and** names) are annotated identically.
.stop_if_repeated_V_alleles_have_incongruent_gaps <- function(dna)
{
    stopifnot(is(dna, "DNAStringSet"))
    ## We only care about the annotations found in the Metadata Columns
    ## Of Interest listed in 'MCOI' when comparing the V allele annotations.
    MCOI <- setdiff(names(IGBLAST_INTDATA_COL2CLASS), "chain_type")
    ok <- .repeated_vector_elts_have_identical_metadata(dna, MCOI=MCOI)
    if (!ok) {
        msg <- c("V alleles with identical ungapped sequences and names ",
                 "must have gaps that result in identical annotations (i.e. ",
                 "in identical FWR/CDR boundaries and \"coding_frame_start\")")
        stop(wmsg(msg))
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .stop_if_allele_sequences_have_gaps()
###

.stop_if_allele_sequences_have_gaps <- function(ngaps)
{
    stopifnot(is.integer(ngaps))
    bad_idx <- which(ngaps != 0L)
    if (length(bad_idx) == 0L)
        return()
    allele_names <- names(ngaps)
    stopifnot(!is.null(allele_names))
    first_bad_allele <- allele_names[[bad_idx[[1L]]]]
    stop(wmsg("Some allele sequences have gaps (", length(bad_idx), " out ",
              "of ", length(ngaps), " input sequences, e.g. allele ",
              first_bad_allele, "). Use 'gapped=TRUE' if the input ",
              "sequences are gapped V allele sequences."))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .drop_repeated_alleles()
###

### Has its own tests in tests/testthat/test-clean_allele_set.R!
### By "repeated" alleles here we mean alleles with identical sequences
### **and** names. Note that:
### - we keep alleles with identical sequences but different names;
### - we also keep alleles with identical names but different sequences.
.drop_repeated_alleles <- function(dna, verbose=FALSE)
{
    stopifnot(is(dna, "XStringSet") || is.character(dna))
    drop_idx <- which(.is_repeated(dna))
    drop_idx_len <- length(drop_idx)
    if (drop_idx_len == 0L)
        return(dna)  # no-op
    if (verbose) {
        in1string <- paste(names(dna)[drop_idx], collapse=", ")
        msg <- c("Deleted ", drop_idx_len, " \"repeated\" alleles from ",
                 "original set of ", length(dna), " alleles: ", in1string, ".")
        note <- c("Note that alleles are considered \"repeated\" when ",
                  "they share the same ungapped sequence and name.")
        message("  o ", wmsg(msg, margin=4L), "\n",
                "    ", wmsg(note, margin=4L), "\n")
    }
    dna[-drop_idx]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_allele_names_unique()
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

.make_allele_names_unique <- function(dna, with.intdata=FALSE, verbose=FALSE)
{
    stopifnot(is(dna, "DNAStringSet"))
    allele_names <- names(dna)
    stopifnot(!is.null(allele_names))
    if (with.intdata) {
        dna_mcols <- mcols(dna, use.names=FALSE)
        stopifnot(is(dna_mcols, "DataFrame"),
                  identical(allele_names, dna_mcols[ , "allele_name"]))
    }
    if (anyDuplicated(allele_names)) {
        new_allele_names <- .make_unique_seqids(allele_names)
        if (verbose) {
            idx <- which(allele_names != new_allele_names)
            in1string <- paste0(allele_names[idx], "->", new_allele_names[idx],
                                collapse=", ")
            msg <- c("Renamed the ", length(idx), " following ",
                     "allele(s): ", in1string)
            message("  o ", wmsg(msg, margin=4L), "\n")
        }
        names(dna) <- new_allele_names
        if (with.intdata) {
            dna_mcols[ , "allele_name"] <- names(dna)
            mcols(dna, use.names=FALSE) <- dna_mcols
        }
    }
    dna
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_allele_set()
###
### The workhorse behind create_region_db().
### Performs an **enhanced** version of the cleaning/editing implemented
### in the edit_imgt_file.pl script included in IgBLAST. More precisely,
### it does the following:
###   (A) basic cleaning/editing as performed by the original
###       edit_imgt_file.pl script;
###   (B) annotates the V alleles;
###   (C) drops repeated alleles;
###   (D) disambiguates allele names.

### Set 'gapped' to TRUE if the sequences in 'dna' are gapped V allele
### sequences  (note that only V allele sequences are allowed to have gaps).
### Set 'with.intdata' to TRUE if the sequences in 'dna' are gapped V allele
### sequences **and** the associated internal data should be computed. Note
### that in this case 'dna' **must** carry a "locus" metadata column.
### Returns a DNAStringSet object that contains the ungapped allele sequences
### and carries the metadata columns of input object 'dna', if any.
### Furthermore, if 'with.intdata' is TRUE, the returned object will also
### carry the annotations obtained with compute_imgt_intdata() in additional
### metadata columns.
### Note that "repeated" alleles (i.e. alleles with identical **ungapped**
### DNA sequences **and** names) are dropped.
clean_allele_set <- function(dna, gapped=FALSE, with.intdata=FALSE,
                             verbose=FALSE)
{
    stopifnot(is(dna, "DNAStringSet"), isTRUEorFALSE(gapped))
    checkarg_with.intdata(with.intdata, gapped)
    dna_names <- names(dna)
    stopifnot(!is.null(dna_names))
    stopifnot(isTRUEorFALSE(verbose))

    ## (A1) Clean possibly messy FASTA headers to keep only allele names.
    names(dna) <- clean_imgt_fasta_header_lines(dna_names)

    ngaps <- setNames(vcountPattern(GAP_LETTER, dna), names(dna))
    if (gapped) {
        if (with.intdata) {
            dna <- .annotate_V_alleles(dna)  # (B)
        } else {
            warn_if_allele_sequences_have_no_gaps(ngaps)
        }
        dna <- remove_gaps(dna)  # (A2)
        if (with.intdata)
            .stop_if_repeated_V_alleles_have_incongruent_gaps(dna)
    } else {
        .stop_if_allele_sequences_have_gaps(ngaps)
    }

    ## (C) Drop "repeated" alleles.
    dna <- .drop_repeated_alleles(dna, verbose=verbose)

    ## (D) Mangle allele names to make them unique if they're not.
    ##     Note that, because we did (C), remaining repeated allele names
    ##     are guaranteed to be associated with distinct DNA sequences.
    .make_allele_names_unique(dna, with.intdata=with.intdata, verbose=verbose)
}

