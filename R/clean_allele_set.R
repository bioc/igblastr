### =========================================================================
### clean_allele_set() and related
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

check_locus <- function(locus, what)
{
    if (!is.character(locus))
        stop(wmsg(what, " must be of type character"))
    if (anyNA(locus))
        stop(wmsg(what, " cannot contain NAs"))
    if (any(nchar(locus) != 3L))
        stop(wmsg(what, " must contain 3-letter strings"))
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
### .drop_repeated_alleles()
###

### Has its own tests in tests/testthat/test-clean_allele_set.R!
### By "repeated" alleles we mean alleles with identical sequences **and**
### names. Note that:
### - we keep alleles with identical sequences but different names;
### - we also keep alleles with identical names but different sequences.
.drop_repeated_alleles <- function(allele_set, verbose=FALSE)
{
    stopifnot(is(allele_set, "XStringSet") || is.character(allele_set))
    drop_idx <- which(.is_repeated(allele_set))
    drop_idx_len <- length(drop_idx)
    if (drop_idx_len == 0L)
        return(allele_set)  # no-op
    if (verbose) {
        in1string <- paste(names(allele_set)[drop_idx], collapse=", ")
        msg <- c("Dropped ", drop_idx_len, " \"repeated\" allele(s) from ",
                 "set of ", length(allele_set), " alleles: ", in1string, ".")
        note <- c("Note that alleles are considered \"repeated\" when ",
                  "they share the same ungapped sequence and name.")
        message("  o ", wmsg(msg, margin=4L), "\n",
                "    ", wmsg(note, margin=4L), "\n")
    }
    allele_set[-drop_idx]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .disambiguate_allele_names()
###

.stop_on_ambiguous_allele_names <- function(allele_names)
{
    stopifnot(is.character(allele_names))
    is_dup <- duplicated(allele_names)
    ambiguous_allele_names <- unique(allele_names[is_dup])
    in1string <- paste(ambiguous_allele_names, collapse=", ")
    stop(wmsg("The following allele names are ambiguous: ", in1string, "."),
         "\n  ",
         wmsg("Use 'disambiguate.allele.names=TRUE' to disambiguate ",
              "them. (Also using 'verbose=TRUE' will display the ",
              "disambiguation details.)"))
}

.disambiguate_allele_names <- function(allele_set, annotated=FALSE,
                                       verbose=FALSE)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    allele_names <- names(allele_set)
    stopifnot(!is.null(allele_names))
    if (annotated) {
        allele_set_mcols <- mcols(allele_set, use.names=FALSE)
        stopifnot(is(allele_set_mcols, "DataFrame"),
                  identical(allele_names, allele_set_mcols[ , "allele_name"]))
    }
    if (anyDuplicated(allele_names)) {
        new_allele_names <- make_unique_allele_names(allele_names)
        if (verbose) {
            idx <- which(allele_names != new_allele_names)
            in1string <- paste0(allele_names[idx], "->", new_allele_names[idx],
                                collapse=", ")
            msg <- c("Disambiguation: renamed the ", length(idx), " ",
                     "following allele(s): ", in1string)
            message("  o ", wmsg(msg, margin=4L), "\n")
        }
        names(allele_set) <- new_allele_names
        if (annotated) {
            allele_set_mcols[ , "allele_name"] <- names(allele_set)
            mcols(allele_set, use.names=FALSE) <- allele_set_mcols
        }
    }
    allele_set
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .drop_or_rename_alleles()
###

### Has its own tests in tests/testthat/test-clean_allele_set.R!
### Make the "Drop Or Rename" summary, a data.frame with one row per allele
### in 'allele_set', and the following columns:
###   - allele_name: the names on 'allele_set'.
###   - suffix: the disambiguation suffix (possible the empty string)
###     or NA if the allele should be dropped. The disambiguation suffix is
###     the suffix to add to the allele name in order to make allele names
###     unique across the allele set.
###   - locus: [optional] the "locus" metadata column from 'allele_set' if
###     the latter is an XStringSet derivative with such metadata column.
.make_DORsummary <- function(allele_set, disambiguate.allele.names=FALSE)
{
    stopifnot(is(allele_set, "XStringSet") || is.character(allele_set))
    suffix <- rep.int(NA_character_, length(allele_set))
    remaining_idx <- which(!.is_repeated(allele_set))
    remaining_allele_names <- names(allele_set)[remaining_idx]
    if (anyDuplicated(remaining_allele_names) && !disambiguate.allele.names)
        .stop_on_ambiguous_allele_names(remaining_allele_names)
    suffix[remaining_idx] <- make_unique_allele_names(remaining_allele_names,
                                                      suffixes.only=TRUE)
    DORsummary <- data.frame(allele_name=names(allele_set), suffix=suffix)
    if (is(allele_set, "XStringSet"))
        DORsummary$locus <- mcols(allele_set)$locus
    DORsummary
}

### If 'summary.only' is FALSE (the default), returns the modified 'allele_set'
### which is not necessarily parallel to the input 'allele_set' because
### some alleles may have been dropped.
### If 'summary.only' is TRUE, returns a 2-col data.frame parallel
### to 'allele_set'. See function .make_DORsummary() above in this file for the
### details. Also in this case 'verbose' is ignored and operations are quiet.
.drop_or_rename_alleles <- function(allele_set, annotated=FALSE,
                                    disambiguate.allele.names=FALSE,
                                    summary.only=FALSE, verbose=FALSE)
{
    if (summary.only)
        return(.make_DORsummary(allele_set,
                     disambiguate.allele.names=disambiguate.allele.names))

    ## (C) Drop "repeated" alleles.
    cleaned_allele_set <- .drop_repeated_alleles(allele_set, verbose=verbose)

    ## (D) Disambiguate allele names.
    allele_names <- names(cleaned_allele_set)
    if (anyDuplicated(allele_names)) {
        ## Note that, because we did (C), alleles with repeated names
        ## are now guaranteed to have distinct DNA sequences.
        if (!disambiguate.allele.names)
            .stop_on_ambiguous_allele_names(allele_names)
        cleaned_allele_set <- .disambiguate_allele_names(cleaned_allele_set,
                                                         annotated=annotated,
                                                         verbose=verbose)
    }

    cleaned_allele_set
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_V_allele_set() helpers
###

### 'allele_set' must have a "locus" metadata column.
### Returns a DNAStringSet object that:
### - is parallel to input DNAStringSet object 'allele_set';
### - carries the names and metadata columns of input object 'allele_set';
### - also carries the annotations obtained with compute_V_gene_delineations()
###   in additional metadata columns.
.annotate_V_alleles <- function(allele_set)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    allele_names <- names(allele_set)
    stopifnot(!is.null(allele_names))

    allele_set_mcols <- mcols(allele_set, use.names=FALSE)
    if (is.null(allele_set_mcols) || is.null(locus <- allele_set_mcols$locus))
        stop(wmsg("'allele_set' must have a \"locus\" metadata column ",
                  "when 'with.intdata' is TRUE"))
    check_locus(locus, "the \"locus\" metadata column on 'allele_set'")

    ## Annotate the V alleles based on their gaps.
    intdata <- compute_V_gene_delineations(allele_set)
    stopifnot(identical(intdata[ , "allele_name"], allele_names))

    ## Add "chain_type" column.
    intdata$chain_type <- make_chain_type("V", locus)

    mcols(allele_set) <- cbind(allele_set_mcols, intdata)
    allele_set
}

### Check that "repeated" V alleles (i.e. V alleles with identical
### **ungapped** DNA sequences **and** names) are annotated identically.
.stop_if_repeated_V_alleles_have_discordant_annotations <- function(allele_set)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    ## We only care about the annotations found in the Metadata Columns
    ## Of Interest listed in 'MCOI' when comparing the V allele annotations.
    MCOI <- setdiff(names(NDM_DATA_COL2CLASS), "chain_type")
    ok <- .repeated_vector_elts_have_identical_metadata(allele_set, MCOI=MCOI)
    if (!ok) {
        msg <- c("V alleles with identical ungapped sequences and names ",
                 "must have gaps that result in identical annotations (i.e. ",
                 "in identical FWR/CDR boundaries and \"coding_frame_start\")")
        stop(wmsg(msg))
    }
}

.stop_if_allele_sequences_have_gaps <- function(ngaps)
{
    stopifnot(is.integer(ngaps))
    bad_idx <- which(ngaps != 0L)
    if (length(bad_idx) == 0L)
        return()
    allele_names <- names(ngaps)
    stopifnot(!is.null(allele_names))
    first_bad_allele <- allele_names[[bad_idx[[1L]]]]
    stop(wmsg("Some V allele sequences have gaps (", length(bad_idx), " ",
              "out of ", length(ngaps), " input V allele sequences, e.g. ",
              "allele ", first_bad_allele, "). Use 'gapped=TRUE' if the ",
              "supplied V allele sequences are gapped."))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_J_allele_set() helpers
###

### Propagate the names.
.extract_codon_starts <- function(headers)
{
    parsed_headers <- parse_imgt_fasta_headers(headers)
    codon_starts <- parsed_headers[ , "codon_start"]
    allele_name  <- parsed_headers[ , "allele_name"]
    organism     <- parsed_headers[ , "organism"]
    ## The codon start reported by IMGT for human allele TRDJ4*01 has changed
    ## from 3 to 1 between IMGT releases 202603-4 and 202611-4, but the DNA
    ## sequence of the allele has not changed. This looks like a mistake as
    ## it completely messes up the translated sequence, which went from
    ## RPLIFGKGTYLEVQQ (with the FWR4 starting at position 5 where motif FGXG
    ## is found) to PDP*SLAKEPIWRYNN. So we revert this back!
    ## TODO: Report this to the IMGT folks.
    bad_idx <- which(organism == "Homo sapiens" & allele_name == "TRDJ4*01")
    if (length(bad_idx) != 0L) {
        stopifnot(length(bad_idx) == 1L)
        codon_starts[[bad_idx]] <- 3L
    }
    ## The codon start reported by IMGT for mouse allele TRAJ7*01 is 2. This is
    ## suspicious because this produces tranlated sequence DYSNNRLTLGKGTQVVVLP
    ## (no FGXG motif) whereas if we start tranlation at position 1 we get
    ## GLQQQQTYFGEGNPGGGVT (FGXG is now seen at position 9). So we fix that!
    ## TODO: Report this to the IMGT folks.
    bad_idx <- which(organism == "Mus musculus_B10.D2-H2dm1" &
                     allele_name == "TRAJ7*01")
    if (length(bad_idx) != 0L) {
        stopifnot(length(bad_idx) == 1L)
        codon_starts[[bad_idx]] <- 1L
    }
    setNames(as.integer(codon_starts), names(codon_starts))
}

### 'ref_auxdata' (if provided) is used to "complete" the computed
### auxiliary data. Note that the completion process will also verify
### that the computed auxiliary data is concordant with 'ref_auxdata'.
.annotate_J_alleles <- function(allele_set, codon_starts=NULL, ref_auxdata=NULL)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    allele_names <- names(allele_set)
    stopifnot(!is.null(allele_names))

    allele_set_mcols <- mcols(allele_set, use.names=FALSE)
    if (is.null(allele_set_mcols) || is.null(locus <- allele_set_mcols$locus))
        stop(wmsg("'allele_set' must have a \"locus\" metadata column ",
                  "when 'with.auxdata' is TRUE"))
    check_locus(locus, "the \"locus\" metadata column on 'allele_set'")

    ## STEP 1: Compute the auxiliary data by searching motifs WGXG
    ## and FGXG.
    auxdata <- compute_auxdata(allele_set, codon_starts=codon_starts,
                                           no.warnings=TRUE)

    ## Sanity checks.
    stopifnot(identical(auxdata[ , "allele_name"], allele_names))
    expected_chain_type <- make_chain_type("J", locus)
    stopifnot(identical(auxdata[ , "chain_type"], expected_chain_type))

    ## STEP 2: Replace missing values (if any) in 'auxdata' with corresponding
    ## values in 'ref_auxdata'. Also raise an error that displays the
    ## differences if 'auxdata' and 'ref_auxdata' are discordant.
    if (!is.null(ref_auxdata))
        auxdata <- fill_missing_auxdata_with_ref(auxdata, ref_auxdata,
                                                 "computed auxdata",
                                                 "IgBLAST auxdata")

    ## STEP 3: Try to get rid of the remaining missing "cdr3_end" values
    ## in 'auxdata' by comparing the unsolved sequences in 'allele_set' (i.e.
    ## the sequences with a missing "cdr3_end") with the set of known FWR4
    ## sequences.
    auxdata <- infer_cdr3_ends_via_fwr4_comparisons(auxdata, allele_set)

    mcols(allele_set) <- cbind(allele_set_mcols, auxdata)
    allele_set
}

### Check that "repeated" J alleles (i.e. J alleles with identical DNA
### sequences **and** names) are annotated identically.
.stop_if_repeated_J_alleles_have_discordant_annotations <- function(allele_set)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    ## We only care about the annotations found in the Metadata Columns
    ## Of Interest listed in 'MCOI' when comparing the V allele annotations.
    MCOI <- setdiff(names(AUXDATA_COL2CLASS), "chain_type")
    ok <- .repeated_vector_elts_have_identical_metadata(allele_set, MCOI=MCOI)
    if (!ok) {
        msg <- c("J alleles with identical sequences and names must have ",
                 "identical annotations")
        stop(wmsg(msg))
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_allele_set()
### clean_V_allele_set()
### clean_J_allele_set()
###
### These are the workhorses behind create_region_db(), create_V_region_db(),
### and create_J_region_db().
###
### They all perform an **enhanced** version of the cleaning/editing
### implemented in the edit_imgt_file.pl script included in IgBLAST.
### More precisely, they all do the following:
###   (A) clean FASTA headers, as with original edit_imgt_file.pl script;
###   (C) drop "repeated" alleles;
###   (D) disambiguate allele names.
###
### In addition to (A) + (C) + (D) above, clean_V_allele_set() performs the
### following steps:
###   (Bv1) annotate the V alleles if 'gapped' and 'with.intdata' are TRUE;
###   (Bv2) remove gaps from sequences if 'gapped' is TRUE, as with original
###         edit_imgt_file.pl script.
###
### In addition to (A) + (C) + (D) above, clean_J_allele_set() performs the
### following step:
###   (Bj) annotate the J alleles if 'with.auxdata' is TRUE.

### Generic cleaning. Performs: (A) -> (C) -> (D).
### Use it to clean a set of germline D or C gene alleles.
### If, after performing the (C) step (i.e. drop the "repeated" alleles),
### the names of the remaining alleles are not unique, then an error will
### be raised, unless 'disambiguate.allele.names' is set to TRUE, in which
### case the ambiguous allele names will be disambiguated.
### If 'summary.only' is TRUE then 'verbose' is ignored and operations
### are quiet.
clean_allele_set <- function(allele_set,
                             disambiguate.allele.names=FALSE,
                             summary.only=FALSE, verbose=FALSE)
{
    stopifnot(is(allele_set, "XStringSet") || is.character(allele_set))
    headers <- names(allele_set)  # typically the FASTA headers
    stopifnot(!is.null(headers),
              isTRUEorFALSE(disambiguate.allele.names),
              isTRUEorFALSE(verbose))

    ## (A) Clean possibly messy FASTA headers to keep only allele names.
    names(allele_set) <- clean_imgt_fasta_headers(headers)

    ## (C) + (D).
    .drop_or_rename_alleles(allele_set,
                    disambiguate.allele.names=disambiguate.allele.names,
                    summary.only=summary.only, verbose=verbose)
}

### Performs: (A) -> (Bv1) -> (Bv2) -> (C) -> (D).
### Set 'gapped' to TRUE if the sequences in 'allele_set' are gapped (note
### that only V allele sequences are allowed to have gaps).
### Set 'with.intdata' to TRUE if the sequences in 'allele_set' are
### gapped **and** the associated internal data should be computed. Note
### that in this case 'allele_set' **must** carry a "locus" metadata column.
### Returns a DNAStringSet object that contains the ungapped allele sequences
### and carries the metadata columns of input object 'allele_set', if any.
### Furthermore, if 'with.intdata' is TRUE, the returned object will also
### carry the annotations obtained with compute_V_gene_delineations() in
### additional metadata columns.
### If 'summary.only' is TRUE then 'verbose' is ignored and operations
### are quiet.
clean_V_allele_set <- function(allele_set,
                               gapped=FALSE, with.intdata=FALSE,
                               disambiguate.allele.names=FALSE,
                               summary.only=FALSE, verbose=FALSE)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    headers <- names(allele_set)  # typically the FASTA headers
    stopifnot(!is.null(headers),
              isTRUEorFALSE(gapped),
              isTRUEorFALSE(disambiguate.allele.names),
              isTRUEorFALSE(verbose))
    checkarg_with.intdata(with.intdata, gapped)

    ## (A) Clean possibly messy FASTA headers to keep only allele names.
    names(allele_set) <- clean_imgt_fasta_headers(headers)

    ngaps <- setNames(vcountPattern(GAP_LETTER, allele_set), names(allele_set))
    if (gapped) {
        if (with.intdata) {
            ## (Bv1) Annotate the V alleles.
            allele_set <- .annotate_V_alleles(allele_set)
        } else {
            warn_if_allele_sequences_have_no_gaps(ngaps)
        }
        allele_set <- remove_gaps(allele_set)  # (Bv2)
        if (with.intdata)
            .stop_if_repeated_V_alleles_have_discordant_annotations(allele_set)
    } else {
        .stop_if_allele_sequences_have_gaps(ngaps)
    }

    ## (C) + (D).
    ### Note that V alleles are considered "repeated" if the have
    ### identical **ungapped** sequences DNA sequences and names.
    .drop_or_rename_alleles(allele_set, annotated=with.intdata,
                    disambiguate.allele.names=disambiguate.allele.names,
                    summary.only=summary.only, verbose=verbose)
}

### Performs: (A) -> (Bj) -> (C) -> (D).
### If 'summary.only' is TRUE then 'verbose' is ignored and operations
### are quiet.
clean_J_allele_set <- function(allele_set,
                               with.auxdata=FALSE, imgt.fasta=FALSE,
                               ref_auxdata=NULL,
                               disambiguate.allele.names=FALSE,
                               summary.only=FALSE, verbose=FALSE)
{
    stopifnot(is(allele_set, "DNAStringSet"))
    headers <- names(allele_set)  # typically the FASTA headers
    stopifnot(!is.null(headers),
              isTRUEorFALSE(with.auxdata),
              isTRUEorFALSE(imgt.fasta),
              isTRUEorFALSE(disambiguate.allele.names),
              isTRUEorFALSE(verbose))

    ## (A) Clean possibly messy FASTA headers to keep only allele names.
    names(allele_set) <- clean_imgt_fasta_headers(headers)

    if (with.auxdata) {
        if (imgt.fasta) {
            names(headers) <- names(allele_set)
            codon_starts <- .extract_codon_starts(headers)
        } else {
            codon_starts <- NULL
        }
        ## (Bj) Annotate the J alleles.
        allele_set <- .annotate_J_alleles(allele_set,
                                          codon_starts=codon_starts,
                                          ref_auxdata=ref_auxdata)
        .stop_if_repeated_J_alleles_have_discordant_annotations(allele_set)
    }

    ## (C) + (D).
    .drop_or_rename_alleles(allele_set, annotated=with.auxdata,
                    disambiguate.allele.names=disambiguate.allele.names,
                    summary.only=summary.only, verbose=verbose)
}

