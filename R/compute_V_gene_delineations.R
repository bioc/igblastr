### =========================================================================
### Compute IgBLAST internal data
### -------------------------------------------------------------------------
###


### Not exported!
warn_if_allele_sequences_have_no_gaps <- function(ngaps)
{
    stopifnot(is.integer(ngaps))
    bad_idx <- which(ngaps == 0L)
    if (length(bad_idx) == 0L)
        return()
    allele_names <- names(ngaps)
    stopifnot(!is.null(allele_names))
    first_bad_allele <- allele_names[[bad_idx[[1L]]]]
    warning(wmsg(length(bad_idx), "/", length(ngaps), " V allele sequences ",
                 "have no gaps (e.g. allele ", first_bad_allele, ")"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_imgt_fasta_header_lines()
###

### Not exported!
### Performs the same FASTA header cleanup as the edit_imgt_file.pl script
### included in IgBLAST. Note that the latter does some funky business with
### IG allele names for Mus spretus, so we mimick it here.
### Returns the germline gene allele names.
clean_imgt_fasta_header_lines  <- function(headers, what="some allele names")
{
    stopifnot(is.character(headers))
    if (anyNA(headers))
        stop(wmsg(what, " are NAs"))
    if (any(is_white_str(headers)))
        stop(wmsg(what, " are empty"))

    header_parts <- CharacterList(strsplit(headers, "|", fixed=TRUE))

    ## Extract 2nd field. Note that this extraction method could also get
    ## the 1st field if a header has no pipe or if it has only one pipe
    ## with nothing after it.
    allele_names <- tails(heads(header_parts, n=2L), n=1L)
    stopifnot(all(lengths(allele_names) == 1L))
    allele_names <- trimws2(as.character(allele_names))
    if (!all(nchar(allele_names) >= 2L))
        stop(wmsg(what, " are less than 2-character long"))

    ## Implement funky business with Mus spretus: append _Mus_spretus suffix
    ## to IG allele names if species reported in 3rd field is Mus spretus.
    funky_idx <- which(substr(allele_names, 1L, 2L) == "IG" &
                       lengths(header_parts) >= 3L)
    if (length(funky_idx) != 0L) {
        ## Extract 3rd field.
        species <- tails(heads(header_parts[funky_idx], n=3L), n=1L)
        stopifnot(all(lengths(species) == 1L))
        idx <- grep("Mus\\s+spretus", as.character(species))
        if (length(idx) != 0L) {
            funky_idx <- funky_idx[idx]
            allele_names[funky_idx] <-
                paste0(allele_names[funky_idx], "_Mus_spretus")
        }
    }

    allele_names
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compute_V_gene_delineations()
###

### The IMGT unique numbering provides a standardized delimitation of
### the FWR and CDR regions. This standard is based on fixed FWR/CDR lengths
### with respect to the germline V gene **gapped** protein sequences.
### See:
###  https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
### and
###  https://www.imgt.org/IMGTScientificChart/Numbering/IMGT-Kabat_part1.html
.IMGT_FWRCDR_MAX_LENGTHS <- c(
    fwr1=26L,
    cdr1=12L,
    fwr2=17L,
    cdr2=10L,
    fwr3=39L
)

### Do the supplied ranges align with the underlying coding frame?
### Returns a logical vector parallel to 'dna_ranges'.
.dna_ranges_align_with_coding_frame <- function(dna_ranges)
{
    stopifnot(is(dna_ranges, "IRanges"))
    (end(dna_ranges) %% 3L == 0L) & (width(dna_ranges) %% 3L == 0L)
}

### Not used at the moment.
### All the ranges in 'dna_ranges' must align with the underlying coding
### frame that starts at position 1. An error will be raised if they don't.
.from_dna_to_aa_ranges <- function(dna_ranges)
{
    nms <- names(dna_ranges)
    if (is(dna_ranges, "PartitioningByWidth")) {
        widths <- setNames(width(dna_ranges), nms)
        ## Check alignment with underlying coding frame.
        stopifnot(all(widths %% 3L == 0L))
        return(PartitioningByWidth(widths %/% 3L))
    }
    if (is(dna_ranges, "PartitioningByEnd")) {
        ends <- setNames(end(dna_ranges), nms)
        ## Check alignment with underlying coding frame.
        stopifnot(all(ends %% 3L == 0L))
        return(PartitioningByEnd(ends %/% 3L))
    }
    if (!is(dna_ranges, "IntegerRanges"))
        stop(wmsg("'dna_ranges' must be an IRanges object ",
                  "or other IntegerRanges derivative"))
    stopifnot(all(.dna_ranges_align_with_coding_frame(dna_ranges)))
    ans_end <- end(dna_ranges) %/% 3L
    ans_width <- width(dna_ranges) %/% 3L
    IRanges(end=ans_end, width=ans_width, names=nms)
}

### Returns an integer vector parallel to 'imgt_bins' with the following
### attributes:
### - starting_gap: size (in number of nucleotides) of gap block located
###   at the very beginning of the sequence if any (0 if no such block);
### - all_gaps_in_frame: indicates whether the gap blocks align with the
###   underlying coding frame or not;
### - all_gaps_contained: TRUE if the gap blocks don't cross the FWR/CDR
###   boundaries.
.compute_fwrcdr_real_lengths <- function(gap_pos, imgt_bins)
{
    stopifnot(is(gap_pos, "IRanges"), is(imgt_bins, "PartitioningByWidth"))
    gap_pos <- as(gap_pos, "StitchedIPos")

    ## Fastest way to check that the gap positions are strictly sorted.
    gap_blocks <- gap_pos@pos_runs  # IRanges object
    stopifnot(isNormal(gap_blocks))

    ## Note that some gapped V allele sequences from IMGT (e.g. human
    ## IGHV1-69-2*02 and IGHV7-34-1*01) can start with a gap block, and most
    ## of the time this block does not align with the underlying coding frame.
    if (length(gap_blocks) >= 1L && start(gap_blocks)[[1L]] == 1L) {
        starting_gap <- width(gap_blocks)[[1L]]
    } else {
        starting_gap <- 0L
    }

    ## Do all gap blocks align with the underlying coding frame?
    all_gaps_in_frame <- all(.dna_ranges_align_with_coding_frame(gap_blocks))

    ## Check that no gap block crosses FWR/CDR boundaries.
    counts <- countOverlaps(gap_blocks, imgt_bins, type="within")
    all_gaps_contained <- all(counts == 1L)

    imgt_bin_widths <- setNames(width(imgt_bins), names(imgt_bins))
    ngap_per_bin <- count_bin_hits(pos(gap_pos), imgt_bin_widths)
    ans <- imgt_bin_widths - ngap_per_bin

    attr(ans, "starting_gap") <- starting_gap
    attr(ans, "all_gaps_in_frame") <- all_gaps_in_frame
    attr(ans, "all_gaps_contained") <- all_gaps_contained
    ans
}

### Returns an ordinary data.frame.
.IRL_to_data_frame <- function(IRL)
{
    stopifnot(is(IRL, "CompressedIRangesList"))
    IRL_len <- length(IRL)
    all_ranges <- unlist(IRL, use.names=FALSE)
    expected_names <- rep.int(names(.IMGT_FWRCDR_MAX_LENGTHS), IRL_len)
    stopifnot(identical(expected_names, names(all_ranges)))

    idx0 <- seq_len(IRL_len) * length(.IMGT_FWRCDR_MAX_LENGTHS)
    df <- data.frame(
        allele_name=names(IRL),
        fwr1_start =start(all_ranges)[idx0 - 4L],
        fwr1_end   =end  (all_ranges)[idx0 - 4L],
        cdr1_start =start(all_ranges)[idx0 - 3L],
        cdr1_end   =end  (all_ranges)[idx0 - 3L],
        fwr2_start =start(all_ranges)[idx0 - 2L],
        fwr2_end   =end  (all_ranges)[idx0 - 2L],
        cdr2_start =start(all_ranges)[idx0 - 1L],
        cdr2_end   =end  (all_ranges)[idx0 - 1L],
        fwr3_start =start(all_ranges)[idx0],
        fwr3_end   =end  (all_ranges)[idx0]
    )
    cbind(df, mcols(IRL, use.names=FALSE))  # ordinary data.frame
}

### 'gapped_V_alleles' can be a named DNAStringSet or BStringSet object,
### or the path to a FASTA file. Note that **all** the sequences
### in 'gapped_V_alleles' are expected to have gaps (the function will
### issue a warning if that's not the case).
### Returns a data.frame with 1 row per sequence in 'gapped_V_alleles'.
compute_V_gene_delineations <- function(gapped_V_alleles, as.IRangesList=FALSE)
{
    if (!isTRUEorFALSE(as.IRangesList))
        stop(wmsg("'as.IRangesList' must be TRUE or FALSE"))
    if (isSingleString(gapped_V_alleles)) {
        what <- paste0("some allele names in ", gapped_V_alleles)
        ## Some IMGT FASTA files (e.g. for Aotus_nancymaae and
        ## Nonhuman_primates) have nucleotide sequences that contain
        ## the letter 'x'. Not sure what that's supposed to represent.
        ## Note that a well established consensus is to use 'n' or 'N' to
        ## represent an unknown nucleotide (wildcard). Anyways, this breaks
        ## readDNAStringSet() so we use readBStringSet() instead.
        gapped_V_alleles <- readBStringSet(gapped_V_alleles)
        allele_names <- names(gapped_V_alleles)
    } else if (is(gapped_V_alleles, "DNAStringSet") ||
               is(gapped_V_alleles, "BStringSet"))
    {
        allele_names <- names(gapped_V_alleles)
        if (is.null(allele_names))
            stop(wmsg("'gapped_V_alleles' must have names"))
        what <- "some of the names on 'gapped_V_alleles'"
    } else {
        stop(wmsg("'gapped_V_alleles' must be a DNAStringSet or BStringSet ",
                  "object, or the path to a FASTA file containing the gapped ",
                  "sequences of the germline V gene alleles"))
    }

    names(gapped_V_alleles) <- clean_imgt_fasta_header_lines(allele_names, what)

    ## IMGT FWR/CDR fixed intervals in nucleotide space.
    imgt_bins <- PartitioningByWidth(.IMGT_FWRCDR_MAX_LENGTHS * 3L)
    midx <- vmatchPattern(GAP_LETTER, gapped_V_alleles)

    ## Note that lengths() should propagate the names by default but it
    ## fails to do so on ByPos_MIndex object 'midx' at the moment (Biostrings
    ## 2.79.4).
    ## TODO: Fix this in Biostrings. Simplest fix is to define a lengths()
    ## method for MIndex objects that does 'lengths(endIndex(midx))' and
    ## add the names to that if 'use.names' is TRUE.
    ngaps <- setNames(lengths(midx), names(midx))
    warn_if_allele_sequences_have_no_gaps(ngaps)
    seq_len <- width(gapped_V_alleles) - ngaps  # lengths of ungapped sequences

    all_real_lens <- lapply(midx, .compute_fwrcdr_real_lengths, imgt_bins)
    tmp <- lapply(all_real_lens, PartitioningByWidth)
    IRL <- as(tmp, "CompressedIRangesList")
    starting_gap <-
        vapply(all_real_lens, attr, integer(1), "starting_gap")
    all_gaps_in_frame <-
        vapply(all_real_lens, attr, logical(1), "all_gaps_in_frame")
    all_gaps_contained <-
        vapply(all_real_lens, attr, logical(1), "all_gaps_contained")
    coding_frame_start <- 2L - (starting_gap + 2L) %% 3L
    mcols(IRL) <- DataFrame(seq_len=seq_len,
                            coding_frame_start=coding_frame_start,
                            starting_gap=starting_gap,
                            all_gaps_in_frame=all_gaps_in_frame,
                            all_gaps_contained=all_gaps_contained)
    if (as.IRangesList)
        return(IRL)
    .IRL_to_data_frame(IRL)
}

compute_imgt_intdata <- function(...)
{
    .Deprecated("compute_V_gene_delineations")
    compute_V_gene_delineations(...)
}

