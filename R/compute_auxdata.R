### =========================================================================
### Compute IgBLAST auxiliary data
### -------------------------------------------------------------------------
###


### For alleles in the IGHJ group (i.e. BCR germline J gene alleles on the
### heavy chain), the FWR4 region is expected to start with AA motif "WGXG.
.WGXG_motif <- "WG.G"
.WGXG_motif_dna <- DNAString("TGGGGNNNNGGN")  # reverse-translation of "WGXG"

### For all other J alleles, that is, for alleles in the IG[KL]J groups
### (i.e. BCR germline J gene alleles on the light chain) and all TCR
### germline J gene alleles, the FWR4 region is expected to start with
### AA motif "FGXG".
.FGXG_motif <- "FG.G"
.FGXG_motif_dna <- DNAString("TTYGGNNNNGGN")  # reverse-translation of "FGXG"

.normarg_codon_starts <- function(codon_starts, allele_names)
{
    if (is.null(codon_starts)) {
        codon_starts <- setNames(rep.int(NA_integer_, length(allele_names)),
                                 allele_names)
        return(codon_starts)
    }
    if (!is.numeric(codon_starts))
        stop(wmsg("'codon_starts' must be NULL or an integer vector"))
    if (is.null(names(codon_starts)))
        stop(wmsg("'codon_starts' must carry the names of the J alleles"))
    if (!identical(names(codon_starts), allele_names))
        stop(wmsg("the names on 'J_alleles' and 'codon_starts' must ",
                  "be identical"))
    if (!is.integer(codon_starts))
        codon_starts <- setNames(as.integer(codon_starts), names(codon_starts))
    if (!(all(codon_starts %in% c(1:3, NA_integer_))))
        stop(wmsg("the non-NA values in 'codon_starts' must ",
                  "be >= 1 and <= 3"))
    codon_starts
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .find_motif()
###

### Returns an integer vector parallel to 'J_alleles' that contains the
### 1-based start position of the motif for each sequence in 'J_alleles'.
### For alleles that don't contain the motif, the position is set to NA.
### For alleles that contain more than one occurence of the motif, we return
### the position of the first occurence.
.find_motif_dna <- function(J_dna, motif_dna)
{
    stopifnot(is(J_dna, "DNAStringSet"), is(motif_dna, "DNAString"))
    m <- vmatchPattern(motif_dna, J_dna, fixed=FALSE)
    as.integer(heads(start(m), n=1L))
}

.find_motif_aa <- function(J_aa, motif)
{
    stopifnot(is(J_aa, "AAStringSet"), isSingleNonWhiteString(motif))
    m <- gregexpr(motif, as.character(J_aa))
    ans <- as.integer(heads(m, n=1L))
    ans[ans == -1L] <- NA_integer_
    ans
}

.find_motif <- function(J_alleles, codon_starts, motif_dna, motif)
{
    stopifnot(is(J_alleles, "DNAStringSet"), is.integer(codon_starts),
              length(J_alleles) == length(codon_starts))
    ans <- rep.int(NA_integer_, length(J_alleles))
    idx1 <- which(is.na(codon_starts))
    if (length(idx1) != 0L)
        ans[idx1] <- .find_motif_dna(J_alleles[idx1], motif_dna)
    idx2 <- which(!is.na(codon_starts))
    if (length(idx2) != 0L) {
        codon_starts2 <- codon_starts[idx2]
        J_aa <- translate_codons(J_alleles[idx2], offset=codon_starts2 - 1L)
        starts2 <- .find_motif_aa(J_aa, motif)
        ans[idx2] <- codon_starts2 + 3L * (starts2 - 1L)
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compute_auxdata()
###

### Not exported!
J_alleles_with_missing_coding_frame_start <- function(auxdata)
{
    stopifnot(is.data.frame(auxdata))
    bad_idx <- which(is.na(auxdata[ , "coding_frame_start"]))
    auxdata[bad_idx, "allele_name"]
}

### Warn user if "coding frame start" could not be determined for some alleles.
.warn_if_missing_coding_frame_starts <- function(auxdata)
{
    bad_alleles <- J_alleles_with_missing_coding_frame_start(auxdata)
    if (length(bad_alleles) != 0L) {
        in1string <- paste(bad_alleles, collapse=", ")
        warning(wmsg("The \"coding frame start\" could not be determined ",
                     "for J allele(s): ", in1string, "."),
                "\n  ",
                wmsg("--> coding_frame_start, cdr3_end, and extra_bps ",
                     "were set to NA for these alleles"))
    }
}

### Not exported!
### Very rare but can happen (we've only seen this for rhesus monkey
### allele IGHJ0-ZXTW*01 from OGRDB so far).
warn_if_negative_cdr3_end <- function(auxdata, what)
{
    bad_idx <- which(auxdata[ , "cdr3_end"] < 0L)
    if (length(bad_idx) != 0L) {
        in1string <- paste(auxdata[bad_idx , "allele_name"], collapse=", ")
        warning(wmsg("The ", what, " is negative for J allele(s): ",
                     in1string, "."),
                "\n  ",
                wmsg("Note that you won't be able to save the returned ",
                     "data.frame with write_auxdata() unless you replace ",
                     "the negative \"cdr3_end\" value(s) with NA(s) or drop ",
                     "the row(s) with a negative \"cdr3_end\" value."))
    }
}

### Compute the auxiliary data for 'J_alleles' by searching motifs WGXG
### and FGXG.
### Returns a data.frame with 1 row per sequence in 'J_alleles' and with
### the same columns as the data.frame returned by load_auxdata() (see
### file R/auxdata-utils.R).
compute_auxdata <- function(J_alleles, codon_starts=NULL, no.warnings=FALSE)
{
    if (!is(J_alleles, "DNAStringSet"))
        stop(wmsg("'J_alleles' must be DNAStringSet object"))
    allele_names <- names(J_alleles)
    if (is.null(allele_names))
        stop(wmsg("'J_alleles' must have names"))

    names(J_alleles) <- allele_names <- clean_imgt_fasta_headers(allele_names)

    allele_loci <- substr(allele_names, 1L, 3L)
    loci_prefix <- extract_loci_prefix(allele_loci)
    if (is.na(loci_prefix))
        stop(wmsg("all allele names must start either ",
                  "with 'IG[HKL]' or with 'TR[ABGD]'"))
    if (!all(substr(allele_names, 4L, 4L) == "J"))
        stop(wmsg("the 4th letter in all allele names must be a J"))

    codon_starts <- .normarg_codon_starts(codon_starts, allele_names)

    if (!isTRUEorFALSE(no.warnings))
        stop(wmsg("'no.warnings' must be TRUE or FALSE"))

    ## Compute 'cdr3_ends'.
    ## We use WGXG motif for J alleles from IGH locus, and FGXG motif for
    ## J alleles from all other loci.
    fwr4_starts <- rep.int(NA_integer_, length(J_alleles))
    idx1 <- which(allele_loci == "IGH")
    fwr4_starts[idx1] <- .find_motif(J_alleles[idx1], codon_starts[idx1],
                                    .WGXG_motif_dna, .WGXG_motif)
    idx2 <- which(allele_loci != "IGH")
    fwr4_starts[idx2] <- .find_motif(J_alleles[idx2], codon_starts[idx2],
                                    .FGXG_motif_dna, .FGXG_motif)
    cdr3_ends <- fwr4_starts - 1L  # 1-based

    ## Compute 'coding_frame_starts'.
    coding_frame_starts <- unname(codon_starts) - 1L
    stopifnot(all(coding_frame_starts == cdr3_ends %% 3L, na.rm=TRUE))
    idx <- which(is.na(coding_frame_starts) & !is.na(cdr3_ends))
    coding_frame_starts[idx] <- cdr3_ends[idx] %% 3L

    ## Compute 'chain_types' and 'extra_bps'.
    chain_types <- make_chain_type("J", allele_loci)
    extra_bps <- (width(J_alleles) - coding_frame_starts) %% 3L

    ans <- data.frame(
        allele_name       =allele_names,
        coding_frame_start=coding_frame_starts,  # 0-based
        chain_type        =chain_types,
        cdr3_end          =cdr3_ends - 1L,  # 0-based
        extra_bps         =extra_bps
    )

    if (!no.warnings) {
        .warn_if_missing_coding_frame_starts(ans)
        what <- "computed CDR3 end position"
        warn_if_negative_cdr3_end(ans, what)
    }

    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### find_discordant_auxdata()
### fill_missing_auxdata_with_ref()
###

### Not exported!
### See find_discordant_rows() in R/utils.R for what it returns.
find_discordant_auxdata <- function(auxdata, ref_auxdata)
{
    check_auxdata_col2class(auxdata)
    check_auxdata_col2class(ref_auxdata)
    find_discordant_rows(auxdata, ref_auxdata, "allele_name")
}

### Not exported!
### See complete_df_with_refdf() in R/utils.R for what it returns.
fill_missing_auxdata_with_ref <- function(auxdata, ref_auxdata,
                                          auxdata_name="auxdata",
                                          ref_auxdata_name="ref_auxdata")
{
    check_auxdata_col2class(auxdata)
    check_auxdata_col2class(ref_auxdata)
    complete_df_with_refdf(auxdata, ref_auxdata, "allele_name",
                           auxdata_name, ref_auxdata_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### infer_cdr3_ends_via_fwr4_comparisons()
###

.normarg_maxdist <- function(maxdist)
{
    if (!isSingleNumber(maxdist))
        stop(wmsg("'maxdist' must be a single integer"))
    if (!is.integer(maxdist))
        maxdist <- as.integer(maxdist)
    if (!(maxdist >= 0L && maxdist <= 5L))
        stop(wmsg("'maxdist' must be >= 0 and <= 5"))
    maxdist
}

.normarg_standout.by <- function(standout.by)
{
    if (!isSingleNumber(standout.by))
        stop(wmsg("'standout.by' must be a single integer"))
    if (!is.integer(standout.by))
        standout.by <- as.integer(standout.by)
    if (!(standout.by >= 1L && standout.by <= 10L))
        stop(wmsg("'standout.by' must be >= 1 and <= 10"))
    standout.by
}

### Most FWR4 we've seen so far have at least 10 amino acids except in
### mouse J allele TRAJ45*02 from IMGT where it has only 9 amino acids
### (this allele sequence is probably truncated).
.MIN_FWR4_LENGTH <- 9L

### Moves a sliding window of 10 amino acids along 'J_allele', and, for
### each position of the window, computes the Hamming distance between the
### sequence in the window and the set of known FWR4 sequences in 'fwr4set'.
### Returns the distances in an integer vector of length 'nchar(J_allele)' - 9.
.compute_sliding_window_dist_to_fwr4set <- function(J_allele, fwr4set)
{
    stopifnot(is(J_allele, "AAString"), nchar(J_allele) >= 12L,
              is(fwr4set, "AAStringSet"),
              all(nchar(fwr4set) >= .MIN_FWR4_LENGTH))
    vapply(seq_len(nchar(J_allele) - 9L),
        function(pos)
            min(neditAt(subseq(J_allele, start=pos, width=10L), fwr4set)),
        integer(1)
    )
}

### Has its own tests in tests/testthat/test-compute_auxdata.R!
### The "best FWR4 start position" in 'J_allele' is defined as follow:
### If D(p) is the distance between the 10-amino acid window in 'J_allele'
### that starts at position 'p' and the set of known FWR4 sequences
### in 'fwr4set', then the "best FWR4 start position" is the position 'P'
### in 'J_allele' that minimize D(P). With the additional constraints that:
###   (a) D(P) <= 2
###   (b) D(p) >= D(P) + 5 for any other position 'p'
### In other words, P must stand out in the crowd!
.find_best_fwr4_start <- function(J_allele, fwr4set,
                                  maxdist=2L, standout.by=5L)
{
    maxdist <- .normarg_maxdist(maxdist)
    standout.by <- .normarg_standout.by(standout.by)
    if (nchar(J_allele) < 12L)
        return(NA_integer_)
    dists <- .compute_sliding_window_dist_to_fwr4set(J_allele, fwr4set)
    P <- which.min(dists)
    DP <- dists[[P]]
    if (DP > maxdist)
        return(NA_integer_)  # no P achieves D(P) <= maxdist
    if (sum(dists < (DP + standout.by)) != 1L)
        return(NA_integer_)  # P does not stand out in the crowd
    P
}

### Returns an integer vector parallel to 'J_alleles' that contains the
### 1-based FWR4 start positions.
.find_best_fwr4_starts <- function(J_alleles, fwr4set,
                                   maxdist=2L, standout.by=5L)
{
    stopifnot(is(J_alleles, "AAStringSet"),
              is(fwr4set, "AAStringSet"),
              all(nchar(fwr4set) >= .MIN_FWR4_LENGTH))
    vapply(seq_along(J_alleles),
        function(i)
            .find_best_fwr4_start(J_alleles[[i]], fwr4set,
                                  maxdist=maxdist, standout.by=standout.by),
        integer(1)
    )
}

infer_cdr3_ends_via_fwr4_comparisons <-
    function(auxdata, J_alleles, maxdist=2L, standout.by=5L)
{
    check_auxdata_col2class(auxdata)
    if (!is(J_alleles, "DNAStringSet") || length(J_alleles) != nrow(auxdata))
        stop(wmsg("'J_alleles' must be a DNAStringSet object ",
                  "with one sequence per row in 'auxdata'"))
    if (!identical(names(J_alleles), auxdata[ , "allele_name"]))
        stop(wmsg("'J_alleles' must have names and they must ",
                  "be identical to 'auxdata$allele_name'"))
    maxdist <- .normarg_maxdist(maxdist)
    standout.by <- .normarg_standout.by(standout.by)

    coding_frame_start <- auxdata[ , "coding_frame_start"]
    if (anyNA(coding_frame_start))
        stop(wmsg("'auxdata$coding_frame_start' cannot contain NAs"))

    cdr3_end <- auxdata[ , "cdr3_end"]  # 0-based
    solve_me <- is.na(cdr3_end)
    if (!any(solve_me))
        return(auxdata)

    cdr3_end <- cdr3_end + 1L  # 1-based
    fwr4set <- translate_codons(J_alleles[!solve_me],
                                offset=cdr3_end[!solve_me])
    unsolved_J_aa <- translate_codons(J_alleles[solve_me],
                                offset=coding_frame_start[solve_me])
    solved_J_aa_fwr4_start <-
        .find_best_fwr4_starts(unsolved_J_aa, fwr4set,
                               maxdist=maxdist, standout.by=standout.by)
    solved_J_aa_cdr3_end <- solved_J_aa_fwr4_start - 1L
    solved_cdr3_end <- coding_frame_start[solve_me] + 3L * solved_J_aa_cdr3_end
    solved_cdr3_end <- solved_cdr3_end - 1L  # 0-based
    auxdata[solve_me, "cdr3_end"] <- solved_cdr3_end
    warn_if_negative_cdr3_end(auxdata, "CDR3 end position")
    auxdata
}

