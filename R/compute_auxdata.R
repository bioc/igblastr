### =========================================================================
### Compute IgBLAST auxiliary data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .find_heavy_fwr4_starts()
### .find_light_fwr4_starts()
###

### For alleles in the IGHJ group (i.e. BCR germline J gene alleles on the
### heavy chain), the FWR4 region is expected to start with AA motif "WGXG.
.WGXG_pattern <- "TGGGGNNNNGGN"  # reverse-translation of "WGXG"

### For all other J alleles, that is, for alleles in the IG[KL]J groups
### (i.e. BCR germline J gene alleles on the light chain) and all TCR
### germline J gene alleles, the FWR4 region is expected to start with
### AA motif "FGXG".
.FGXG_pattern <- "TTYGGNNNNGGN"  # reverse-translation of "FGXG"

### EXPERIMENTAL!
### The "FGXG" motif is not found for 4 J alleles in
### IMGT-202531-1.Mus_musculus.IGH+IGK+IGL: IGKJ3*01, IGKJ3*02,
### IGLJ2P*01, IGLJ3P*01. However, except for IGLJ2P*01, these alleles
### are annotated in mouse_gl.aux with a CDR3 end reported at position 6
### (0-based). Turns out that for the 3 alleles annotated in mouse_gl.aux,
### the two first codons of the FWR4 region translate to AA sequence "FS".
### Is this a coincidence or does the FS sequence actually play a role on
### the light chain? What do biologists say about this? In particular, does
### it make sense to use this alternative motif to identify the start of
### the FWR4 region on the light chain when the "FGXG" motif is not found?
### Note that all the possible reverse-translations of FS cannot be
### represented with a single DNA pattern (even with the use of IUPAC
### ambiguity codes).
.FS_pattern1 <- "TTYTCN"
.FS_pattern2 <- "TTYAGY"

### UPDATE on using the "FS" motif to identify the start of the FWR4
### region on the light chain when the "FGXG" motif is not found:
### Works well for IMGT-202531-1.Mus_musculus.IGH+IGK+IGL (well, it was
### specifically designed for that so no surprise here), but not
### so well for IMGT-202531-1.Rattus_norvegicus.IGH+IGK+IGL or
### IMGT-202531-1.Oryctolagus_cuniculus.IGH+IGK+IGL (rabbit)
### or IMGT-202531-1.Macaca_mulatta.IGH+IGK+IGL (rhesus monkey).
### So we disabled this feature in .find_light_fwr4_starts() below.

### .find_heavy_fwr4_starts() and .find_light_fwr4_starts() both return
### a named integer vector parallel to 'J_alleles' that contains
### the **0-based** FWR4 start position for each sequence in 'J_alleles'.
### Th FWR4 start will be set to NA for alleles that don't have a match.
### For alleles with more than one match, we keep the first match only.
### The names on the returned vector indicate the AA motif that was used
### to determine the start of the FWR4 region.

.find_heavy_fwr4_starts <- function(J_alleles)
{
    stopifnot(is(J_alleles, "DNAStringSet"))
    m <- vmatchPattern(.WGXG_pattern, J_alleles, fixed=FALSE)
    ans <- as.integer(heads(start(m), n=1L)) - 1L
    names(ans) <- ifelse(is.na(ans), NA_character_, "WGXG")
    ans
}

.find_light_fwr4_starts <- function(J_alleles)
{
    stopifnot(is(J_alleles, "DNAStringSet"))
    m <- vmatchPattern(.FGXG_pattern, J_alleles, fixed=FALSE)
    FGXG_starts <- as.integer(heads(start(m), n=1L))
    names(FGXG_starts) <- ifelse(is.na(FGXG_starts), NA_character_, "FGXG")
    ## Disabling search for alternative "FS" motif for now.
    #na_idx <- which(is.na(FGXG_starts))
    #if (length(na_idx) != 0L) {
    #    dangling_alleles <- J_alleles[na_idx]
    #    m <- vmatchPattern(.FS_pattern1, dangling_alleles, fixed=FALSE)
    #    FS_starts1 <- as.integer(heads(start(m), n=1L))
    #    m <- vmatchPattern(.FS_pattern2, dangling_alleles, fixed=FALSE)
    #    FS_starts2 <- as.integer(heads(start(m), n=1L))
    #    FS_starts <- pmin(FS_starts1, FS_starts2, na.rm=TRUE)
    #    FGXG_starts[na_idx] <- FS_starts
    #    names(FGXG_starts)[na_idx] <-
    #        ifelse(is.na(FS_starts), NA_character_, "FS")
    #}
    FGXG_starts - 1L
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compute_auxdata()
###

### Returns a data.frame with the same column names as the data.frame
### returned by load_auxdata() (see file R/auxdata-utils.R), minus
### the "extra_bps" column and plus the "fwr4_start_motif" column.
### NOTE: We set coding_frame_start/cdr3_end/fwr4_start_motif
### to NA for alleles for which the FWR4 start cannot be determined.
.compute_auxdata_for_locus <- function(J_alleles, locus)
{
    stopifnot(is(J_alleles, "DNAStringSet"), locus %in% c(IG_LOCI , TR_LOCI))
    allele_names <- names(J_alleles)
    stopifnot(!is.null(allele_names))
    if (length(J_alleles) == 0L) {
        chain_type <- character(0)
    } else {
        allele_loci <- substr(allele_names, 1L, 3L)
        stopifnot(all(allele_loci == locus))
        chain_type <- make_chain_type("J", locus)
    }
    if (locus == "IGH") {
        fwr4_starts <- .find_heavy_fwr4_starts(J_alleles)
    } else {
        ## Used for loci IGK and IGL (light chain) and all TR* loci.
        fwr4_starts <- .find_light_fwr4_starts(J_alleles)
    }
    coding_frame_starts <- unname(fwr4_starts) %% 3L
    data.frame(
        allele_name       =allele_names,
        coding_frame_start=coding_frame_starts,      # 0-based
        chain_type        =chain_type,
        cdr3_end          =unname(fwr4_starts) - 1L  # 0-based
        ## Returning this column only made sense when we were using "FS"
        ## motif as a 2nd-chance motif on the light chain.
        #fwr4_start_motif  =names(fwr4_starts)
    )
}

.normarg_codon_starts <- function(codon_starts, allele_names)
{
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

.refine_with_codon_starts <- function(auxdata, codon_starts)
{
    stopifnot(is.data.frame(auxdata), is.integer(codon_starts),
              identical(auxdata[ , "allele_name"], names(codon_starts)))
    coding_frame_start <- auxdata[ , "coding_frame_start"]
    expected_coding_frame_start <- codon_starts - 1L
    bad_idx <- which(coding_frame_start != expected_coding_frame_start)
    if (length(bad_idx) != 0L) {
        in1string <- paste(auxdata[bad_idx, "allele_name"], collapse=", ")
        stop(wmsg("the supplied \"codon start\" is in disagreement with ",
                  "the computed \"coding frame start\" for allele(s): ",
                  in1string))
    }
    na_idx <- which(is.na(coding_frame_start))
    if (length(na_idx) != 0L) {
        coding_frame_start[na_idx] <- expected_coding_frame_start[na_idx]
        auxdata$coding_frame_start <- coding_frame_start
    }
    auxdata
}

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
### Returns a data.frame with 1 row per sequence in 'J_alleles'.
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

    if (!is.null(codon_starts))
        codon_starts <- .normarg_codon_starts(codon_starts, allele_names)

    if (!isTRUEorFALSE(no.warnings))
        stop(wmsg("'no.warnings' must be TRUE or FALSE"))

    if (loci_prefix == "IG") {
        JH_alleles <- J_alleles[allele_loci == "IGH"]
        JK_alleles <- J_alleles[allele_loci == "IGK"]
        JL_alleles <- J_alleles[allele_loci == "IGL"]
        JH_df <- .compute_auxdata_for_locus(JH_alleles, "IGH")
        JK_df <- .compute_auxdata_for_locus(JK_alleles, "IGK")
        JL_df <- .compute_auxdata_for_locus(JL_alleles, "IGL")
        ans <- rbind(JH_df, JK_df, JL_df)
    } else {
        JA_alleles <- J_alleles[allele_loci == "TRA"]
        JB_alleles <- J_alleles[allele_loci == "TRB"]
        JG_alleles <- J_alleles[allele_loci == "TRG"]
        JD_alleles <- J_alleles[allele_loci == "TRD"]
        JA_df <- .compute_auxdata_for_locus(JA_alleles, "TRA")
        JB_df <- .compute_auxdata_for_locus(JB_alleles, "TRB")
        JG_df <- .compute_auxdata_for_locus(JG_alleles, "TRG")
        JD_df <- .compute_auxdata_for_locus(JD_alleles, "TRD")
        ans <- rbind(JA_df, JB_df, JG_df, JD_df)
    }
    rownames(ans) <- NULL

    i <- match(allele_names, ans[ , "allele_name"])
    ans <- S4Vectors:::extract_data_frame_rows(ans, i)

    if (!is.null(codon_starts))
        ans <- .refine_with_codon_starts(ans, codon_starts)

    ## Add "extra_bps" column.
    ans$extra_bps <- (width(J_alleles) - ans[ , "coding_frame_start"]) %% 3L

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

### Moves a sliding window of 10 amino acids across 'J_allele', and, for
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### infer_cdr3_ends_via_full_J_sequence_comparisons()
###
### Not used and superseded by infer_cdr3_ends_via_fwr4_comparisons() above!
###
### For each J allele with an unknown CDR3 end, we infer it from the J
### alleles with a known CDR3 end that are "close". Proximity is defined by
### the Hamming distance between the amino acid sequences.
### Based on Biostrings::matchPattern().

.find_cdr3_end <- function(x, y, y_cdr3_end, max.mismatch)
{
    stopifnot(is(x, "AAString"),
              is(y, "AAStringSet"),
              is.integer(y_cdr3_end),
              length(y) == length(y_cdr3_end),
              identical(names(y), names(y_cdr3_end)),
              isSingleNumber(max.mismatch))
    delta <- rep.int(NA_integer_, length(y))
    for (j in seq_along(y)) {
        subject <- y[[j]]
        if (length(x) <= length(subject)) {
            v <- matchPattern(x, subject, max.mismatch=max.mismatch)
            d <- start(v) - 1L
        } else {
            v <- matchPattern(subject, x, max.mismatch=max.mismatch)
            d <- 1L - start(v)
        }
        if (length(d) >= 2L)
            return(NA_integer_)
        if (length(d) == 0L)
            next
        delta[[j]] <- d
    }
    x_cdr3_end <- y_cdr3_end - delta
    ans <- unique(x_cdr3_end[!is.na(x_cdr3_end)])
    if (length(ans) == 1L) ans else NA_integer_
}

### Returns an integer vector parallel to 'unsolved_J_aa' that contains the
### 1-based CDR3 end positions w.r.t. to the sequences in 'unsolved_J_aa'.
.find_cdr3_ends <- function(unsolved_J_aa,
                            ref_J_aa, ref_J_aa_cdr3_end, max.mismatch=2L)
{
    stopifnot(is(unsolved_J_aa, "AAStringSet"),
              is(ref_J_aa, "AAStringSet"),
              is.integer(ref_J_aa_cdr3_end),
              length(ref_J_aa) == length(ref_J_aa_cdr3_end),
              identical(names(ref_J_aa), names(ref_J_aa_cdr3_end)),
              isSingleNumber(max.mismatch))
    vapply(seq_along(unsolved_J_aa),
        function(i)
            .find_cdr3_end(unsolved_J_aa[[i]],
                           ref_J_aa, ref_J_aa_cdr3_end,
                           max.mismatch=max.mismatch),
        integer(1)
    )
}

### Not exported!
infer_cdr3_ends_via_full_J_sequence_comparisons <-
    function(auxdata, J_alleles, max.mismatch=2L)
{
    check_auxdata_col2class(auxdata)
    coding_frame_start <- auxdata[ , "coding_frame_start"]
    stopifnot(!anyNA(coding_frame_start),
              is(J_alleles, "DNAStringSet"),
              identical(names(J_alleles), auxdata[ , "allele_name"]),
              isSingleNumber(max.mismatch))

    cdr3_end <- auxdata[ , "cdr3_end"]  # 0-based
    solve_me <- is.na(cdr3_end)
    if (!any(solve_me))
        return(auxdata)

    cdr3_end <- setNames(cdr3_end + 1L, names(J_alleles))  # 1-based
    cdr3_coding_frame_width <- cdr3_end - coding_frame_start
    stopifnot(all(cdr3_coding_frame_width %% 3L == 0L, na.rm=TRUE))
    J_aa_cdr3_end <- cdr3_coding_frame_width %/% 3L

    J_aa <- translate_codons(J_alleles, offset=coding_frame_start)
    unsolved_J_aa <- J_aa[solve_me]
    ref_J_aa <- J_aa[!solve_me]
    ref_J_aa_cdr3_end <- J_aa_cdr3_end[!solve_me]
    solved_J_aa_cdr3_end <- .find_cdr3_ends(unsolved_J_aa,
                                            ref_J_aa, ref_J_aa_cdr3_end,
                                            max.mismatch=max.mismatch)
    solved_cdr3_end <- coding_frame_start[solve_me] + 3L * solved_J_aa_cdr3_end
    solved_cdr3_end <- solved_cdr3_end - 1L  # 0-based
    auxdata[solve_me, "cdr3_end"] <- solved_cdr3_end
    warn_if_negative_cdr3_end(auxdata, "CDR3 end position")
    auxdata
}

