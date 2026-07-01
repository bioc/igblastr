### =========================================================================
### percent_mutation()
### -------------------------------------------------------------------------


.stop_on_missing_AIRR_cols <- function(missing_colnames)
{
    in1string <- paste(paste0("\"", missing_colnames, "\""), collapse=", ")
    what <- if (length(missing_colnames) == 1L) " is" else "s are"
    msg <- c("The following column", what, " missing in the supplied ",
             "data.frame or tibble: ", in1string, ". Was it obtained ",
             "with igblastn()?")
    stop(wmsg(msg))
}

### In some rare situations, the sequences in 'v_sequence_alignment_aa'
### don't have the same lengths as the corresponding sequences
### in 'v_germline_alignment_aa'. This can happen for example when, for
### some query, 'v_sequence_alignment' has a 1-nucleotide insertion or
### deletion w.r.t. 'v_germline_alignment'. For example:
###
###   v_sequence_alignment GATTCTGCAACTTATTACTGCCAACAGTCTAATA-TTATTCT
###                        |||| ||||||||||||||||||||||| ||||| |||||||
###   v_germline_alignment GATTTTGCAACTTATTACTGCCAACAGTATAATAGTTATTCT
###
### Because of this 1-nucleotide deletion, 'v_sequence_alignment' has one
### less codon than 'v_germline_alignment' (13 codons vs 14 codons). So the
### sequences translate to:
###
###   v_sequence_alignment_aa DSATYYCQQSNII
###                           | ||||||| |
###   v_germline_alignment_aa DFATYYCQQYNSYS
###
### with the consequence that 'v_sequence_alignment_aa'
### and 'v_germline_alignment_aa' don't have the same length.
### Note that this deletion messes up the coding frame, with the consequence
### that the downstream codons on both nucledotide sequences translate to
### completely different amino acid sequences.

### Note that .compute_hamming_dist() is 100x faster or more than
### .compute_hamming_dist_OLD().
.compute_hamming_dist <- function(sequence_aln, germline_aln, for.aa=FALSE)
{
    stopifnot(is.character(sequence_aln), !anyNA(sequence_aln),
              is.character(germline_aln), !anyNA(germline_aln),
              length(sequence_aln) == length(germline_aln),
              isTRUEorFALSE(for.aa))
    nc1 <- nchar(sequence_aln)
    nc2 <- nchar(germline_aln)
    if (for.aa) {
        ## When 'for.aa' is TRUE, 'nc1' and 'nc2' are not guaranteed to
        ## be identical. More precisely, the lengths of some sequences
        ## in 'sequence_aln' and their corresponding sequences in 'germline_aln'
        ## can differ by 1. See comment above why.
        ## To make them the same lengths, we truncate the longest to the
        ## lengths of the shortest.
        nc <- pmin(nc1, nc2)
        sequence_aln <- substr(sequence_aln, 1L, nc)
        germline_aln <- substr(germline_aln, 1L, nc)
    } else {
        stopifnot(identical(nc1, nc2))
    }
    vapply(seq_along(sequence_aln),
        function(i) {
            r1 <- charToRaw(sequence_aln[[i]])
            r2 <- charToRaw(germline_aln[[i]])
            sum(r1 != r2)
        },
        integer(1))
}

### No longer used. Superseded by much faster .compute_hamming_dist() above.
.compute_hamming_dist_OLD <- function(sequence_aln, germline_aln, for.aa=FALSE)
{
    stopifnot(is.character(sequence_aln), !anyNA(sequence_aln),
              is.character(germline_aln), !anyNA(germline_aln),
              length(sequence_aln) == length(germline_aln),
              isTRUEorFALSE(for.aa))
    if (for.aa) {
        sequence_aln <- AAStringSet(sequence_aln)
        germline_aln <- AAStringSet(germline_aln)
    } else {
        ## We only do this sanity check when 'for.aa' is FALSE. See comment
        ## above why.
        stopifnot(identical(nchar(sequence_aln), nchar(germline_aln)))
        sequence_aln <- DNAStringSet(sequence_aln)
        germline_aln <- DNAStringSet(germline_aln)
    }
    vapply(seq_along(sequence_aln),
        function(i) neditAt(sequence_aln[[i]], germline_aln[[i]]),
        integer(1))
}

.region_percent_mutation <- function(AIRR_df, region_type=VDJC_REGION_TYPES,
                                     for.aa=FALSE)
{
    stopifnot(is.data.frame(AIRR_df),
              isSingleNonWhiteString(region_type),
              isTRUEorFALSE(for.aa))
    region_type <- match.arg(region_type)

    ## Extract the Columns Of Interest.
    COI <- paste0(tolower(region_type),
                  c("_sequence_alignment", "_germline_alignment"))
    if (for.aa)
        COI <- paste0(COI, "_aa")
    missing_idx <- which(COI %notin% colnames(AIRR_df))
    if (length(missing_idx) != 0L)
        .stop_on_missing_AIRR_cols(COI[missing_idx])
    sequence_aln <- AIRR_df[[COI[[1L]]]]
    germline_aln <- AIRR_df[[COI[[2L]]]]
    notna_idx <- which(!is.na(sequence_aln))
    stopifnot(identical(notna_idx, which(!is.na(germline_aln))))

    ## If a Column Of Interest is filled with NAs only (can happen for the
    ## d_* columns if all the query sequences are from the light chain),
    ## then there's no guarantee that is will be of type character (could
    ## be of type logical or integer).
    if (length(notna_idx) == 0L) {
        sequence_aln <- as.character(sequence_aln)
        germline_aln <- as.character(germline_aln)
    }

    d <- rep.int(NA_integer_, nrow(AIRR_df))
    d[notna_idx] <- .compute_hamming_dist(sequence_aln[notna_idx],
                                          germline_aln[notna_idx],
                                          for.aa=for.aa)
    ## See .compute_hamming_dist() for why we do pmin() here.
    nc <- pmin(nchar(sequence_aln), nchar(germline_aln))
    100 * d / nc
}

### Computes percent mutation in V, D, J segments at the nucleotide or
### amino acid levels.
### Returns a data.frame or tibble parallel to 'AIRR_df'.
percent_mutation <- function(AIRR_df, for.aa=FALSE, as.matrix=FALSE)
{
    if (!is.data.frame(AIRR_df))
        stop(wmsg("'AIRR_df' must be data.frame or tibble ",
                  "as returned by igblastn()"))
    if (!isTRUEorFALSE(for.aa))
        stop(wmsg("'for.aa' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(as.matrix))
        stop(wmsg("'as.matrix' must be TRUE or FALSE"))
    sequence_id <- AIRR_df[["sequence_id"]]
    if (is.null(sequence_id))
        .stop_on_missing_AIRR_cols("sequence_id")
    v_perc_mut <- .region_percent_mutation(AIRR_df, "V", for.aa=for.aa)
    d_perc_mut <- .region_percent_mutation(AIRR_df, "D", for.aa=for.aa)
    j_perc_mut <- .region_percent_mutation(AIRR_df, "J", for.aa=for.aa)
    if (as.matrix) {
        ans <- cbind(v_perc_mut=v_perc_mut,
                     d_perc_mut=d_perc_mut,
                     j_perc_mut=j_perc_mut)
        rownames(ans) <- sequence_id
    } else {
        ans <- data.frame(v_perc_mut=v_perc_mut,
                          d_perc_mut=d_perc_mut,
                          j_perc_mut=j_perc_mut)
    }
    if (for.aa)
        colnames(ans) <- paste0(colnames(ans), "_aa")
    if (as.matrix)
        return(ans)
    locus <- AIRR_df[["locus"]]
    if (is.null(locus))
        .stop_on_missing_AIRR_cols("locus")
    ans <- cbind(data.frame(sequence_id=sequence_id, locus=locus), ans)
    if (is_tibble(AIRR_df))
        ans <- as_tibble(ans)
    ans
}

