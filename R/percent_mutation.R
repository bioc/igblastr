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

.compute_hamming_dist <- function(sequence_aln, germline_aln, for.aa=FALSE)
{
    stopifnot(is.character(sequence_aln), !anyNA(sequence_aln),
              is.character(germline_aln), !anyNA(germline_aln),
              identical(nchar(sequence_aln), nchar(germline_aln)),
              isTRUEorFALSE(for.aa))
    if (for.aa) {
        sequence_aln <- AAStringSet(sequence_aln)
        germline_aln <- AAStringSet(germline_aln)
    } else {
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
    100 * d / nchar(sequence_aln)  # guaranteed to be = nchar(germline_aln)
}

### Computes percent mutation in V, D, J segments at the nucleotide or
### amino acid levels.
### Returns a data.frame or tibble parallel to 'AIRR_df'.
percent_mutation <- function(AIRR_df, for.aa=FALSE)
{
    if (!is.data.frame(AIRR_df))
        stop(wmsg("'AIRR_df' must be data.frame or tibble ",
                  "as returned by igblastn()"))
    sequence_id <- AIRR_df[["sequence_id"]]
    if (is.null(sequence_id))
        .stop_on_missing_AIRR_cols("sequence_id")
    locus <- AIRR_df[["locus"]]
    if (is.null(locus))
        .stop_on_missing_AIRR_cols("sequence_id")
    v_perc_mut <- .region_percent_mutation(AIRR_df, "V", for.aa=for.aa)
    d_perc_mut <- .region_percent_mutation(AIRR_df, "D", for.aa=for.aa)
    j_perc_mut <- .region_percent_mutation(AIRR_df, "J", for.aa=for.aa)
    ans <- data.frame(sequence_id=sequence_id,
                      locus=locus,
                      v_perc_mut=v_perc_mut,
                      d_perc_mut=d_perc_mut,
                      j_perc_mut=j_perc_mut)
    if (is_tibble(AIRR_df))
        ans <- as_tibble(ans)
    ans
}

