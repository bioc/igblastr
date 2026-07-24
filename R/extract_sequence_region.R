### =========================================================================
### Extract region of interest from a set of analyzed sequences
### -------------------------------------------------------------------------


stop_on_missing_AIRR_cols <- function(missing_colnames)
{
    in1string <- paste(paste0("\"", missing_colnames, "\""), collapse=", ")
    what <- if (length(missing_colnames) == 1L) " is" else "s are"
    msg <- c("The following column", what, " missing in the supplied ",
             "data.frame or tibble: ", in1string, ". Was it obtained ",
             "with igblastn()?")
    stop(wmsg(msg))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_sequence_region()
### extract_region_length()
###

.check_region_type <- function(region_type)
{
    if (!isSingleNonWhiteString(region_type))
        stop(wmsg("'region_type' must be a single (non-empty) string"))
    valid_region_types <- c(VDJC_REGION_TYPES, FWRCDR_NAMES)
    if (!(tolower(region_type) %in% tolower(valid_region_types))) {
        in1string <- paste0("\"", valid_region_types, "\"", collapse=", ")
        stop(wmsg("'region_type' must be one of ", in1string))
    }
}

### Returns a character vector with one nucleotide sequence per row
### in 'AIRR_df'.
.extract_sequence_region_as_nuc <- function(AIRR_df, region_type)
{
    stopifnot(is.data.frame(AIRR_df))
    .check_region_type(region_type)
    sequences <- AIRR_df[["sequence"]]
    if (is.null(sequences))
        stop_on_missing_AIRR_cols("sequence")
    prefix <- tolower(region_type)
    if (prefix %in% tolower(VDJC_REGION_TYPES))
        prefix <- paste0(prefix, "_sequence")
    COI <- paste0(prefix, c("_start", "_end"))
    missing_idx <- which(!(COI %in% colnames(AIRR_df)))
    if (length(missing_idx) != 0L)
        stop_on_missing_AIRR_cols(COI[missing_idx])
    region_starts <- AIRR_df[[COI[[1L]]]]
    region_ends <- AIRR_df[[COI[[2L]]]]
    substr(sequences, region_starts, region_ends)
}

### Returns an AAStringSet object with one amino acid sequence per row
### in 'AIRR_df'.
.extract_sequence_region_as_aa <- function(AIRR_df, region_type)
{
    stopifnot(is.data.frame(AIRR_df))
    .check_region_type(region_type)
    prefix <- tolower(region_type)
    if (prefix %in% FWRCDR_NAMES) {
        ## 'aa' is a character vector that should never contain NAs.
        aa <- AIRR_df[[paste0(prefix, "_aa")]]
    } else {
        ## 'aa' is a character vector that can contain NAs. Since AASringSet
        ## objects don't support them, we replace them with empty strings.
        aa <- AIRR_df[[paste0(prefix, "_sequence_alignment_aa")]]
        aa[is.na(aa)] <- ""
        aa <- gsub("-", "", aa)  # remove the gaps if any
    }
    AAStringSet(aa)
}

### Returns a character vector (when 'as.aa' is FALSE) or an AAStringSet
### object (when 'as.aa' is TRUE) with one element per row in 'AIRR_df'.
extract_sequence_region <- function(AIRR_df, region_type, as.aa=FALSE)
{
    if (!is.data.frame(AIRR_df))
        stop(wmsg("'AIRR_df' must be data.frame or tibble ",
                  "as returned by igblastn()"))
    if (!isTRUEorFALSE(as.aa))
        stop(wmsg("'as.aa' must be TRUE or FALSE"))
    if (as.aa)
        return(.extract_sequence_region_as_aa(AIRR_df, region_type))
    ans <- .extract_sequence_region_as_nuc(AIRR_df, region_type)
    region_type <- tolower(region_type)
    if (region_type %in% FWRCDR_NAMES) {
        ## The sequences in 'ans' were obtained by extracting the substrings
        ## from 'AIRR_df$sequence' that are defined by the start/end positions
        ## found in the appropriate '*_start' and '*_end' columns.
        ## However, in the case of the FWR/CDR regions, the sequences of the
        ## regions are also provided in their own dedicated 'AIRR_df' columns.
        ## As a sanity check, we verify that these are the same as the
        ## sequences in 'ans'.
        stopifnot(identical(ans, AIRR_df[[region_type]]))
    }
    ans
}

extract_region_length <- function(AIRR_df, region_type)
{
    nchar(extract_sequence_region(AIRR_df, region_type))
}

