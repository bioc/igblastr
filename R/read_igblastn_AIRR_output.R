### =========================================================================
### read_igblastn_AIRR_output()
### -------------------------------------------------------------------------


.sanitize_AIRR_df <- function(AIRR_df)
{
    stopifnot(is.data.frame(AIRR_df),
              is.character(AIRR_df[ , "sequence_id"]),
              is.character(AIRR_df[ , "sequence"]),
              is.character(AIRR_df[ , "locus"]))

    ## AIRR field "sequence_aa" was added in igblast 1.21.0.
    sequence_aa <- AIRR_df$sequence_aa
    if (!is.null(sequence_aa))
        stopifnot(is.character(sequence_aa))

    CHARACTER_COLS <- paste0(c("v", "d", "j"), "_call")
    for (colname in CHARACTER_COLS)
        AIRR_df[ , colname] <- as.character(AIRR_df[ , colname])
    c_call <- AIRR_df$c_call
    if (!is.null(c_call))
        AIRR_df$c_call <- as.character(c_call)

    LOGICAL_COLS <- c("stop_codon", "vj_in_frame",
                      "productive", "rev_comp", "complete_vdj")
    for (colname in LOGICAL_COLS)
        AIRR_df[ , colname] <- as.logical(AIRR_df[ , colname])

    ## AIRR field "v_frameshift" was added in igblast 1.17.0.
    v_frameshift <- AIRR_df$v_frameshift
    if (!is.null(v_frameshift))
        AIRR_df$v_frameshift <- as.logical(v_frameshift)

    ## AIRR field "d_frame" was added in igblast 1.21.0.
    d_frame <- AIRR_df$d_frame
    if (!is.null(d_frame))
        AIRR_df$d_frame <- as.logical(d_frame)

    AIRR_df
}

### Read igblastn output format 19 (AIRR format). This is a .tsv file.
read_igblastn_AIRR_output <- function(out)
{
    ## Make sure to use 'tryLogical=FALSE'. By default, i.e.
    ## when 'tryLogical=TRUE', there's the risk that read.table()
    ## will erroneously decide to interpret some columns as logical!
    ## For example this can happen to column "d_sequence_alignment_aa"
    ## if it contains only single-letter amino acid sequences T (Thr)
    ## or F (Phe).
    AIRR_df <- read.table(out, header=TRUE, sep="\t",
                          tryLogical=FALSE, na.strings="")
    .sanitize_AIRR_df(AIRR_df)
}

