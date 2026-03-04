### =========================================================================
### Low-level functions to read/write data in "IgBLAST internal data" format
### -------------------------------------------------------------------------
###


### Not exported!
### Not the true colnames used in IgBLAST intdata files: ours are all
### lowercase and we've replaced spaces with underscores.
### Note that many columns are redundant:
### - columns 'cdr1_start', 'fwr2_start', 'cdr2_start', and 'fwr3_start'
###   are redundant with columns 'fwr1_end', 'cdr1_end', 'fwr2_end',
###   and 'cdr2_end', respectively;
### - columns 'fwr1_start' and 'coding_frame_start' are redundant (and
###   column 'fwr1_start' is a dumb column anyways because it should always
###   be set to 1).
IGBLAST_INTDATA_COL2CLASS <- c(
    allele_name="character",
    fwr1_start="integer",
    fwr1_end="integer",
    cdr1_start="integer",
    cdr1_end="integer",
    fwr2_start="integer",
    fwr2_end="integer",
    cdr2_start="integer",
    cdr2_end="integer",
    fwr3_start="integer",
    fwr3_end="integer",
    chain_type="character",
    coding_frame_start="integer"
)

### Not exported!
V_GENE_SEGMENTS <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")
V_GENE_DELINEATION_COLNAMES <- paste0(rep(V_GENE_SEGMENTS, each=2L),
                                      c("_start", "_end"))
stopifnot(all(V_GENE_DELINEATION_COLNAMES %in%
              names(IGBLAST_INTDATA_COL2CLASS)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### check_V_ndm_data_col2class()
###

### Not exported!
check_V_ndm_data_col2class <- function(V_ndm_data, what="'V_ndm_data'")
{
    if (!is.data.frame(V_ndm_data))
        stop(wmsg(what, " must be a data.frame"))
    expected_colnames <- names(IGBLAST_INTDATA_COL2CLASS)
    if (!identical(colnames(V_ndm_data), expected_colnames)) {
        in1string <- paste(expected_colnames, collapse=", ")
        stop(wmsg(what, " must have the following columns ",
                  "(in this order): ", in1string))
    }
    col2class <- vapply(V_ndm_data, function(x) class(x)[[1L]], character(1))
    if (!identical(col2class, IGBLAST_INTDATA_COL2CLASS)) {
        in1string <- paste0("    ", names(IGBLAST_INTDATA_COL2CLASS), " -> ",
                            IGBLAST_INTDATA_COL2CLASS, collapse="\n")
        stop(wmsg(what, " must have the following column types:"), "\n",
             in1string)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### check_V_ndm_data()
###

.check_region_boundaries <- function(V_ndm_data, region, prev_end)
{
    starts <- V_ndm_data[ , paste0(region, "_start")]
    ends   <- V_ndm_data[ , paste0(region, "_end")]
    (starts == prev_end + 1L) & (ends > starts) & (ends %% 3L == 0L)
}

### Does not check the "chain_type" column at the moment.
check_V_ndm_data <- function(V_ndm_data, allow.dup.entries=FALSE)
{
    check_V_ndm_data_col2class(V_ndm_data)
    if (!isTRUEorFALSE(allow.dup.entries))
        stop(wmsg("'allow.dup.entries' must be TRUE or FALSE"))
    if (allow.dup.entries) {
        ## We allow duplicated entries in 'V_ndm_data' as long as they
        ## tell the same story.
        if (!rows_with_same_key_are_identical(V_ndm_data, "allele_name"))
            stop(wmsg("rows in 'V_ndm_data' with same \"allele_name\" ",
                      "must be identical"))
    } else {
        if (anyDuplicated(V_ndm_data[ , "allele_name"]))
            stop(wmsg("'V_ndm_data$allele_name' cannot contain duplicates"))
    }

    fwr1_ok <- .check_region_boundaries(V_ndm_data, "fwr1", 0L)
    cdr1_ok <- .check_region_boundaries(V_ndm_data, "cdr1", V_ndm_data$fwr1_end)
    fwr2_ok <- .check_region_boundaries(V_ndm_data, "fwr2", V_ndm_data$cdr1_end)
    cdr2_ok <- .check_region_boundaries(V_ndm_data, "cdr2", V_ndm_data$fwr2_end)
    fwr3_ok <- .check_region_boundaries(V_ndm_data, "fwr3", V_ndm_data$cdr2_end)
    coding_frame_start_ok <- V_ndm_data[ , "coding_frame_start"] == 0L
    fwr1_ok & cdr1_ok & fwr2_ok & cdr2_ok & fwr3_ok & coding_frame_start_ok
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_V_ndm_data()
### write_V_ndm_data()
###

read_V_ndm_data <- function(filepath)
{
    read_broken_table(filepath, IGBLAST_INTDATA_COL2CLASS)
}

write_V_ndm_data <- function(V_ndm_data, file="", check.data=FALSE)
{
    if (!isTRUEorFALSE(check.data))
        stop(wmsg("'check.data' must be TRUE or FALSE"))
    check_V_ndm_data_col2class(V_ndm_data)
    if (check.data) {
        ok <- check_V_ndm_data(V_ndm_data)
        if (!all(ok))
            stop(wmsg("'V_ndm_data' contains invalid rows. ",
                      "Use 'check_V_ndm_data()' to identify them."))
    }
    header <- paste0("#", paste(colnames(V_ndm_data), collapse=", "))
    cat(header, "\n", sep="", file=file)
    write.table(V_ndm_data, file, append=TRUE, quote=FALSE,
                sep="\t", row.names=FALSE, col.names=FALSE)
}

