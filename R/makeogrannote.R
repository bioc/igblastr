### =========================================================================
### makeogrannote()
### -------------------------------------------------------------------------


### Not the true colnames used in IgBLAST intdata files: ours are all
### lowercase and we've replaced spaces with underscores.
### Note that many columns are redundant:
### - columns 'cdr1_start', 'fwr2_start', 'cdr2_start', and 'fwr3_start'
###   are redundant with columns 'fwr1_end', 'cdr1_end', 'fwr2_end',
###   and 'cdr2_end', respectively;
### - columns 'fwr1_start' and 'coding_frame_start' are redundant (and
###   column 'fwr1_start' is a dumb column anyways because it should always
###   be set to 1).
.IGBLAST_INTDATA_COL2CLASS <- c(
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

V_GENE_SEGMENTS <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")
V_GENE_DELINEATION_COLNAMES <- paste0(rep(V_GENE_SEGMENTS, each=2L),
                                      c("_start", "_end"))
stopifnot(all(V_GENE_DELINEATION_COLNAMES %in%
              names(.IGBLAST_INTDATA_COL2CLASS)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeogrannote()
###

.extract_v_gene_delineation <- function(allele_description)
{
    stopifnot(is.list(allele_description))
    sequence_type <- allele_description$sequence_type
    if (!identical(sequence_type, "V"))
        return(NULL)
    allele_name <- allele_description$label
    stopifnot(is.character(allele_name))
    locus <- allele_description$locus
    stopifnot(is.character(locus))
    v_gene_delineations <- allele_description$v_gene_delineations
    stopifnot(is.list(v_gene_delineations),
              is.null(names(v_gene_delineations)))
    for (v_gene_delineation in v_gene_delineations) {
        if (!identical(v_gene_delineation$delineation_scheme, "IMGT"))
            next
        stopifnot(all(V_GENE_DELINEATION_COLNAMES %in%
                      names(v_gene_delineation)))
        starts_ends <- v_gene_delineation[V_GENE_DELINEATION_COLNAMES]
        starts_ends <- setNames(as.character(starts_ends), names(starts_ends))
        nc <- nchar(locus)
        chain_type <- paste0(sequence_type, substr(locus, nc, nc))
        return(c(allele_name=allele_name, starts_ends, chain_type=chain_type))
    }
    stop(wmsg("no IMGT v_gene_delineation information found ",
              "for allele ", allele_name))
}

### A reimplementation in R of Python script makeogrannote.py included in
### IgBLAST.
### Returns a data.frame with 1 row per V allele in JSON file 'germline_file'
### and the same columns as the data.frame returned by load_intdata() (see
### R/intdata-utils.R).
makeogrannote <- function(germline_file)
{
    if (!isSingleNonWhiteString(germline_file))
        stop(wmsg("'germline_file' must be a single (non-empty) string"))
    parsed_json <- fromJSON(germline_file, simplifyDataFrame=FALSE)
    stopifnot(is.list(parsed_json), length(parsed_json) == 1L,
              identical(names(parsed_json), "GermlineSet"))
    germline_set <- parsed_json[[1L]]
    stopifnot(is.list(germline_set), length(germline_set) == 1L)
    allele_descriptions <- germline_set[[1L]]$allele_descriptions
    ## 'allele_descriptions' is an unnamed list with 1 list element per allele.
    ## Each element is itself a named list.
    stopifnot(is.list(allele_descriptions),
              is.null(names(allele_descriptions)))
    data <- lapply(allele_descriptions, .extract_v_gene_delineation)
    data <- unlist(data, use.names=FALSE)
    col2class <- head(.IGBLAST_INTDATA_COL2CLASS, n=-1L)
    m <- matrix(data, ncol=length(col2class), byrow=TRUE)
    cbind(matrix2df(m, col2class), coding_frame_start=0L)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### check_V_ndm_data()
###

### Used in tests/testthat/test-auxdata-utils.R!
.rows_with_same_keys_are_identical <- function(df, key)
{
    stopifnot(is.data.frame(df), isSingleNonWhiteString(key))
    keys <- df[ , key]
    m <- match(keys, keys)
    df2 <- df[m, ]
    rownames(df2) <- NULL
    identical(df, df2)
}

.check_region_boundaries <- function(V_ndm_data, region, prev_end)
{
    starts <- V_ndm_data[ , paste0(region, "_start")]
    ends   <- V_ndm_data[ , paste0(region, "_end")]
    (starts == prev_end + 1L) & (ends > starts) & (ends %% 3L == 0L)
}

### Does not check the "chain_type" column at the moment.
check_V_ndm_data <- function(V_ndm_data, allow.dup.entries=FALSE)
{
    if (!is.data.frame(V_ndm_data))
        stop(wmsg("'V_ndm_data' must be a data.frame"))
    if (!isTRUEorFALSE(allow.dup.entries))
        stop(wmsg("'allow.dup.entries' must be TRUE or FALSE"))
    expected_colnames <- names(.IGBLAST_INTDATA_COL2CLASS)
    if (!identical(colnames(V_ndm_data), expected_colnames)) {
        in1string <- paste(expected_colnames, collapse=", ")
        stop(wmsg("'V_ndm_data' must have the following columns ",
                  "(in this order): ", in1string))
    }
    col2classes <- vapply(V_ndm_data, function(x) class(x)[[1L]], character(1))
    if (!identical(col2classes, .IGBLAST_INTDATA_COL2CLASS)) {
        in1string <- paste0("  ", names(.IGBLAST_INTDATA_COL2CLASS), ": ",
                            .IGBLAST_INTDATA_COL2CLASS, collapse="\n")
        stop("'V_ndm_data' must have the following column types:\n", in1string)
    }
    if (allow.dup.entries) {
	## We allow duplicated entries in 'V_ndm_data' as long as they
        ## tell the same story.
        if (!.rows_with_same_keys_are_identical(V_ndm_data, "allele_name"))
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
    read_broken_table(filepath, .IGBLAST_INTDATA_COL2CLASS)
}

write_V_ndm_data <- function(V_ndm_data, file="")
{
    colnames <- names(.IGBLAST_INTDATA_COL2CLASS)
    cat("#", paste(colnames, collapse=", "), "\n", sep="", file=file)
    write.table(V_ndm_data, file, append=TRUE, quote=FALSE,
                sep="\t", row.names=FALSE, col.names=FALSE)
}

