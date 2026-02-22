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

V_GENE_SEGMENTS <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")
V_GENE_DELINEATION_COLNAMES <- paste0(rep(V_GENE_SEGMENTS, each=2L),
                                      c("_start", "_end"))
stopifnot(all(V_GENE_DELINEATION_COLNAMES %in%
              names(IGBLAST_INTDATA_COL2CLASS)))


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

### An R reimplementation of Python script makeogrannote.py included in
### IgBLAST. But with the important difference that we don't adjust the
### FR1/CDR1/FR2/CDR2/FR3 boundaries found in the JSON file by subtracting
### the number of gaps that precede them in the corresponding gapped V allele
### sequence, like they do in makeogrannote.py. The reason we don't do this
### is because this adjustment is meant to transform boundaries that are
### relative to the **gapped** V allele sequences into boundaries relative
### to the **ungapped** sequences, which is not necessary because the
### boundaries found in the JSON file are relative to the **ungapped** V
### allele sequences. At least that's how they seem to be in 2026.
### Note that the makeogrannote.py script is from 2022 when OGRDB was still
### in its infancy, so maybe the early JSON files that they generated at the
### time had the FR1/CDR1/FR2/CDR2/FR3 boundaries relative to the **gapped**
### V allele sequences? Is this the reason why makeogrannote.py adjusts them?
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
    col2class <- head(IGBLAST_INTDATA_COL2CLASS, n=-1L)
    m <- matrix(data, ncol=length(col2class), byrow=TRUE)
    cbind(matrix2df(m, col2class), coding_frame_start=0L)
}


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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### add_V_ndm_data_to_germline_db()
###

.write_V_ndm_data_to_germline_db <- function(V_ndm_data, destfile)
{
    check_V_ndm_data_col2class(V_ndm_data)
    stopifnot(isSingleNonWhiteString(destfile))
    internal_data_path <- dirname(destfile)

    ## Check that the set of V alleles annotated in 'V_ndm_data' is the
    ## same as the set of V alleles in germline db file V.fasta.
    db_path <- dirname(internal_data_path)
    db_V_fasta_file <- get_db_fasta_file(db_path, "V")
    allele_names1 <- names(fasta.seqlengths(db_V_fasta_file))
    allele_names2 <- V_ndm_data[ , "allele_name"]
    stopifnot(length(allele_names1) == length(allele_names2),
              setequal(allele_names1, allele_names2))

    ## Reorder the rows in the 'V_ndm_data' data.frame to make it
    ## parallel to 'db_V_fasta_file'.
    m <- match(allele_names1, allele_names2)
    V_ndm_data <- S4Vectors:::extract_data_frame_rows(V_ndm_data, m)

    if (!dir.exists(internal_data_path))
        dir.create(internal_data_path)
    write_V_ndm_data(V_ndm_data, destfile)
}

### Not exported!
add_V_ndm_data_to_germline_db <- function(db_path, V_ndm_store,
                                          domain_system=c("imgt", "kabat"))
{
    stopifnot(isSingleNonWhiteString(db_path), dir.exists(db_path),
              isSingleNonWhiteString(V_ndm_store), dir.exists(V_ndm_store))
    domain_system <- match.arg(domain_system)

    internal_data_path <- file.path(db_path, "internal_data")
    destfile <- file.path(internal_data_path, paste0("V.ndm.", domain_system))
    if (file.exists(destfile))
        return()

    ## Prepare 'V_ndm_data'.
    V_ndm_files <- file.path(V_ndm_store, paste0(IG_LOCI, "V.ndm.imgt"))
    V_ndm_data <- do.call(rbind, lapply(V_ndm_files, read_V_ndm_data))

    .write_V_ndm_data_to_germline_db(V_ndm_data, destfile)
}

