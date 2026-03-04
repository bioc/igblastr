### =========================================================================
### makeogrannote() and related
### -------------------------------------------------------------------------
###
### Only makeogrannote() in this file is exported.
###


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
### validate_OGRDB_intdata()
###

### Returns the intdata in a data.frame.
.infer_intdata_from_OGRDB_gapped_V_sequences <-
    function(organism, germline_set, source_set=FALSE, recache=FALSE, ...)
{
    tmp_dir <- tempfile()
    dir.create(tmp_dir)
    on.exit(nuke_file(tmp_dir))
    download_OGRDB_germline_sequences(organism, germline_set,
                                      source_set=source_set,
                                      destdir=tmp_dir,
                                      recache=recache, ...)
    V_fasta_file <- list_fasta_files(tmp_dir, pattern="V\\.fasta$")
    if (length(V_fasta_file) == 0L) {
        msg <- c(organism, " germline set \"", names(germline_set), "\" ",
                 "(version ", germline_set, ") does not seem to contain ",
                 "any V allele")
        stop(wmsg(msg))
    }
    stopifnot(length(V_fasta_file) == 1L)

    intdata <- compute_imgt_intdata(V_fasta_file)
    locus <- infer_loci_from_OGRDB_set_names(names(germline_set))
    intdata$chain_type <- paste0("V", substr(locus, 3L, 3L))
    intdata[ , names(IGBLAST_INTDATA_COL2CLASS)]
}

### Returns the intdata in a data.frame.
.extract_intdata_from_OGRDB_json_file <-
    function(organism, germline_set, source_set=FALSE, recache=FALSE, ...)
{
    ## Download OGRDB germline set to local store if it's not
    ## already there.
    local_file <- download_OGRDB_germline_set_to_OGRDB_store(organism,
                                 names(germline_set), germline_set,
                                 format="airr", source_set=source_set,
                                 recache=recache, ...)
    makeogrannote(local_file)
}

### Returns TRUE if the 2 data.frames contain the same intdata (possibly
### with their rows in different order).
.compare_intdata <- function(intdata1, intdata2)
{
    if (!identical(dim(intdata1), dim(intdata2)))
        return(FALSE)
    allele_names1 <- intdata1[ , "allele_name"]
    allele_names2 <- intdata2[ , "allele_name"]
    if (!setequal(allele_names1, allele_names2))
        return(FALSE)
    m <- match(allele_names1, allele_names2)
    intdata2 <- S4Vectors:::extract_data_frame_rows(intdata2, m)
    identical(intdata1, intdata2)
}

### Not exported!
### Validates the intdata associated with a given OGRDB germline set by
### comparing the two methods of acquisition:
###   1. Infer intdata from the gaps in the V allele sequences.
###   2. Extract intdata from OGRDB json file.
### Returns TRUE if the two methods produce exactly the same intdata, or
### FALSE if they don't.
validate_OGRDB_intdata <- function(organism, germline_set, source_set=FALSE,
                                   recache=FALSE, ...)
{
    organism <- normalize_OGRDB_organism(organism)
    germline_set <- normalize_OGRDB_germline_sets(germline_set)
    stopifnot(length(germline_set) == 1L)
    if (!isTRUEorFALSE(source_set))
        stop(wmsg("'source_set' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))

    ## Method 1: Infer intdata from the gaps in the V allele sequences.
    intdata1 <- .infer_intdata_from_OGRDB_gapped_V_sequences(
                                      organism, germline_set,
                                      source_set=source_set,
                                      recache=recache, ...)


    ## Method 2: Extract intdata from OGRDB json file.
    intdata2 <- .extract_intdata_from_OGRDB_json_file(
                                      organism, germline_set,
                                      source_set=source_set,
                                      recache=recache, ...)

    ## Compare.
    .compare_intdata(intdata1, intdata2)
}

