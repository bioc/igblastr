### =========================================================================
### download_OGRDB_germline_json() and related
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_OGRDB_germline_json()
###

.stop_if_duplicated_IGcomps <- function(IGcomps)
{
    m <- match(IGcomps, IGcomps)
    bad_idx <-which(m != seq_along(m))
    if (length(bad_idx) != 0L) {
        idx1 <- bad_idx[[1L]]
        set_name1 <- names(IGcomps)[[m[[idx1]]]]
        set_name2 <- names(IGcomps)[[idx1]]
        msg1 <- c("More than one germline set with the same \"IG* ",
                  "components\": \"", set_name1, "\" and \"", set_name2, "\" ",
                  "(both have the \"", IGcomps[[idx1]], "\" component ",
                  "in their name).")
        msg2 <- c("The specified germline sets must have unique ",
                  "\"IG* components\".")
        stop(wmsg(msg1), "\n  ", wmsg(msg2))
    }
}

download_OGRDB_germline_json <- function(organism, germline_sets,
                                         source_set=FALSE,
                                         destdir=".", overwrite=FALSE,
                                         recache=FALSE, ...)
{
    organism <- normalize_OGRDB_organism(organism)
    germline_sets <- normalize_OGRDB_germline_sets(germline_sets)
    set_names <- names(germline_sets)
    if (!isTRUEorFALSE(source_set))
        stop(wmsg("'source_set' must be TRUE or FALSE"))
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir)) {
        if (file.exists(destdir))
            stop(wmsg(destdir, ": not a directory"))
        stop(wmsg(destdir, ": no such directory"))
    }
    if (!isTRUEorFALSE(overwrite))
        stop(wmsg("'overwrite' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))

    ## We will download one JSON file per supplied germline set and will
    ## use the "IG* component" of each germline set name in 'set_names' to
    ## name the corresponding JSON file. Therefore, the "IG* components"
    ## of the supplied germline sets must be unique.
    IGcomps <- extract_IGcomps_from_OGRDB_set_names(organism, set_names)
    .stop_if_duplicated_IGcomps(IGcomps)

    local_files <- vapply(seq_along(germline_sets),
        function(i) {
            set_name <- set_names[[i]]
            set_version <- germline_sets[[i]]
            ## Download OGRDB germline set to local store if it's not
            ## already there.
            download_OGRDB_germline_set_to_OGRDB_store(organism,
                                    set_name, set_version,
                                    format="airr", source_set=source_set,
                                    recache=recache, ...)
        },
        character(1)
    )

    filenames <- paste0(IGcomps, ".json")  # guaranteed to be unique!
    names(local_files) <- filenames
    copy_files_to_dir(local_files, destdir, overwrite=overwrite)

    invisible(setNames(filenames, set_names))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeogrannote()
###

### Returns an unnamed list with 1 list element per allele.
### Each list element is itself a named list (nested lists).
.extract_allele_descriptions <- function(ogrdb_json_file)
{
    if (!isSingleNonWhiteString(ogrdb_json_file))
        stop(wmsg("'ogrdb_json_file' must be a single (non-empty) string"))
    parsed_json <- fromJSON(ogrdb_json_file, simplifyDataFrame=FALSE)
    stopifnot(is.list(parsed_json), length(parsed_json) == 1L,
              identical(names(parsed_json), "GermlineSet"))
    germline_set <- parsed_json[[1L]]
    stopifnot(is.list(germline_set), length(germline_set) == 1L)
    allele_descriptions <- germline_set[[1L]]$allele_descriptions
    stopifnot(is.list(allele_descriptions),
              is.null(names(allele_descriptions)))
    allele_descriptions
}

.extract_IMGT_v_gene_delineation <- function(allele_description)
{
    stopifnot(is.list(allele_description), !is.null(names(allele_description)))
    sequence_type <- allele_description$sequence_type
    stopifnot(isSingleNonWhiteString(sequence_type),
              sequence_type %in% VDJC_REGION_TYPES)
    if (sequence_type != "V")
        return(NULL)
    allele_name <- allele_description$label
    stopifnot(isSingleNonWhiteString(allele_name))
    locus <- allele_description$locus
    stopifnot(isSingleNonWhiteString(locus),
              locus %in% IG_LOCI)
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
### boundaries found in the JSON file are already relative to the **ungapped**
### V allele sequences. At least that's how they seem to be in 2026.
### Note that the makeogrannote.py script is from 2022 when OGRDB was still
### in its infancy, so maybe the early JSON files that they generated at the
### time had the FR1/CDR1/FR2/CDR2/FR3 boundaries relative to the **gapped**
### V allele sequences? Is this the reason why makeogrannote.py adjusts them?
### Returns a data.frame with 1 row per V allele in JSON file 'ogrdb_json_file'
### and the same columns as the data.frame returned by load_intdata() (see
### R/intdata-utils.R).
makeogrannote <- function(ogrdb_json_file)
{
    allele_descriptions <- .extract_allele_descriptions(ogrdb_json_file)
    data <- lapply(allele_descriptions, .extract_IMGT_v_gene_delineation)
    data <- unlist(data, use.names=FALSE)
    if (is.null(data)) {
        warning(wmsg("no V allele descriptions found ",
                     "in JSON file: ", ogrdb_json_file))
        data <- character(0)
    }
    col2class <- head(NDM_DATA_COL2CLASS, n=-1L)
    m <- matrix(data, ncol=length(col2class), byrow=TRUE)
    df <- matrix2df(m, col2class)
    cbind(df, coding_frame_start=integer(nrow(df)))
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

    intdata <- compute_V_gene_delineations(V_fasta_file)
    locus <- extract_loci_from_OGRDB_set_names(organism, names(germline_set))
    intdata$chain_type <- paste0("V", substr(locus, 3L, 3L))
    intdata[ , names(NDM_DATA_COL2CLASS)]
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
    same_ndm_data(intdata1, intdata2)
}

