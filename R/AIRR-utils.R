### =========================================================================
### Low-level utilities to retrieve data from the AIRR-community/OGRDB
### download site
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


.OGRDB_URL <- "https://ogrdb.airr-community.org/"

### See short introduction to OGRDB REST API at
###   https://wordpress.vdjbase.org/index.php/ogrdb_news/downloading-germline-sets-from-the-command-line-or-api/
.OGRDB_API_URL <- paste0(.OGRDB_URL, "api")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .fetch_germline_set_from_OGRDB()
###

.encode_OGRDB_API_query <- function(query)
{
    stopifnot(is.character(query))
    if (length(query) != 0L)
        stopifnot(!is.null(names(query)))
    vapply(query, function(x) gsub("/", "%25252f", URLencode(x)), character(1))
}

### Retrieves:
###     /germline/set/
###       {species}/{set_name}/{release_version}/{format}
###     or
###       {species}/{species_subgroup}/{set_name}/{release_version}/{format}
### as documented at https://ogrdb.airr-community.org/api/
###
### Typical usage:
###
###   .fetch_germline_set_from_OGRDB("Human", set_name="IGH_VDJ")
###   .fetch_germline_set_from_OGRDB("Mouse", species_subgroup="C57BL/6",
###                                  set_name="C57BL/6 IGH")
###
### Note that passing "Homo_sapiens" or "Mus musculus" instead of "Human"
### or "Mouse" above also works and produces the same results.
### Returns a character vector containing the nucleotide sequences in FASTA
### format, except when 'format' is set to "airr" or "airr_ex".
###
### About the "airr" and "airr_ex" formats
### --------------------------------------
###
### When 'format' is set to "airr" or "airr_ex", the function returns a
### named list containing a bunch of information about the germline set
### in adition to the nucleotide sequences.
### If 'res' is the results obtained with 'format="airr"' then all the
### sequences in the germline set are described in 'res$allele_descriptions'
### which is a data.frame with 1 row per sequence and dozens of columns.
### However, for some mysterious reason, this data.frame has more rows than
### the number of sequences returned when 'format' is not "airr" or "airr_ex".
### For example, 'res$allele_descriptions' has 245 rows for Human/IGH_VDJ
### but using 'format="ungapped"' returns only 236 sequences!
###
### Assuming 'alleles' is the data.frame obtained by removing the extra
### sequences from 'res$allele_descriptions', then the seq ids, species,
### species subgroups, loci, nucleotide sequences (ungapped and gapped),
### and region types, can be extracted with:
###
###   - seqids <- alleles$label
###
###   - species <- alleles$species$label
###
###   - species_subgroups <- alleles$species_subgroup
###
###   - loci <- alleles$locus
###
###   - region_types <- alleles$sequence_type
###
###   - ungapped_seqs <- alleles$coding_sequence
###     Note that for sequences of type V, this should match:
###       sapply(alleles$v_gene_delineations,
###         function(x) if (is.null(x)) NA_character_ else x$unaligned_sequence)
###
###   - gapped_seqs <- sapply(alleles$v_gene_delineations,
###         function(x) if (is.null(x)) NA_character_ else x$aligned_sequence)
###     Note that only sequences of type V can have gaps.
.fetch_germline_set_from_OGRDB <-
    function(species, species_subgroup=NULL, set_name,
             release_version="published",
             format=c("ungapped", "gapped", "airr",
                      "ungapped_ex", "gapped_ex", "airr_ex"))
{
    stopifnot(isSingleNonWhiteString(species),
              is.null(species_subgroup) ||
                  isSingleNonWhiteString(species_subgroup),
              isSingleNonWhiteString(set_name),
              isSingleNonWhiteString(release_version))
    format <- match.arg(format)

    query <- c(species=species, species_subgroup=species_subgroup,
               set_name=set_name, release_version=release_version,
               format=format)
    query <- .encode_OGRDB_API_query(query)
    url <- paste(c(.OGRDB_API_URL, "germline/set", query), collapse="/")
    content <- getUrlContent(url, type="text", encoding="UTF-8")
    ### Is it possible that 'content' will be NA? See R/REST_API.R in
    ### the UCSC.utils package for more info.
    stopifnot(isSingleString(content))

    if (!(format %in% c("airr", "airr_ex"))) {
        fasta_lines <- strsplit(content, split="\n", fixed=TRUE)[[1L]]
        return(fasta_lines)  # nucleotide sequences in FASTA format
    }

    parsed_json <- fromJSON(content, simplifyDataFrame=FALSE)
    ## Sanity checks.
    stopifnot(is.list(parsed_json),
              length(parsed_json) == 1L,
              identical(names(parsed_json), "GermlineSet"),
              is.list(parsed_json[[1L]]),
              length(parsed_json[[1L]]) == 1L)
    ans <- parsed_json[[1L]][[1L]]
    stopifnot(is.list(ans), !is.null(names(ans)),
              "allele_descriptions" %in% names(ans))
    ans$allele_descriptions <-
        jsonlite:::simplifyDataFrame(ans$allele_descriptions,
                                     flatten=FALSE, simplifyMatrix=TRUE)
    ans  # named list
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_germline_sequences_from_OGRDB()
###

.infer_region_type_from_seqid <- function(seqid, locus)
{
    stopifnot(is.character(seqid), isSingleNonWhiteString(locus),
              all(has_prefix(seqid, locus)))
    pos <- nchar(locus) + 1L
    region_type <- substr(seqid, pos, pos)
    stopifnot(all(region_type %in% VDJ_REGION_TYPES))
    region_type
}

.check_ungapped_seqs <- function(ungapped_seqs,
                                 species, species_subgroup=NULL, set_name,
                                 release_version="published",
                                 extended=FALSE)
{
    format <- if (extended) "airr_ex" else "airr"
    parsed_json <- .fetch_germline_set_from_OGRDB(species,
                        species_subgroup=species_subgroup,
                        set_name=set_name,
                        release_version=release_version,
                        format=format)
    stop("paranoid.mode not ready yet")
}

### Downloads the germline sequences (ungapped) for the specified
### species/species_subgroup/set_name/release_version.
### Produces a subset of the following 7 FASTA files: IGHV.fasta,
### IGHD.fasta, IGHJ.fasta, IGKV.fasta, IGKJ.fasta, IGLV.fasta, and
### IGLJ.fasta. The exact subset produced depends on the species/set_name.
### Returns the number of files produced.
download_germline_sequences_from_OGRDB <-
    function(species, species_subgroup=NULL, set_name,
             locus=c("IGH", "IGK", "IGL"),
             release_version="published",
             extended=FALSE, destdir=".",
             paranoid.mode=FALSE)
{
    stopifnot(isSingleNonWhiteString(release_version),
              isTRUEorFALSE(extended),
              isSingleNonWhiteString(destdir),
              isTRUEorFALSE(paranoid.mode))
    locus <- match.arg(locus)

    format <- if (extended) "ungapped_ex" else "ungapped"
    fasta_lines <- .fetch_germline_set_from_OGRDB(species,
                        species_subgroup=species_subgroup,
                        set_name=set_name,
                        release_version=release_version,
                        format=format)
    filepath <- tempfile(fileext=".fasta")
    writeLines(fasta_lines, filepath)
    ungapped_seqs <- readDNAStringSet(filepath)
    seq_region_types <- .infer_region_type_from_seqid(names(ungapped_seqs),
                                                      locus)

    if (paranoid.mode)
        .check_ungapped_seqs(ungapped_seqs,
                             species, species_subgroup=species_subgroup,
                             set_name=set_name,
                             release_version=release_version,
                             extended=extended)

    from <- paste0("germline set ", set_name, " (", species, ")")
    file_count <- 0L
    for (region_type in VDJ_REGION_TYPES) {
        selected_seqs <- ungapped_seqs[seq_region_types == region_type]
        if (length(selected_seqs) == 0L)
            next
        filename <- paste0(locus, region_type, ".fasta")
        destfile <- file.path(destdir, filename)
        message("Write ", length(selected_seqs), " ", region_type, " regions ",
                "from ", from, " to ", filename, " ... ", appendLF=FALSE)
        if (file.exists(destfile))
            stop(wmsg(filename, " file already exists in ", destdir, "/"))
        writeXStringSet(selected_seqs, destfile)
        message("ok")
        file_count <- file_count + 1L
    }
    file_count
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_human_germline_sequences_from_OGRDB()
###

### Downloads the germline sequences (ungapped) for all the Human germline
### sets. See list of Human germline sets here:
###   https://ogrdb.airr-community.org/germline_sets/Homo%20sapiens
### Returns the number of files produced (should be 7).
download_human_germline_sequences_from_OGRDB <-
    function(release_version="published",
             extended=FALSE, destdir=".",
             paranoid.mode=FALSE)
{
    stopifnot(isSingleNonWhiteString(release_version),
              isTRUEorFALSE(extended),
              isSingleNonWhiteString(destdir),
              isTRUEorFALSE(paranoid.mode))

    ## Produces files: IGHV.fasta, IGHD.fasta, IGHJ.fasta
    file_count <-
        download_germline_sequences_from_OGRDB("Human",
                    set_name="IGH_VDJ", locus="IGH",
                    release_version=release_version, extended=extended,
                    destdir=destdir, paranoid.mode=paranoid.mode)
    ## Produces files: IGKV.fasta, IGKJ.fasta
    file_count <- file_count +
        download_germline_sequences_from_OGRDB("Human",
                    set_name="IGKappa_VJ", locus="IGK",
                    release_version=release_version, extended=extended,
                    destdir=destdir, paranoid.mode=paranoid.mode)
    ## Produces files: IGLV.fasta, IGLJ.fasta
    file_count <- file_count +
        download_germline_sequences_from_OGRDB("Human",
                    set_name="IGLambda_VJ", locus="IGL",
                    release_version=release_version, extended=extended,
                    destdir=destdir, paranoid.mode=paranoid.mode)
    file_count
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_mouse_germline_sequences_from_OGRDB()
###

### 'set_names' is expected to be a vector of Set Names from
### the Germline Sets table displayed at
###   https://ogrdb.airr-community.org/germline_sets/Mus%20musculus
### except special Set Names "IGKJ (all strains)" and "IGLJ (all strains)".
### Assumes that all the Set Names in 'set_names' are of the
### form "<Species subgroup> <group>", where <group> is
### one of IGH, IGKV, an IGLV.
### Returns a 4-column matrix with colum names: set_name, species_subgroup,
### group, and locus.
.set_names_as_matrix <- function(set_names)
{
    if (!is.character(set_names) || anyNA(set_names))
        stop(wmsg("'set_names' must be a character vector with no NAs"))
    set_names <- trimws2(set_names)
    if (anyDuplicated(set_names))
        stop(wmsg("'set_names' cannot contain duplicates"))
    FORBIDDEN_SET_NAMES <- c("IGKJ (all strains)", "IGLJ (all strains)")
    if (any(tolower(set_names) %in% tolower(FORBIDDEN_SET_NAMES))) {
        in1string <- paste0('"', FORBIDDEN_SET_NAMES, '"', collapse=" or ")
        stop(wmsg("'set_names' cannot contain ", in1string))
    }
    split_set_names <- strsplit(set_names, split=" ", fixed=TRUE)
    if (!all(lengths(split_set_names) == 2L))
        stop(wmsg("each Set Name in 'set_names' must contain exactly 1 space"))
    data <- unlist(split_set_names)
    if (is.null(data))
        data <- character(0)
    m <- matrix(data, ncol=2L, byrow=TRUE,
                dimnames=list(NULL, c("species_subgroup", "group")))
    loci <- substr(m[ , "group"], 1L, 3L)
    bad_idx <- which(!(loci %in% c("IGH", "IGK", "IGL")))
    if (length(bad_idx) != 0L) {
        in1string <- paste(set_names[bad_idx], collapse=", ")
        stop(wmsg("Bad Set Name(s): ", in1string, "."),
             "\n  ",
             wmsg("The first 3 letters of the 2nd part of the name ",
                  "is not a valid locus name."))
    }
    set_name <- matrix(set_names, dimnames=list(NULL, "set_name"))
    locus    <- matrix(loci, dimnames=list(NULL, "locus"))
    ans <- cbind(set_name, m, locus)
    species_subgroup <- unique(ans[ , "species_subgroup"])
    if (length(species_subgroup) > 1L)
        warning(wmsg("the supplied Set Names are from ",
                     "more than one Species subgroup"))
    ans
}

### Downloads the germline sequences (ungapped) for the specified Mouse
### germline sets. See list of Mouse germline sets here:
###   https://ogrdb.airr-community.org/germline_sets/Mus%20musculus
### Returns the number of files produced (will be typically between 4 and 7).
###
### To download all the germline sequences for Mouse strain A/J:
###
###   set_names <- c("A/J IGKV", "A/J IGLV")
###   download_mouse_germline_sequences_from_OGRDB(set_names)
###
### --> produces 4 FASTA files (partial germline db).
###
### To download all the germline sequences for Mouse strain C57BL/6:
###
###   download_mouse_germline_sequences_from_OGRDB("C57BL/6 IGH")
###
### --> produces 5 FASTA files (partial germline db).
###
### To download all the germline sequences for Mouse strain C57BL/6J:
###
###   set_names <- c("C57BL/6J IGKV", "C57BL/6J IGLV")
###   download_mouse_germline_sequences_from_OGRDB(set_names)
###
### --> produces 4 FASTA files (partial germline db).
###
### To download all the germline sequences for Mouse strain CAST/EiJ:
###
###   set_names <- c("CAST/EiJ IGH", "CAST/EiJ IGKV", "CAST/EiJ IGLV")
###   download_mouse_germline_sequences_from_OGRDB(set_names)
###
### --> produces 7 FASTA files (full germline db).
download_mouse_germline_sequences_from_OGRDB <-
    function(set_names,
             release_version="published",
             extended=FALSE, destdir=".",
             paranoid.mode=FALSE)
{
    stopifnot(isSingleNonWhiteString(release_version),
              isTRUEorFALSE(extended),
              isSingleNonWhiteString(destdir),
              isTRUEorFALSE(paranoid.mode))

    m <- .set_names_as_matrix(set_names)
    set_names         <- m[ , "set_name"]
    species_subgroups <- m[ , "species_subgroup"]
    loci              <- m[ , "locus"]
    file_count <- 0L
    for (i in seq_len(nrow(m))) {
        file_count <- file_count +
            download_germline_sequences_from_OGRDB("Mouse",
                        species_subgroup=species_subgroups[[i]],
                        set_name=set_names[[i]], locus=loci[[i]],
                        release_version=release_version, extended=extended,
                        destdir=destdir, paranoid.mode=paranoid.mode)
    }
    file_count <- file_count +
        download_germline_sequences_from_OGRDB("Mouse",
                    set_name="IGKJ (all strains)", locus="IGK",
                    release_version=release_version, extended=extended,
                    destdir=destdir, paranoid.mode=paranoid.mode)
    file_count <- file_count +
        download_germline_sequences_from_OGRDB("Mouse",
                    set_name="IGLJ (all strains)", locus="IGL",
                    release_version=release_version, extended=extended,
                    destdir=destdir, paranoid.mode=paranoid.mode)
    file_count
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_V_ndm_data_from_OGRDB()
###

### 'germline_sets' must be a named integer vector where the names are
### the germline sets to download and the values are their versions.
download_V_ndm_data_from_OGRDB <-
    function(organism, germline_sets,
             json_file=c("auto", "airr_ex", "airr"), check.data=FALSE, ...)
{
    stopifnot(isSingleNonWhiteString(organism),
              is.numeric(germline_sets),
              isTRUEorFALSE(check.data))
    sets <- names(germline_sets)
    if (is.null(sets))
        stop(wmsg("'germline_sets' must have names on it"))
    json_file <- match.arg(json_file)
    if (json_file == "auto") {
        is_human <- grepl("Homo.sapiens", organism, ignore.case=TRUE)
        json_file <- ifelse(is_human, "airr_ex", "airr")
    }
    base_url <- paste0(.OGRDB_URL, "download_germline_set/")
    tmp_json_file <- tempfile()
    for (i in seq_along(germline_sets)) {
        set <- sets[[i]]  # name of germline set
        version <- germline_sets[[i]]
        url <- sprintf("%s%s/%s/%s/%s", base_url, organism, set,
                                        version, json_file)
        download.file(URLencode(url), tmp_json_file, ...)
        V_ndm_data <- makeogrannote(tmp_json_file)
        ## We infer the locus from the name of the germline set.
        destfile <- paste0(substr(set, 1L, 3L), "V.ndm.imgt")
        message("Writing ", destfile, " (", nrow(V_ndm_data), " rows) ... ",
                appendLF=FALSE)
        write_V_ndm_data(V_ndm_data, destfile, check.data=check.data)
        message("ok\n")
    }
}

