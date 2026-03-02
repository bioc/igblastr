### =========================================================================
### Low-level utilities to retrieve data from the AIRR-community/OGRDB
### download site
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


.OGRDB_URL <- "https://ogrdb.airr-community.org/"

### For some reason OGRDB URLs are expected to have the forward slash ("/")
### in their components encoded with "%25252f" instead of the usual "%2f".
### This means that we cannot simply use 'URLencode(x, reserved=TRUE)' to
### encode "/" in the components.
.encode_OGRDB_URL_component <- function(component)
{
    stopifnot(isSingleString(component))
    gsub("/", "%25252f", URLencode(component))
}

.OGRDB_organism_url <- function(organism, for.download=FALSE)
{
    stopifnot(isSingleNonWhiteString(organism), isTRUEorFALSE(for.download))
    top <- if (for.download) "download_germline_set" else "germline_sets"
    paste0(.OGRDB_URL, top, "/", .encode_OGRDB_URL_component(organism))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normalize_OGRDB_organism()
### infer_OGRDB_species_subgroup()
### normalize_OGRDB_format()
### OGRDB_format2fileext()
###

### Organisms available at OGRDB as of February 2026.
.OGRDB_ORGANISMS <- c(
    "Homo sapiens",
    "Macaca mulatta",
    "Mus musculus"
)

normalize_OGRDB_organism <- function(organism)
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    organism <- chartr("_", " ", organism)
    m <- match(tolower(organism), tolower(.OGRDB_ORGANISMS))
    if (is.na(m)) {
        in1string <- paste0("\"", .OGRDB_ORGANISMS, "\"", collapse=", ")
        stop(wmsg("'organism' should be one of ", in1string))
    }
    .OGRDB_ORGANISMS[[m]]
}

### Infer the "Species subgroup" from the organism/set_name combination.
### Note that OGRDB only defines the "Species subgroup" for various mouse
### strains at the moment. In this case the "Species subgroup" is simply
### the name of the strain. See:
###   https://ogrdb.airr-community.org/germline_sets/Mus%20musculus
### Returns the empty string ("") if the organism/set_name combination does
### not correspond to any subgroup.
### Performs various sanity checks that should never fail.
infer_OGRDB_species_subgroup <- function(organism, set_name)
{
    stopifnot(isSingleNonWhiteString(organism))
    if (!isSingleNonWhiteString(set_name))
        stop(wmsg("'set_name' must be a single (non-empty) string"))
    if (organism != "Mus musculus")
        return("")
    stop_if_invalid_set_name <- function() {
        org_url <- .OGRDB_organism_url(organism)
        msg1 <- c("Don't know how to infer \"Species subgroup\" ",
                  "for ", organism, " germline set \"", set_name, "\".")
        msg2 <- c("Please visit ", org_url, " and make sure that the ",
                  "name of the germline set is spelled correctly.")
        stop(wmsg(msg1), "\n  ", wmsg(msg2))
    }
    parts <- strsplit(set_name, " ", fixed=TRUE)[[1L]]
    if (length(parts) == 3L) {
        if (parts[[2L]] != "(all" || parts[[3L]] != "strains)")
            stop_if_invalid_set_name()
        return("")
    }
    if (length(parts) == 2L) {
        if (substr(parts[[2L]], 1L, 2L) != "IG")
            stop_if_invalid_set_name()
        return(parts[[1L]])
    }
    stop_if_invalid_set_name()
}

normalize_OGRDB_format <- function(format, organism, source_set=FALSE)
{
    if (!isSingleNonWhiteString(format))
        stop(wmsg("'format' must be a single (non-empty) string"))
    stopifnot(isSingleNonWhiteString(organism))
    if (!isTRUEorFALSE(source_set))
        stop(wmsg("'source_set' must be TRUE or FALSE"))

    if (source_set) {
        if (organism != "Homo sapiens")
            stop(wmsg("'source_set=TRUE' is only supported for Homo sapiens"))
    } else {
        if (organism == "Homo sapiens")
            format <- paste0(format, "_ex")
    }
    format
}

OGRDB_format2fileext <- function(format)
{
    stopifnot(isSingleNonWhiteString(format))
    switch(format, airr=".json", gapped=, ungapped=".fasta",
           stop(wmsg("unknown OGRDB format: ", format)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_OGRDB_germline_set()
###

.build_OGRDB_germline_set_url <- function(organism, set_name, set_version,
                                          format)
{
    subgroup <- infer_OGRDB_species_subgroup(organism, set_name)
    org_url <- .OGRDB_organism_url(organism, for.download=TRUE)
    if (subgroup != "")
        org_url <- paste0(org_url, "/", .encode_OGRDB_URL_component(subgroup))
    sprintf("%s/%s/%s/%s", org_url, .encode_OGRDB_URL_component(set_name),
                           set_version, format)
}

### Note that the OGRDB website does NOT return an "HTTP 404 Not Found" error
### if the requested germline set does not exist. Instead it redirects to
### the website home page. However, looking at the headers of a HEAD requests
### gives us a clue.
.OGRDB_germline_set_exists <- function(url, ...)
{
    if (!has_internet())
        stop(wmsg("no internet"))
    config <- config(...)
    response <- try(HEAD(url, config, user_agent("igblastr")), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(wmsg(as.character(response)))
    if (response$status_code == 404L)
        stop(wmsg("Not Found (HTTP 404): ", url))
    stop_for_status(response)
    headers <- headers(response)
    content_type <- headers[["content-type"]]
    content_disposition <- headers[["content-disposition"]]
    if (is.null(content_type) || is.null(content_disposition))
        return(FALSE)
    identical(content_type, "application/octet-stream")
}

### Similar to .fetch_germline_set_from_OGRDB() below but doesn't use
### OGRDB API.
### Returns the filename of the downloaded germline set.
download_OGRDB_germline_set <-
    function(organism, set_name, set_version,
             format=c("airr", "gapped", "ungapped"), source_set=FALSE,
             destdir=".", ...)
{
    ## Check arguments.
    organism <- normalize_OGRDB_organism(organism)
    if (!isSingleNonWhiteString(set_name))
        stop(wmsg("'set_name' must be a single (non-empty) string"))
    if (!isSingleNumber(set_version))
        stop(wmsg("'set_version' must be a single number"))
    if (!is.integer(set_version))
        set_version <- as.integer(set_version)
    format <- match.arg(format)
    fileext <- OGRDB_format2fileext(format)
    format <- normalize_OGRDB_format(format, organism, source_set=source_set)
    if (!isSingleNonWhiteString(destdir))
        stop(wmsg("'destdir' must be a single (non-empty) string"))
    if (!dir.exists(destdir)) {
        if (file.exists(destdir))
            stop(wmsg(destdir, ": not a directory"))
        stop(wmsg(destdir, ": no such directory"))
    }

    url <- .build_OGRDB_germline_set_url(organism, set_name, set_version,
                                         format)
    if (!.OGRDB_germline_set_exists(url)) {
        org_url <- .OGRDB_organism_url(organism)
        msg1 <- c(organism, " germline set \"", set_name, "\" ",
                  "(version ", set_version, ") not found.")
        msg2 <- c("Please visit ", org_url, " and make sure that the ",
                  "name of the germline set is spelled correctly and that ",
                  "the version is valid.")
        stop(wmsg(msg1), "\n  ", wmsg(msg2))
    }
    filename <- paste0(format, fileext)
    download_file(url, file.path(destdir, filename), ...)
    filename
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .fetch_germline_set_from_OGRDB()
###

### See short introduction to OGRDB REST API at
###   https://wordpress.vdjbase.org/index.php/ogrdb_news/downloading-germline-sets-from-the-command-line-or-api/
.OGRDB_API_URL <- paste0(.OGRDB_URL, "api")

.encode_OGRDB_API_query <- function(query)
{
    stopifnot(is.character(query))
    if (length(query) != 0L)
        stopifnot(!is.null(names(query)))
    vapply(query, .encode_OGRDB_URL_component, character(1))
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
             format=c("gapped", "ungapped", "airr",
                      "gapped_ex", "ungapped_ex", "airr_ex"))
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
### NOTE: Superseded by download_OGRDB_germline_sequences()!
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
###
### NOTE: No longer needed! Use download_OGRDB_germline_sequences() instead.
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
###
### NOTE: No longer needed! Use download_OGRDB_germline_sequences() instead.
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

