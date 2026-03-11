### =========================================================================
### Download data from OGRDB thru the OGRDB API
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### See short introduction to OGRDB REST API at
###   https://wordpress.vdjbase.org/index.php/ogrdb_news/downloading-germline-sets-from-the-command-line-or-api/
.OGRDB_API_URL <- paste0(OGRDB_URL, "api")

.encode_OGRDB_API_query <- function(query)
{
    stopifnot(is.character(query))
    if (length(query) != 0L)
        stopifnot(!is.null(names(query)))
    vapply(query, encode_OGRDB_URL_component, character(1))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetch_germline_set_from_OGRDB()
###
### No longer used. Superseded by download_OGRDB_germline_set() defined
### in R/OGRDB-utils.R.
###
### Retrieves:
###     /germline/set/
###       {species}/{set_name}/{release_version}/{format}
###     or
###       {species}/{species_subgroup}/{set_name}/{release_version}/{format}
### as documented at https://ogrdb.airr-community.org/api/
###
### Typical usage:
###
###   fetch_germline_set_from_OGRDB("Human", set_name="IGH_VDJ")
###   fetch_germline_set_from_OGRDB("Mouse", species_subgroup="C57BL/6",
###                                 set_name="C57BL/6 IGH")
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

fetch_germline_set_from_OGRDB <-
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

