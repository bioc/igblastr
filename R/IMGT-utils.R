### =========================================================================
### Low-level utilities to retrieve data from the IMGT/V-QUEST download site
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


IMGT_URL <- "https://www.imgt.org"

### Do not remove the trailing slash.
.VQUEST_DOWNLOAD_ROOT_URL <- paste0(IMGT_URL, "/download/V-QUEST/")

### .VQUEST_REFERENCE_DIRECTORY
.VQUEST_REFERENCE_DIRECTORY <- "IMGT_V-QUEST_reference_directory"

.VQUEST_RELEASE_FILE_URL <-
    paste0(.VQUEST_DOWNLOAD_ROOT_URL, "IMGT_vquest_release.txt")

### Do not remove the trailing slash.
.VQUEST_ARCHIVES_URL <- paste0(.VQUEST_DOWNLOAD_ROOT_URL, "archives/")

get_IMGT_connecttimeout <- function() getOption("IMGT_connecttimeout")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_latest_IMGT_release()
### list_archived_IMGT_zips()
###

.IMGT_cache <- new.env(parent=emptyenv())

.fetch_latest_IMGT_release <- function()
{
    content <- getUrlContent(.VQUEST_RELEASE_FILE_URL,
                             connecttimeout=get_IMGT_connecttimeout())
    sub("^([^ ]*)(.*)$", "\\1", content)
}

get_latest_IMGT_release <- function(recache=FALSE)
{
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    release <- .IMGT_cache[["LATEST_RELEASE"]]
    if (is.null(release) || recache) {
        release <- .fetch_latest_IMGT_release()
        .IMGT_cache[["LATEST_RELEASE"]] <- release
    }
    release
}

### Returns a data.frame with 3 columns (Name, Last modified, Size)
### and 1 row per .zip file.
.fetch_list_of_archived_IMGT_zips <- function()
{
    scrape_html_dir_index(.VQUEST_ARCHIVES_URL,
                          css="body section", suffix=".zip",
                          connecttimeout=get_IMGT_connecttimeout())
}

### If 'as.df' is TRUE then the listing is returned as a data.frame
### with 3 columns (Name, Last modified, Size) and 1 row per .zip file.
list_archived_IMGT_zips <- function(as.df=FALSE, recache=FALSE)
{
    if (!isTRUEorFALSE(as.df))
        stop(wmsg("'as.df' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    listing <- .IMGT_cache[["ARCHIVES_TABLE"]]
    if (is.null(listing) || recache) {
        listing <- .fetch_list_of_archived_IMGT_zips()
        .IMGT_cache[["ARCHIVES_TABLE"]] <- listing
    }
    if (!as.df)
        listing <- listing[ , "Name"]
    listing
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_and_unzip_IMGT_release()
###

.download_and_unzip_latest_IMGT_zip <- function(exdir, ...)
{
    zip_filename <- paste0(.VQUEST_REFERENCE_DIRECTORY, ".zip")

    ## Sometimes, after a new release, the IMGT people forget to make
    ## the zip file of the new release available. We're trying to detect
    ## this and fail graciously when it's the case.
    zip_url <- paste0(.VQUEST_DOWNLOAD_ROOT_URL, zip_filename)
    zip_exists <- urlExists(zip_url, connecttimeout=get_IMGT_connecttimeout())
    if (!zip_exists) {
        release <- get_latest_IMGT_release()
        stop(wmsg("It looks like the zip of the latest IMGT/V-QUEST ",
                  "release (", release, ") is not available (yet?) at ",
                  .VQUEST_DOWNLOAD_ROOT_URL),
             "\n  ",
             wmsg("Please install an older release in the meantime."))
    }

    local_zip <- download_as_tempfile(.VQUEST_DOWNLOAD_ROOT_URL, zip_filename,
                                      ...)
    nuke_file(exdir)
    unzip(local_zip, exdir=exdir)
}

.get_archived_IMGT_zip <- function(release)
{
    stopifnot(isSingleNonWhiteString(release))
    all_zips <- list_archived_IMGT_zips()
    idx <- grep(release, all_zips, fixed=TRUE)
    if (length(idx) == 0L)
        stop(wmsg("Anomaly: no .zip file found at ",
                  .VQUEST_ARCHIVES_URL, " for release ", release))
    if (length(idx) > 1L)
        stop(wmsg("Anomaly: more that one .zip file found at ",
                  .VQUEST_ARCHIVES_URL, " for release ", release))
    all_zips[[idx]]
}

.unzip_archived_IMGT_zip <- function(zipfile, release, exdir)
{
    nuke_file(exdir)
    unzip(zipfile, exdir=exdir, junkpaths=TRUE)
    zip_filename <- paste0(.VQUEST_REFERENCE_DIRECTORY, ".zip")
    local_zip <- file.path(exdir, zip_filename)
    unzip(local_zip, exdir=exdir)
    unlink(local_zip)
}

.download_and_unzip_archived_IMGT_zip <- function(release, exdir, ...)
{
    archived_zip_filename <- .get_archived_IMGT_zip(release)
    archived_zipfile <- download_as_tempfile(.VQUEST_ARCHIVES_URL,
                                             archived_zip_filename, ...)
    .unzip_archived_IMGT_zip(archived_zipfile, release, exdir)
}

### Download and unzip in 'exdir'.
download_and_unzip_IMGT_release <- function(release, exdir, ...)
{
    if (dir.exists(exdir))
        nuke_file(exdir)
    if (release == get_latest_IMGT_release()) {
        .download_and_unzip_latest_IMGT_zip(exdir, ...)
    } else {
        .download_and_unzip_archived_IMGT_zip(release, exdir, ...)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### find_organism_in_IMGT_local_store()
###

list_organisms_in_IMGT_local_store <- function(local_store)
{
    refdir <- file.path(local_store, .VQUEST_REFERENCE_DIRECTORY)
    if (!dir.exists(refdir))
        stop(wmsg("Anomaly: directory ", refdir, " not found"))
    sort(list.files(refdir))
}

### 'local_store' must be the path to the local store of a given IMGT release.
### Returns the path to the subdir of 'local_store' that corresponds to the
### specified organism. For example, for IMGT release 202449-1 and Homo
### sapiens, this path is:
###     <igblastr-cache>
###     └── store
###         └── IMGT-releases
###             └── 202449-1
###                 └── IMGT_V-QUEST_reference_directory
###                     └──  Homo_sapiens
find_organism_in_IMGT_local_store <- function(organism, local_store)
{
    all_organisms <- list_organisms_in_IMGT_local_store(local_store)
    idx <- match(tolower(organism), tolower(all_organisms))
    if (!is.na(idx)) {
        refdir <- file.path(local_store, .VQUEST_REFERENCE_DIRECTORY)
        return(file.path(refdir, all_organisms[[idx]]))
    }
    all_in_1string <- paste0("\"", all_organisms, "\"", collapse=", ")
    stop(wmsg(organism, ": organism not found in ",
              "IMGT/V-QUEST release ", basename(local_store), "."),
         "\n  ",
         wmsg("Available organisms: ", all_in_1string, "."))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level utilities to query IMGT/GENE-DB
###

### IMGT/GENE-DB Query page.
.IMGT_GENE_DB_URL <- "https://www.imgt.org/genedb/"

### IMPORTANT NOTE: Used to map 'seqset_nb' to 'seqset_internal_nb' so
### order is important! See download_C_sequence_set_from_IMGT() below
### for more information about these "sequence set internal numbers".
.SEQSET_SET_INTERNAL_NUMBERS <- c("7.2", "7.5", "7.1", "14.1")

### According to our findings, IMGT/GENE-DB can be queried using an URL
### of the form:
###
###   https://www.imgt.org/genedb/GENElect?query=<query>
###
### where <query> is something like:
###
###   7.2+IGHC&species=Homo+sapiens
###
### The number at the beginning of the query (e.g. 7.2) is an internal
### number used by IMGT to refer to a particular set of sequences.
### See download_C_sequence_set_from_IMGT() below for more information
### about this.
### Returns an ugly HTML page in a character vector and a nucleotide sequence
### embedded in it. Use .scrape_IMGT_GENE_DB_result() below to extract that
### sequence.
.query_IMGT_GENE_DB <- function(species, seqset_internal_nb, group)
{
    stopifnot(isSingleNonWhiteString(species),
              isSingleNonWhiteString(seqset_internal_nb),
              seqset_internal_nb %in% .SEQSET_SET_INTERNAL_NUMBERS,
              isSingleNonWhiteString(group))
    query <- list(query=paste(seqset_internal_nb, group), species=species)
    ## Querying IMGT/GENE-DB can be very slow so we increase the allowed
    ## time by 50%.
    getUrlContent(paste0(.IMGT_GENE_DB_URL, "GENElect"), query=query,
                  type="text", encoding="UTF-8",
                  connecttimeout=get_IMGT_connecttimeout() * 1.5)
}

### 'fasta_lines' must be a character vector.
### Returns FALSE if no FASTA records or if all records are empty.
.is_dna_fasta <- function(fasta_lines)
{
    stopifnot(is.character(fasta_lines))
    header_idx <- grep("^>", fasta_lines)
    if (length(header_idx) == 0L)
        return(FALSE)  # no FASTA records
    dna_lines <- fasta_lines[-header_idx]
    dna <- paste(toupper(dna_lines), collapse="")
    if (!nzchar(dna))
        return(FALSE)  # all records are empty
    all(safeExplode(dna) %in% DNA_ALPHABET)
}

### 'html' is expected to be a character vector containing the HTML document
### returned by .query_IMGT_GENE_DB(). It's expected to contain 2 <pre></pre>
### sections:
###   1. The first one is a section that describes the 15 fields of the
###      FASTA headers.
###   2. The second one contains our nucleotide sequences in FASTA format.
### Instead of assuming that our nucleotide sequences are in the 2nd
### <pre></pre> section, the .scrape_IMGT_GENE_DB_result() function returns
### the content of the first <pre></pre> section that contains valid FASTA.
### This should be the content of the 2nd <pre></pre> section but the hope
### is that this approach is a little bit more robust.
.scrape_IMGT_GENE_DB_result <- function(html)
{
    stopifnot(is.character(html))
    xml <- read_html(html)
    all_pre_elts <- html_text(html_elements(xml, "pre"))
    for (pre_elt in all_pre_elts) {
        fasta_lines <- strsplit(pre_elt, split="\n", fixed=TRUE)[[1L]]
        fasta_lines <- fasta_lines[nzchar(fasta_lines)]
        if (length(fasta_lines) == 0L)
            stop(wmsg("IMGT/GENE-DB returned 0 sequences"))
        if (.is_dna_fasta(fasta_lines))
            return(fasta_lines)
    }
    stop(wmsg("failed to scrape the results returned by IMGT/GENE-DB"))
}

.fetch_C_sequence_set_from_IMGT <-
    function(species, seqset_internal_nb, group=c("IGHC", "IGKC", "IGLC"))
{
    group <- match.arg(group)
    html <- .query_IMGT_GENE_DB(species, seqset_internal_nb, group)
    .scrape_IMGT_GENE_DB_result(html)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_C_sequence_set_from_IMGT()
###

### All kinds of conventions are used across the IMGT website to name
### organisms. I guess picking one and sticking to it would be kind of
### boring...
.map_organism_to_IMGT_species <- function(organism)
{
    stopifnot(isSingleNonWhiteString(organism))
    org <- tolower(organism)
    IMGT_species <- c(human="Homo sapiens",
                      mouse="Mus",
                      rat="Rattus norvegicus",
                      alpaca="Vicugna pacos",
                      rabbit="Oryctolagus cuniculus",
                      rhesus_monkey="Macaca mulatta")
    m <- match(chartr("_", " ", org), tolower(IMGT_species))
    if (!is.na(m))
        return(IMGT_species[[m]])
    m <- match(chartr(" ", "_", org), names(IMGT_species))
    if (!is.na(m))
        return(IMGT_species[[m]])

    shortname <- find_organism_shortname(org)
    ## find_organism_shortname() is guaranteed to return one of
    ## 'names(LATIN_NAMES)', and 'names(LATIN_NAMES)' is currently
    ## a subset of 'names(IMGT_species)' (see file LATIN_NAMES.R).
    ## So 'shortname' is guaranteed to be in 'names(IMGT_species)'.
    ## This means that 'IMGT_species[[shortname]]' should never fail.
    ## So why don't we just return that? Because this could change in
    ## the future e.g. if new entries get added to 'LATIN_NAMES'.
    ## Hence the extra work below.
    m <- match(shortname, names(IMGT_species))
    if (is.na(m))
        stop(wmsg("unrecognized organism: ", organism))
    IMGT_species[[m]]
}

### Fetch the C-region sequences from the links provided in the tables
### displayed at:
###
###   https://www.imgt.org/vquest/refseqh.html
###
### The C-region sequences are split in 3 IMGT groups: IGHC, IGKC, and IGLC.
### For each group, depending on the organism, 4 different sets of the
### C-region sequences are provided:
###
###   set | description                                         | note
###   --- | --------------------------------------------------- | ----
###    #1 | C-GENE exons: F+ORF+all P                           | (a)
###    #2 | C-GENE exons: F+ORF+in-frame P                      | (b)
###    #3 | C-GENE exons: F+ORF+in-frame P with IMGT gaps       | (c)
###    #4 | C-GENE artificially spliced exons: F+ORF+in-frame P | (d)
###
###   (a) The exon sequences in set #1 can contain N's.
###   (b) Set #2 is a subset of set #1. TO BE CONFIRMED: Sequences in this
###       set tend to be "cleaner" i.e. they have no N's (confirmed for Human,
###       still to be confirmed for other organisms). Note that, for Rhesus
###       monkey, one IGHC sequence in set #2 has a Y.
###   (c) Same exon sequences as in set #2 but with IMGT gaps. Note that
###       removing the gaps produces exactly the same sequences as in set #2.
###   (d) Set #4 is only available for a very limited number of organisms:
###       Human, Mouse, Rat, Alpaca, and Rabbit (as of Aug 19, 2025).
###       In particular rhesus monkey (Macaca mulatta) is missing.
###       TO BE CONFIRMED: It seems that the artificial gene sequences in
###       this set can be constructed by concatenating the corresponding
###       exon sequences from set #2, possibly after dropping the last exon.
###
### Sequence sets #1-3 correspond to the yellow columns in the big table of
### the 'IG "V-REGION", "D-REGION", "J-REGION", "C-GENE exon" sets' section.
### Sequence set #4 corresponds to the yellow column in the left table of
### the 'Constant gene artificially spliced exons sets' section (located
### at the bottom of the page).
### Note that:
### - Not all sequence sets are available for all organisms.
### - The link to a given set is a query to IMGT/GENE-DB.
###   See .query_IMGT_GENE_DB() above in this file for more information.
### - The 'seqset_nb' argument below must be an integer between 1 and 4
###   that specifies which set to fetch.
### - As of Aug 19, 2025, the mapping between 'seqset_nb'
###   and 'seqset_internal_nb' is as follow:
###
###     seqset_nb | seqset_internal_nb
###     --------- | ------------------
###             1 |              "7.2"
###             2 |              "7.5"
###             3 |              "7.1"
###             4 |             "14.1"
###
### 'seqset_nb' can also be a "sequence set internal number".
download_C_sequence_set_from_IMGT <-
    function(organism, destfile, group=c("IGHC", "IGKC", "IGLC"),
             seqset_nb=1L)
{
    species <- .map_organism_to_IMGT_species(organism)
    group <- match.arg(group)
    if (isSingleNonWhiteString(seqset_nb)) {
        seqset_internal_nb <- seqset_nb
    } else {
        stopifnot(isSingleInteger(seqset_nb), seqset_nb %in% 1:4)
        seqset_internal_nb <- .SEQSET_SET_INTERNAL_NUMBERS[[seqset_nb]]
    }
    sequences <- .fetch_C_sequence_set_from_IMGT(species, seqset_internal_nb,
                                                 group=group)
    writeLines(sequences, destfile)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_C_sequence_sets_from_IMGT()
###

### Reflects the C-region sequence sets available in IMGT/GENE-DB for
### the 5 official IgBLAST organisms as of Aug 21, 2025.
### See file LATIN_NAMES.R for more information.
.IMGT_ORGANISM_TO_C_SEQUENCE_SETS <- list(
    human=list(
        `7.2`=c("IGHC", "IGKC", "IGLC"),
        `7.5`=c("IGHC", "IGKC", "IGLC"),
        `7.1`=c("IGHC", "IGKC", "IGLC"),
       `14.1`=c("IGHC", "IGKC", "IGLC")
    ),
    mouse=list(
        `7.2`=c("IGHC", "IGKC", "IGLC"),
        `7.5`=c("IGHC", "IGKC", "IGLC"),
        `7.1`=c("IGHC", "IGKC", "IGLC"),
       `14.1`=c("IGHC")
    ),
    rabbit=list(
        `7.2`=c("IGHC", "IGKC", "IGLC"),
        `7.5`=c("IGHC", "IGKC", "IGLC"),
        `7.1`=c("IGHC", "IGKC", "IGLC"),
       `14.1`=c("IGHC")
    ),
    rat=list(
        `7.2`=c("IGHC", "IGKC", "IGLC"),
        `7.5`=c("IGHC", "IGKC", "IGLC"),
        `7.1`=c("IGHC", "IGKC", "IGLC"),
       `14.1`=c("IGHC")
    ),
    rhesus_monkey=list(
        `7.2`=c("IGHC", "IGKC", "IGLC"),
        `7.5`=c("IGHC", "IGKC", "IGLC"),
        `7.1`=c("IGHC", "IGKC", "IGLC")
    )
)

.normarg_IMGT_organisms <- function(organisms)
{
    supported_organisms <- names(.IMGT_ORGANISM_TO_C_SEQUENCE_SETS)
    if (is.null(organisms))
        return(supported_organisms)
    if (!is.character(organisms))
        stop(wmsg("'organisms' must be NULL or a character vector"))
    if (!all(organisms %in% supported_organisms)) {
        in1string <- paste(supported_organisms, collapse=", ")
        stop(wmsg("'organisms' must be a subset of: ", in1string))
    }
    if (anyDuplicated(organisms))
        stop(wmsg("'organisms' cannot contain duplicates"))
    organisms
}

### Use this to (re-)populate the igblastr/inst/extdata/constant_regions/IMGT/
### folder. The function must be called from within the folder.
### To fully (re-)populate it:
###
###     igblastr:::download_C_sequence_sets_from_IMGT()
###
### To (re-)populate only for a given organism:
###
###     igblastr:::download_C_sequence_sets_from_IMGT("rhesus_monkey")
###
### 'organisms' should be NULL or a character vector of organism names for
### which to download the sequence sets. If set to NULL, then the sequence
### sets for all the organisms listed in .IMGT_ORGANISM_TO_C_SEQUENCE_SETS
### get downloaded.
download_C_sequence_sets_from_IMGT <- function(organisms=NULL)
{
    organisms <- .normarg_IMGT_organisms(organisms)
    for (organism in organisms) {
        sequence_sets <- .IMGT_ORGANISM_TO_C_SEQUENCE_SETS[[organism]]
        for (seqset_internal_nb in names(sequence_sets)) {
            destdir <- file.path(organism, seqset_internal_nb)
            groups <- sequence_sets[[seqset_internal_nb]]
            for (group in groups) {
                filename <- paste0(group, ".fasta")
                destfile <- file.path(destdir, filename)
                seqset_label <- paste0(seqset_internal_nb, "/", group)
                message("Download sequence set ", seqset_label, " ",
                        "for ", organism, " ",
                        "to ", destfile, " ... ", appendLF=FALSE)
                download_C_sequence_set_from_IMGT(organism, destfile,
                                                  group=group,
                                                  seqset_nb=seqset_internal_nb)
                message("ok")
                nregion <- length(readDNAStringSet(destfile))
                message("  (", nregion, " region(s) downloaded)")
            }
        }
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normalize_IMGT_organism()
### form_IMGT_germline_db_name()

normalize_IMGT_organism <- function(organism)
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    chartr(" ", "_", organism)
}

form_IMGT_germline_db_name <- function(release, organism="Homo sapiens")
{
    if (!isSingleNonWhiteString(release))
        stop(wmsg("'relesase' must be a single (non-empty) string"))
    organism <- normalize_IMGT_organism(organism)
    sprintf("IMGT-%s.%s.%s", release, organism, "IGH+IGK+IGL")
}

