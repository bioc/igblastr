### =========================================================================
### Low-level utilities to retrieve C regions from IMGT
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level utilities to query IMGT/GENE-DB
###

### IMGT/GENE-DB Query page.
.IMGT_GENE_DB_URL <- "https://www.imgt.org/genedb/"

### IMPORTANT NOTE: Used to map 'seqset_nb' to 'seqset_internal_nb' so
### order is important! See .download_C_sequence_set_from_IMGT() below
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
### See .download_C_sequence_set_from_IMGT() below for more information
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .download_C_sequence_set_from_IMGT()
###

### The workhorse behind .download_C_sequence_set_from_IMGT().
.fetch_C_sequence_set_from_IMGT <-
    function(species, seqset_internal_nb, group=c("IGHC", "IGKC", "IGLC"))
{
    group <- match.arg(group)
    html <- .query_IMGT_GENE_DB(species, seqset_internal_nb, group)
    .scrape_IMGT_GENE_DB_result(html)
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
.download_C_sequence_set_from_IMGT <-
    function(organism, destfile, group=c("IGHC", "IGKC", "IGLC"),
             seqset_nb=1L)
{
    species <- find_organism_latin_name(organism)
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
        `7.1`=c("IGHC", "IGKC", "IGLC"),
       `14.1`=c("IGHC", "IGKC", "IGLC")
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
                .download_C_sequence_set_from_IMGT(organism, destfile,
                                                   group=group,
                                                   seqset_nb=seqset_internal_nb)
                message("ok")
                nregion <- length(readDNAStringSet(destfile))
                message("  (", nregion, " region(s) downloaded)")
            }
        }
    }
}

