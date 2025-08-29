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
.download_C_sequence_set_from_IMGT <-
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### PROBLEM: The IMGT folks do not provide the 14.1 IGHC sequences (a.k.a.
### "Constant gene artificially spliced exons") for Rhesus monkey.
###
### QUESTION: Can the 14.1 IGHC sequences for Human be programatically
### inferred from the 7.5 IGHC sequences? If that turns out to be possible,
### then we should be able to infer the 14.1 IGHC sequences for Rhesus monkey
### from the 7.5 IGHC sequences.
###
### After a quick examination, this seems to be more complicated than it
### sounds. Because:
### 1. Only a small subset of the genes represented in 7.5/IGHC for Human
###    are found in 14.1/IGHC. So it seems that the IMGT folks are filtering
###    out some C-genes when they generate the 14.1/IGHC set for Human.
###    The question is: what criteria do they use for this filtering?
### 2. There seems to be some alternate splicing involved in the IGHC genes
###    for Human. This hypothesis is based on the following observations:
###    - The first two sequences in human/14.1/IGHC.fasta are both for gene
###      allele IGHA1*01.
###    - There are 4 exons associated with IGHA1*01 in human/7.5/IGHC.fasta
###      Splicing the first 3 exons produces exactly the first sequence in
###      human/14.1/IGHC.fasta.
###    - However, it's not clear how to splice the exons associated with
###      IGHA1*01 to produce the second sequence in human/14.1/IGHC.fasta.
###      The sequence seems to be the result of concatenating the first 2
###      exons plus something else but it's not clear what this something
###      else is (it's not the 3rd exons and it's not the 4th exon either).
###    So quite confusing but this seems to suggest that some sort of alternate
###    splicing is involved.
###
### So because of these complications, my attempt below to programatically
### infer the 14.1 IGHC sequences for Human from the 7.5 IGHC sequences didn't
### go very far. Before trying to go any further, I sent an email to the IMGT
### folks on Aug 27, 2025, to ask for help.

if (FALSE) {
  IGHC_7.5 <- readDNAStringSet("7.5/IGHC.fasta")
  IGHC_14.1 <- readDNAStringSet("14.1/IGHC.fasta")

  ### Returns the split FASTA headers in a 15-col matrix.
  .parse_IMGT_FASTA_headers <- function(dna)
  {
    stopifnot(is(dna, "DNAStringSet"))
    split_names <- strsplit(names(dna), "|", fixed=TRUE)
    matrix(unlist(split_names), ncol=15L, byrow=TRUE)
  }

  ### 'heads' must be the matrix returned by calling
  ### .parse_IMGT_FASTA_headers() on IGHC_14.1.
  .extract_ranges_of_exon_to_splice <- function(heads)
  {
    stopifnot(is.matrix(heads), ncol(heads) == 15L)
    heads[ , 6L]
  }

  IGHC_7.5_heads <- .parse_IMGT_FASTA_headers(IGHC_7.5)
  IGHC_14.1_heads <- .parse_IMGT_FASTA_headers(IGHC_14.1)

  ### The genes in IGHC_14.1 must be a subset of the genes in IGHC_7.5.
  stopifnot(all(IGHC_14.1_heads[ , 2L] %in% IGHC_7.5_heads[ , 2L]))
}

