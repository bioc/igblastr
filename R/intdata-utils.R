### =========================================================================
### Access, manipulate, and generate IgBLAST internal data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_intdata_path()
###

get_intdata_path <- function(organism, for.aa=FALSE,
                             domain_system=c("imgt", "kabat"),
                             which=c("live", "original"))
{
    organism <- normalize_igblast_organism(organism)
    if (!isTRUEorFALSE(for.aa))
        stop(wmsg("'for.aa' must be TRUE or FALSE"))
    domain_system <- match.arg(domain_system)
    which <- match.arg(which)
    intdata_dir <- file.path(path_to_igdata(which), "internal_data", organism)
    if (!dir.exists(intdata_dir))
        stop(wmsg("no internal data found in ",
                  dirname(intdata_dir), " for ", organism))
    intdata_filename <- sprintf("%s.%s.%s", organism,
                                if (for.aa) "pdm" else "ndm", domain_system)
    intdata_path <- file.path(intdata_dir, intdata_filename)
    if (!file.exists(intdata_path))
        stop(wmsg("internal data file ", intdata_filename, " ",
                  "not found in ", intdata_dir))
    intdata_path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_intdata()
###

### Not the true colnames used in IgBLAST intdata files: ours are all
### lowercase and we've replaced spaces with underscores.
### Note that many columns are redundant:
### - columns 'cdr1_start', 'fwr2_start', 'cdr2_start', and 'fwr3_start'
###   are redundant with columns 'fwr1_end', 'cdr1_end', 'fwr2_end',
###   and 'cdr2_end', respectively;
### - columns 'fwr1_start' and 'coding_frame_start' are redundant.
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

### IMPORTANT NOTE: The FWR/CDR positions in the returned data.frame are
### 1-based while the coding frame start positions are 0-based!
load_intdata <- function(organism, for.aa=FALSE,
                         domain_system=c("imgt", "kabat"),
                         which=c("live", "original"))
{
    domain_system <- match.arg(domain_system)
    which <- match.arg(which)
    intdata_path <- get_intdata_path(organism, for.aa=for.aa,
                                     domain_system=domain_system,
                                     which=which)
    read_broken_table(intdata_path, .IGBLAST_INTDATA_COL2CLASS)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### translate_V_alleles()
### V_allele_has_stop_codon()
###

### Extracts the specified column from the 'indata' data.frame, and
### subset/reorder it to keep only the column values that correspond
### to the alleles in 'V_alleles'. Returns them in a named vector that
### is parallel to 'V_alleles' and has the allele names on it.
### The returned vector will have NAs for alleles that are not annotated
### in 'indata' or when 'indata[[colname]]' reports an NA for the allele.
.query_intdata <- function(intdata, V_alleles, colname)
{
    if (!is.data.frame(intdata))
        stop(wmsg("'intdata' must be a data.frame as returned ",
                  "by load_intdata()"))
    intdata_allele_name <- intdata$allele_name
    if (is.null(intdata_allele_name))
        stop(wmsg("'intdata' has no \"allele_name\" column. Make sure ",
                  "that it's a data.frame as returned by load_intdata() ",
                  "or compute_intdata()."))
    if (!is(V_alleles, "DNAStringSet"))
        stop(wmsg("'V_alleles' must be DNAStringSet object"))
    V_names <- names(V_alleles)
    if (is.null(V_names))
        stop(wmsg("'V_alleles' must have names"))
    if (!isSingleNonWhiteString(colname))
        stop(wmsg("'colname' must be a single (non-empty) string"))
    intdata_col <- intdata[[colname]]
    if (is.null(intdata_col))
        stop(wmsg("'intdata' has no \"", colname, "\" column. Make sure ",
                  "that it's a data.frame as returned by load_intdata()."))
    setNames(intdata_col[match(V_names, intdata_allele_name)], V_names)
}

.translate_V_codons <- function(V_alleles, offsets, with.init.codon)
{
    stopifnot(is(V_alleles, "DNAStringSet"), is.integer(offsets),
              length(V_alleles) == length(offsets))
    ans <- rep.int(NA_character_, length(V_alleles))
    selection_idx <- which(!is.na(offsets))
    if (length(selection_idx) != 0L) {
        dna <- V_alleles[selection_idx]
        off <- offsets[selection_idx]
        aa <- translate_codons(dna, offset=off, with.init.codon=with.init.codon)
        ans[selection_idx] <- as.character(aa)
    }
    setNames(ans, names(V_alleles))
}

### Translates the coding frame contained in the V allele sequence.
### Only needs access to the "coding_frame_start" column in 'intdata'.
### Returns the amino acid sequences in a named character vector that
### is parallel to 'V_alleles' and has the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'intdata' or for which 'intdata$coding_frame_start' has
### an NA.
.translate_V_coding_frame <- function(V_alleles, intdata)
{
    offsets <- .query_intdata(intdata, V_alleles, "coding_frame_start")
    .translate_V_codons(V_alleles, offsets, with.init.codon=TRUE)
}

### Only needs access to the "<region>_start" and "<region>_end" columns
### of the 'intdata' data.frame.
### Returns the amino acid sequences in a named character vector that
### is parallel to 'V_alleles' and has the allele names on it.
### The returned vector will contain an NA for any allele that is
### not annotated in 'intdata' or for which 'intdata$<region>_start'
### or 'intdata$<region>_end' has an NA.
.translate_V_region <- function(V_alleles, intdata, region)
{
    valid_regions <- c(paste0(c("fwr", "cdr"), rep(1:2, each=2)), "fwr3")
    if (!(isSingleNonWhiteString(region) && (region %in% valid_regions))) {
        in1string <- paste0("\"", valid_regions, "\"", collapse=", ")
        stop(wmsg("'region' must be one of ", in1string))
    }
    start_colname <- paste0(region, "_start")
    end_colname <- paste0(region, "_end")
    starts <- .query_intdata(intdata, V_alleles, start_colname)  # 1-based
    ends <- .query_intdata(intdata, V_alleles, end_colname)  # 1-based
    offsets <- starts - 1L
    with.init.codon <- region == "fwr1"
    ans <- .translate_V_codons(V_alleles, offsets, with.init.codon)
    ncodons <- (ends - offsets) %/% 3L
    substr(ans, 1L, ncodons)
}

translate_V_alleles <- function(V_alleles, intdata, region=NULL)
{
    if (is.null(region))
        return(.translate_V_coding_frame(V_alleles, intdata))
    .translate_V_region(V_alleles, intdata, region)
}

### Only needs access to the "coding_frame_start" column in 'intdata'.
### Returns a named logical vector that is parallel to 'V_alleles' and has
### the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'intdata' or for which 'intdata$coding_frame_start' has an NA.
V_allele_has_stop_codon <- function(V_alleles, intdata)
{
    V_aa <- translate_V_alleles(V_alleles, intdata)
    ans <- setNames(grepl("*", V_aa, fixed=TRUE), names(V_aa))
    ans[is.na(V_aa)] <- NA
    ans
}

