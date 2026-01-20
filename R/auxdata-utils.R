### =========================================================================
### Access, manipulate, and generate IgBLAST auxiliary data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_auxdata_path()
###

get_auxdata_path <- function(organism, which=c("live", "original"))
{
    organism <- normalize_igblast_organism(organism)
    which <- match.arg(which)
    auxdata_dir <- file.path(path_to_igdata(which), "optional_file")
    auxdata_filename <- paste0(organism, "_gl.aux")
    auxdata_path <- file.path(auxdata_dir, auxdata_filename)
    if (!file.exists(auxdata_path))
        stop(wmsg("no auxiliary data found in ",
                  auxdata_dir, " for ", organism))
    auxdata_path
}

get_igblast_auxiliary_data <- function(...)
{
    .Deprecated("get_auxdata_path")
    get_auxdata_path(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_auxdata()
###

### Not the true colnames used in IgBLAST auxdata files. However, ours are
### shorter, all lowercase, and contain underscores instead of spaces.
.IGBLAST_AUXDATA_COL2CLASS <- c(
    allele_name="character",
    coding_frame_start="integer",
    chain_type="character",
    cdr3_end="integer",
    extra_bps="integer"
)

### IMPORTANT NOTE: Unlike with the data.frame returned by load_intdata(),
### all the positions in the data.frame returned by load_auxdata() (that is,
### the positions reported in columns 'coding_frame_start' and 'cdr3_end')
### are 0-based!
load_auxdata <- function(organism, which=c("live", "original"))
{
    which <- match.arg(which)
    auxdata_path <- get_auxdata_path(organism, which=which)
    read_broken_table(auxdata_path, .IGBLAST_AUXDATA_COL2CLASS)
}

load_igblast_auxiliary_data <- function(...)
{
    .Deprecated("load_auxdata")
    load_auxdata(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### translate_J_alleles()
### J_allele_has_stop_codon()
### translate_fwr4()
###

.get_auxdata_col <- function(auxdata, colname)
{
    if (!is.data.frame(auxdata))
        stop(wmsg("'auxdata' must be a data.frame as returned ",
                  "by load_auxdata() or compute_auxdata()"))
    if (!isSingleNonWhiteString(colname))
        stop(wmsg("'colname' must be a single (non-empty) string"))
    auxdata_col <- auxdata[[colname]]
    if (is.null(auxdata_col))
        stop(wmsg("'auxdata' has no \"", colname, "\" column. Make sure ",
                  "that it's a data.frame as returned by load_auxdata() ",
                  "or compute_auxdata()."))
    auxdata_col
}

### Extracts the specified column from the 'auxdata' data.frame, and
### subset/reorder it to keep only the column values that correspond
### to the alleles in 'J_alleles'. Returns them in a named vector that
### is parallel to 'J_alleles' and has the allele names on it.
### The returned vector will have NAs for alleles that are not annotated
### in 'auxdata' or when 'auxdata[[colname]]' reports an NA for the allele.
.query_auxdata <- function(auxdata, J_alleles, colname)
{
    allele_names <- .get_auxdata_col(auxdata, "allele_name")
    if (!is(J_alleles, "DNAStringSet"))
        stop(wmsg("'J_alleles' must be DNAStringSet object"))
    J_names <- names(J_alleles)
    if (is.null(J_names))
        stop(wmsg("'J_alleles' must have names"))
    auxdata_col <- .get_auxdata_col(auxdata, colname)
    setNames(auxdata_col[match(J_names, allele_names)], J_names)
}

### Translates the coding frame contained in the J allele sequence.
### Only needs access to the "coding_frame_start" column in 'auxdata'.
### Returns the amino acid sequences in a named character vector that
### is parallel to 'J_alleles' and has the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'auxdata' or for which 'auxdata$coding_frame_start' has an NA.
translate_J_alleles <- function(J_alleles, auxdata)
{
    coding_frame_start <- .query_auxdata(auxdata, J_alleles,
                                         "coding_frame_start")
    ans <- rep.int(NA_character_, length(J_alleles))
    selection_idx <- which(!is.na(coding_frame_start))
    if (length(selection_idx) != 0L) {
        dna <- J_alleles[selection_idx]
        offset <- coding_frame_start[selection_idx]
        aa <- translate_codons(dna, offset=offset)
        ans[selection_idx] <- as.character(aa)
    }
    setNames(ans, names(J_alleles))
}

### Only needs access to the "coding_frame_start" column in 'auxdata'.
### Returns a named logical vector that is parallel to 'J_alleles' and has
### the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'auxdata' or for which 'auxdata$coding_frame_start' has an NA.
J_allele_has_stop_codon <- function(J_alleles, auxdata)
{
    J_aa <- translate_J_alleles(J_alleles, auxdata)
    ans <- setNames(grepl("*", J_aa, fixed=TRUE), names(J_aa))
    ans[is.na(J_aa)] <- NA
    ans
}

### Only needs access to the "cdr3_end" column of the 'auxdata' data.frame.
### Returns the amino acid sequences in a named character vector that
### is parallel to 'J_alleles' and has the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'auxdata' or for which 'auxdata$cdr3_end' has an NA.
translate_fwr4 <- function(J_alleles, auxdata, max.codons=NA)
{
    if (!isSingleNumberOrNA(max.codons))
        stop(wmsg("'max.codons' must be a single number or NA"))
    if (!is.integer(max.codons))
        max.codons <- as.integer(max.codons)

    cdr3_end <- .query_auxdata(auxdata, J_alleles, "cdr3_end")  # 0-based
    ans <- rep.int(NA_character_, length(J_alleles))
    selection_idx <- which(!is.na(cdr3_end))
    if (length(selection_idx) != 0L) {
        dna <- J_alleles[selection_idx]
        offset <- cdr3_end[selection_idx] + 1L  # 0-based FWR4 start
        aa <- translate_codons(dna, offset=offset)
        ans[selection_idx] <- as.character(aa)
    }
    if (!is.na(max.codons))
        ans <- substr(ans, 1L, max.codons)
    setNames(ans, names(J_alleles))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .find_heavy_fwr4_starts()
### .find_light_fwr4_starts()
###

### For alleles in the IGHJ group (i.e. BCR germline J gene alleles on the
### heavy chain), the FWR4 region is expected to start with AA motif "WGXG.
.WGXG_pattern <- "TGGGGNNNNGGN"  # reverse-translation of "WGXG"

### For all other J alleles, that is, for alleles in the IG[KL]J groups
### (i.e. BCR germline J gene alleles on the light chain) and all TCR
### germline J gene alleles, the FWR4 region is expected to start with
### AA motif "FGXG".
.FGXG_pattern <- "TTYGGNNNNGGN"  # reverse-translation of "FGXG"

### EXPERIMENTAL!
### The "FGXG" motif is not found for 4 J alleles in
### IMGT-202531-1.Mus_musculus.IGH+IGK+IGL: IGKJ3*01, IGKJ3*02,
### IGLJ2P*01, IGLJ3P*01. However, except for IGLJ2P*01, these alleles
### are annotated in mouse_gl.aux with a CDR3 end reported at position 6
### (0-based). Turns out that for the 3 alleles annotated in mouse_gl.aux,
### the two first codons of the FWR4 region translate to AA sequence "FS".
### Is this a coincidence or does the FS sequence actually play a role on
### the light chain? What do biologists say about this? In particular, does
### it make sense to use this alternative motif to identify the start of
### the FWR4 region on the light chain when the "FGXG" motif is not found?
### Note that all the possible reverse-translations of FS cannot be
### represented with a single DNA pattern (even with the use of IUPAC
### ambiguity codes).
.FS_pattern1 <- "TTYTCN"
.FS_pattern2 <- "TTYAGY"

### UPDATE on using the "FS" motif to identify the start of the FWR4
### region on the light chain when the "FGXG" motif is not found:
### Works well for IMGT-202531-1.Mus_musculus.IGH+IGK+IGL (well, it was
### specifically designed for that so no surprise here), but not
### so well for IMGT-202531-1.Rattus_norvegicus.IGH+IGK+IGL or
### IMGT-202531-1.Oryctolagus_cuniculus.IGH+IGK+IGL (rabbit)
### or IMGT-202531-1.Macaca_mulatta.IGH+IGK+IGL (rhesus monkey).
### So we disabled this feature in .find_light_fwr4_starts() below.

### .find_heavy_fwr4_starts() and .find_light_fwr4_starts() both return
### a named integer vector parallel to 'J_alleles' that contains
### the **0-based** FWR4 start position for each sequence in 'J_alleles'.
### Th FWR4 start will be set to NA for alleles that don't have a match.
### For alleles with more than one match, we keep the first match only.
### The names on the returned vector indicate the AA motif that was used
### to determine the start of the FWR4 region.

.find_heavy_fwr4_starts <- function(J_alleles)
{
    stopifnot(is(J_alleles, "DNAStringSet"))
    m <- vmatchPattern(.WGXG_pattern, J_alleles, fixed=FALSE)
    ans <- as.integer(heads(start(m), n=1L)) - 1L
    names(ans) <- ifelse(is.na(ans), NA_character_, "WGXG")
    ans
}

.find_light_fwr4_starts <- function(J_alleles)
{
    stopifnot(is(J_alleles, "DNAStringSet"))
    m <- vmatchPattern(.FGXG_pattern, J_alleles, fixed=FALSE)
    FGXG_starts <- as.integer(heads(start(m), n=1L))
    names(FGXG_starts) <- ifelse(is.na(FGXG_starts), NA_character_, "FGXG")
    ## Disabling search for alternative "FS" motif for now.
    #na_idx <- which(is.na(FGXG_starts))
    #if (length(na_idx) != 0L) {
    #    dangling_alleles <- J_alleles[na_idx]
    #    m <- vmatchPattern(.FS_pattern1, dangling_alleles, fixed=FALSE)
    #    FS_starts1 <- as.integer(heads(start(m), n=1L))
    #    m <- vmatchPattern(.FS_pattern2, dangling_alleles, fixed=FALSE)
    #    FS_starts2 <- as.integer(heads(start(m), n=1L))
    #    FS_starts <- pmin(FS_starts1, FS_starts2, na.rm=TRUE)
    #    FGXG_starts[na_idx] <- FS_starts
    #    names(FGXG_starts)[na_idx] <-
    #        ifelse(is.na(FS_starts), NA_character_, "FS")
    #}
    FGXG_starts - 1L
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compute_auxdata()
###

.VALID_J_GROUPS <- paste0("IG", c("H", "K", "L"), "J")

### Returns a data.frame with the same column names as the data.frame
### returned by load_auxdata() above, plus the "fwr4_start_motif" column.
### NOTE: We set coding_frame_start/cdr3_end/extra_bps/fwr4_start_motif
### to NA for alleles for which the FWR4 start cannot be determined.
.compute_auxdata_for_J_group <- function(J_alleles, J_group)
{
    stopifnot(is(J_alleles, "DNAStringSet"), J_group %in% .VALID_J_GROUPS)
    allele_names <- names(J_alleles)
    stopifnot(!is.null(allele_names))
    if (length(J_alleles) == 0L) {
        chain_type <- character(0)
    } else {
        allele_groups <- substr(allele_names, 1L, 4L)
        stopifnot(all(allele_groups == J_group))
        chain_type <- paste0("J", substr(J_group, 3L, 3L))
    }
    if (J_group == "IGHJ") {
        fwr4_starts <- .find_heavy_fwr4_starts(J_alleles)
    } else {
        fwr4_starts <- .find_light_fwr4_starts(J_alleles)
    }
    coding_frame_starts <- fwr4_starts %% 3L
    extra_bps <- (width(J_alleles) - coding_frame_starts) %% 3L
    data.frame(
        allele_name       =allele_names,
        coding_frame_start=coding_frame_starts,  # 0-based
        chain_type        =chain_type,
        cdr3_end          =fwr4_starts - 1L,     # 0-based
        extra_bps         =extra_bps
        ## Returning this column only made sense when we were using "FS"
        ## motif as a 2nd-chance motif on the light chain.
        #fwr4_start_motif  =names(fwr4_starts)
    )
}

### Returns a data.frame with 1 row per sequence in 'J_alleles'.
compute_auxdata <- function(J_alleles)
{
    if (!is(J_alleles, "DNAStringSet"))
        stop(wmsg("'J_alleles' must be DNAStringSet object"))
    allele_names <- names(J_alleles)
    if (is.null(allele_names))
        stop(wmsg("'J_alleles' must have names"))
    allele_groups <- substr(allele_names, 1L, 4L)
    if (!all(allele_groups %in% .VALID_J_GROUPS))
        stop(wmsg("all allele names must start with 'IG[HKL]J'"))

    JH_alleles <- J_alleles[allele_groups == "IGHJ"]
    JK_alleles <- J_alleles[allele_groups == "IGKJ"]
    JL_alleles <- J_alleles[allele_groups == "IGLJ"]
    JH_df <- .compute_auxdata_for_J_group(JH_alleles, "IGHJ")
    JK_df <- .compute_auxdata_for_J_group(JK_alleles, "IGKJ")
    JL_df <- .compute_auxdata_for_J_group(JL_alleles, "IGLJ")
    ans <- rbind(JH_df, JK_df, JL_df)

    i <- match(allele_names, ans[ , "allele_name"])
    ans <- S4Vectors:::extract_data_frame_rows(ans, i)
    rownames(ans) <- NULL

    ## Warn user if CDR3 end not found for some alleles.
    bad_idx <- which(is.na(ans[ , "cdr3_end"]))
    if (length(bad_idx) != 0L) {
        in1string <- paste(ans[bad_idx, "allele_name"], collapse=", ")
        warning(wmsg("CDR3 end not found for allele(s): ", in1string),
                "\n  ",
                wmsg("--> coding_frame_start, cdr3_end, and extra_bps ",
                     "were set to NA for these alleles"))
    }

    ans
}

