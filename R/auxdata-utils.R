### =========================================================================
### Access or generate IgBLAST auxiliary data
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_auxdata_path()
###

get_auxdata_path <- function(organism, which=c("live", "original"))
{
    organism <- normalize_igblast_organism(organism)
    which <- match.arg(which)
    auxdir <- file.path(path_to_igdata(which), "optional_file")
    auxfile <- file.path(auxdir, paste0(organism, "_gl.aux"))
    if (!file.exists(auxfile))
        stop(wmsg("no auxiliary data found in ", auxdir, " for ", organism))
    auxfile
}

get_igblast_auxiliary_data <- function(...)
{
    .Deprecated("get_auxdata_path")
    get_auxdata_path(...)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_auxdata()
###

### 'x' must be a list of character vectors of variable length.
### Conceptually right-pads the list elements with empty strings ("")
### to make the list "constant-width" before unlisting it.
### Returns a character vector of length 'length(x) * width'.
.right_pad_with_empty_strings_and_unlist <- function(x, width=NA)
{
    stopifnot(is.list(x), isSingleNumberOrNA(width))
    x_len <- length(x)
    if (x_len == 0L)
        return(character(0))
    x_lens <- lengths(x)
    max_x_lens <- max(x_lens)
    if (is.na(width)) {
        width <- max_x_lens
    } else {
        width <- as.integer(width)
        stopifnot(width >= max_x_lens)
    }
    y_lens <- width - x_lens
    x_seqalong <- seq_along(x)
    f <- rep.int(x_seqalong, y_lens)
    attributes(f) <- list(levels=as.character(x_seqalong), class="factor")
    y <- split(character(length(f)), f)
    collate_subscript <- rep(x_seqalong, each=2L)
    collate_subscript[2L * x_seqalong] <- x_seqalong + x_len
    unlist(c(x, y)[collate_subscript], recursive=FALSE, use.names=FALSE)
}

### Returns the table in a character matrix.
.read_jagged_table_as_character_matrix <- function(file)
{
    lines <- readLines(file)
    lines <- lines[nzchar(lines) & !has_prefix(lines, "#")]
    data <- strsplit(lines, split="[ \t]+")
    data <- .right_pad_with_empty_strings_and_unlist(data)
    matrix(data, nrow=length(lines), byrow=TRUE)
}

.make_auxiliary_data_df_from_matrix <- function(m, filepath)
{
    if (ncol(m) != 5L)
        stop(wmsg("error loading ", filepath, ": unexpected number of fields"))
    data.frame(
        allele_name       =           m[ , 1L],
        coding_frame_start=as.integer(m[ , 2L]),
        chain_type        =           m[ , 3L],
        CDR3_stop         =as.integer(m[ , 4L]),
        extra_bps         =as.integer(m[ , 5L])
    )
}

### IgBLAST *.aux files are supposedly "tab-delimited" but they are
### broken in various ways:
###   - they use a mix of whitespaces for the field separators;
###   - each line contains a variable number of fields;
###   - some lines contain trailing whitespaces.
### So we cannot read them with read.table().
load_auxdata <- function(organism, which=c("live", "original"))
{
    which <- match.arg(which)
    filepath <- get_auxdata_path(organism, which)
    m <- .read_jagged_table_as_character_matrix(filepath)
    .make_auxiliary_data_df_from_matrix(m, filepath)
}

load_igblast_auxiliary_data <- function(...)
{
    .Deprecated("")
    load_auxdata(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compute_auxdata()
###

### The FWR4 region is expected to **always** start with AA pattern "WGXG"
### on the heavy chain.
.WGXG_pattern <- "TGGGGNNNNGGN"  # reverse-translation of "WGXG"

### Returns a named integer vector parallel to 'J_alleles' that contains
### the **0-based** FWR4 start position for each sequence in 'J_alleles'.
### 'FWR4_starts' will be set to NA for alleles that don't have a match.
### For alleles with more than one match, we keep the first match only.
### The names on the returned vector indicate the AA pattern that was used
### to determine the start of the FWR4 region.
.find_heavy_FWR4_starts <- function(J_alleles)
{
    stopifnot(is(J_alleles, "DNAStringSet"))
    m <- vmatchPattern(.WGXG_pattern, J_alleles, fixed=FALSE)
    ans <- as.integer(heads(start(m), n=1L)) - 1L
    names(ans) <- ifelse(is.na(ans), NA_character_, "WGXG")
    ans
}

### The FWR4 region is expected to **always** start with AA pattern "FGXG"
### on the light chain.
.FGXG_pattern <- "TTYGGNNNNGGN"  # reverse-translation of "FGXG"

### EXPERIMENTAL!
### The "FGXG" pattern is not found for 4 J alleles in
### IMGT-202531-1.Mus_musculus.IGH+IGK+IGL: IGKJ3*01, IGKJ3*02,
### IGLJ2P*01, IGLJ3P*01. However, except for IGLJ2P*01, these alleles
### are annotated in mouse_gl.aux with a CDR3 stop reported at position 6
### (0-based). Turns out that for the 3 alleles annotated in mouse_gl.aux,
### the two first codons of the FWR4 region translate to AA sequence "FS".
### Is this a coincidence or does the FS sequence actually play a role on
### the light chain? What do biologists say about this? In particular, does
### it make sense to use this alternative pattern to identify the start of
### the FWR4 region on the light chain when the "FGXG" pattern is not found?
### Note that all the possible reverse-translations of FS cannot be
### represented with a single DNA pattern (even with the use of IUPAC
### ambiguity codes).
.FS_pattern1 <- "TTYTCN"
.FS_pattern2 <- "TTYAGY"

### UPDATE on using the "FS" pattern to identify the start of the FWR4
### region on the light chain when the "FGXG" pattern is not found:
### Works well for IMGT-202531-1.Mus_musculus.IGH+IGK+IGL (well, it was
### specifically designed for that so no surprise here), but not
### so well for IMGT-202531-1.Rattus_norvegicus.IGH+IGK+IGL or
### IMGT-202531-1.Oryctolagus_cuniculus.IGH+IGK+IGL (rabbit)
### or IMGT-202531-1.Macaca_mulatta.IGH+IGK+IGL (rhesus monkey).
### So we disabled this feature in .find_light_FWR4_starts() below.

.find_light_FWR4_starts <- function(J_alleles)
{
    stopifnot(is(J_alleles, "DNAStringSet"))
    m <- vmatchPattern(.FGXG_pattern, J_alleles, fixed=FALSE)
    FGXG_starts <- as.integer(heads(start(m), n=1L))
    names(FGXG_starts) <- ifelse(is.na(FGXG_starts), NA_character_, "FGXG")
    ## Disabling search for alternative "FS" pattern for now.
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

.VALID_J_GROUPS <- paste0("IG", c("H", "K", "L"), "J")

### Returns a data.frame with the same column names as the data.frame
### returned by .make_auxiliary_data_df_from_matrix() above, plus
### the "FWR4_start_pattern" column.
### NOTE: We set coding_frame_start/CDR3_stop/extra_bps/FWR4_start_pattern
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
        FWR4_starts <- .find_heavy_FWR4_starts(J_alleles)
    } else {
        FWR4_starts <- .find_light_FWR4_starts(J_alleles)
    }
    coding_frame_starts <- FWR4_starts %% 3L
    extra_bps <- (width(J_alleles) - coding_frame_starts) %% 3L
    data.frame(
        allele_name       =allele_names,
        coding_frame_start=coding_frame_starts,
        chain_type        =chain_type,
        CDR3_stop         =FWR4_starts - 1L,  # 0-based
        extra_bps         =extra_bps
        ## Returning this column only made sense when we were using "FS"
        ## pattern as a 2nd-chance pattern on the light chain.
        #FWR4_start_pattern=names(FWR4_starts)
    )
}

.has_stop_codon <- function(J_alleles, auxdata)
{
    stopifnot(is(J_alleles, "DNAStringSet"), is.data.frame(auxdata),
              identical(names(J_alleles), auxdata[ , "allele_name"]))
    ans <- rep.int(NA, length(J_alleles))
    coding_frame_start <- auxdata[ , "coding_frame_start"]
    searchme_idx <- which(!is.na(coding_frame_start))
    if (length(searchme_idx) != 0L) {
        alleles <- J_alleles[searchme_idx]
        Ltrim <- coding_frame_start[searchme_idx]
        Rtrim <- (width(alleles) - Ltrim) %% 3L
        trimmed <- narrow(alleles, start=Ltrim+1L, end=-Rtrim-1L)
        alleles_aa <- translate(trimmed, no.init.codon=TRUE)
        ans[searchme_idx] <- vcountPattern("*", alleles_aa) != 0L
    }
    ans
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
    ans$has_stop_codon <- .has_stop_codon(J_alleles, ans)

    ## Warn user if CDR3 stop not found for some alleles.
    bad_idx <- which(is.na(ans[ , "CDR3_stop"]))
    if (length(bad_idx) != 0L) {
        in1string <- paste(ans[bad_idx, "allele_name"], collapse=", ")
        warning(wmsg("CDR3 stop not found for allele(s): ", in1string),
                "\n  ",
                wmsg("--> coding_frame_start/CDR3_stop/extra_bps ",
                     "were set to NA for these alleles"))
    }

    ## Warn user if stop codon found in some alleles.
    bad_idx <- which(ans[ , "has_stop_codon"])
    if (length(bad_idx) != 0L) {
        in1string <- paste(ans[bad_idx, "allele_name"], collapse=", ")
        warning(wmsg("stop codon found in allele(s): ", in1string))
    }

    ans
}

