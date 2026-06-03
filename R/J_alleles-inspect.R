### =========================================================================
### Basic inspection of J allele sequences
### -------------------------------------------------------------------------
###


.get_auxdata_col <- function(auxdata, colname)
{
    what <- c("a data.frame as returned by load_auxdata(), ",
              "compute_auxdata(), or compute_germline_db_auxdata()")
    if (!is.data.frame(auxdata))
        stop(wmsg("'auxdata' must be ", what))
    if (!isSingleNonWhiteString(colname))
        stop(wmsg("'colname' must be a single (non-empty) string"))
    auxdata_col <- auxdata[[colname]]
    if (is.null(auxdata_col))
        stop(wmsg("'auxdata' has no \"", colname, "\" column. ",
                  "Make sure that it's ", what, "."))
    auxdata_col
}

### Extracts the specified column from the 'auxdata' data.frame, and
### subset/reorder it to keep only the column values that correspond
### to the alleles in 'J_alleles'. Returns them in a named vector that
### is parallel to 'J_alleles' and has the allele names on it.
### The returned vector will have NAs for alleles that are not annotated
### in 'auxdata' or when 'auxdata[[colname]]' reports an NA for the allele.
.query_auxdata <- function(auxdata, J_alleles, colname, no.NAs=FALSE)
{
    stopifnot(isTRUEorFALSE(no.NAs))
    if (!is(J_alleles, "DNAStringSet"))
        stop(wmsg("'J_alleles' must be DNAStringSet object"))
    J_names <- names(J_alleles)
    if (is.null(J_names))
        stop(wmsg("'J_alleles' must have names"))
    auxdata_allele_names <- .get_auxdata_col(auxdata, "allele_name")
    m <- match(J_names, auxdata_allele_names)
    if (no.NAs) {
        bad_idx <- which(is.na(m))
        if (length(bad_idx) != 0L) {
            in1string <- paste(J_names[bad_idx], collapse=", ")
            msg <- c("the following alleles don't have ",
                     "an entry in 'auxdata': ", in1string)
            stop(wmsg(msg))
        }
    }
    auxdata_col <- .get_auxdata_col(auxdata, colname)
    ans <- auxdata_col[m]
    if (no.NAs) {
        bad_idx <- which(is.na(ans))
        if (length(bad_idx) != 0L) {
            in1string <- paste(J_names[bad_idx], collapse=", ")
            msg <- c("the \"", colname, "\" column in 'auxdata' is ",
                     "set to NA for the following alleles: ", in1string)
            stop(wmsg(msg))
        }
    }
    setNames(ans, J_names)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### print_J_alleles()
###

.get_aa_parts <- function(J_aa, cdr3_coding_frame_width)
{
    cdr3_ncodon <- cdr3_coding_frame_width %/% 3L
    cdr3 <- subseq(J_aa, end=cdr3_ncodon)
    fwr4 <- subseq(J_aa, start=cdr3_ncodon+1L)
    DataFrame(cdr3=cdr3, fwr4=fwr4)
}

.get_dna_parts <- function(J_dna, cdr3_end)
{
    cdr3 <- subseq(J_dna, end=cdr3_end)
    fwr4 <- subseq(J_dna, start=cdr3_end+1L)
    DataFrame(cdr3=cdr3, fwr4=fwr4)
}

.extract_cdr3_fwr4_parts_from_J_alleles <-
    function(J_alleles, auxdata, translate=FALSE, igblast_organism=NA)
{
    if (!isTRUEorFALSE(translate))
        stop(wmsg("'translate' must be TRUE or FALSE"))
    coding_frame_start <-
            .query_auxdata(auxdata, J_alleles, "coding_frame_start",
                           no.NAs=TRUE)
    cdr3_end <- .query_auxdata(auxdata, J_alleles, "cdr3_end")  # 0-based
    splitidx <- which(!is.na(cdr3_end))
    cdr3_end <- cdr3_end + 1L  # 1-based
    cdr3_coding_frame_width <- cdr3_end - coding_frame_start
    stopifnot(all(cdr3_coding_frame_width %% 3L == 0L, na.rm=TRUE))
    if (translate) {
        full_seq <- translate_codons(J_alleles, offset=coding_frame_start)
        aa0 <- rep.int(AAStringSet(""), length(full_seq))
        parts <- DataFrame(cdr3=aa0, fwr4=aa0)
        parts[splitidx, ] <- .get_aa_parts(full_seq[splitidx],
                                           cdr3_coding_frame_width[splitidx])
    } else {
        full_seq <- J_alleles
        ## Drop 'full_seq' metadata cols so they don't do strange things
        ## when we pass 'full_seq' to the DataFrame() call below.
        mcols(full_seq) <- NULL
        dna0 <- rep.int(DNAStringSet(""), length(full_seq))
        parts <- DataFrame(cdr3=dna0, fwr4=dna0)
        parts[splitidx, ] <- .get_dna_parts(full_seq[splitidx],
                                            cdr3_end[splitidx])
    }
    ans <- cbind(DataFrame(allele_name=names(J_alleles)),
                 parts,
                 DataFrame(full_seq=full_seq))
    if (!identical(igblast_organism, NA)) {
        if (!isSingleNonWhiteString(igblast_organism))
            stop(wmsg("'igblast_organism' must be NA or ",
                      "a single (non-empty) string"))
        igblast_organism <- normalize_igblast_organism(igblast_organism)
        igblast_auxdata <- load_auxdata(igblast_organism)
        ans$also_in_IgBLAST_auxdata <-
                    names(J_alleles) %in% igblast_auxdata[ , "allele_name"]
    }
    extra_cols <- mcols(J_alleles, use.names=FALSE)
    if (length(extra_cols) != 0L)
        ans <- cbind(ans, extra_cols)
    ans
}

.format_extra_cols <- function(extra_cols)
{
    stopifnot(is(extra_cols, "DataFrame"))
    fcolnames <- vapply(seq_along(extra_cols),
        function(j) {
            colname <- colnames(extra_cols)[[j]]
            colname_nc <- nchar(colname)
            col <- extra_cols[[j]]
            fcol <- format(col)
            ## We count on 'format(col)' to return a character vector parallel
            ## to 'col' but unfortunately it won't always do that. For example,
            ## if 'col' is an S4 object then there's a high chance that it will
            ## return a single string. For example, if 'col' is a DNAStringSet
            ## object then 'format(col)' will return:
            ##  "<S4 class ‘DNAStringSet’ [package “Biostrings”] with 5 slots>"
            ## Note that the check below can produce false negatives e.g. if
            ## it will evaluate to FALSE if 'col' is a DNAStringSet object of
            ## length 1!
            if (!is.character(fcol) || length(fcol) != length(col))
                stop(wmsg("format() failed to return a character vector ",
                          "of the right length on metadata column \"",
                          colname, "\""))
            colW <- nchar(fcol[[1L]])
            if (colW <= colname_nc)
                return(colname)
            justify <- if (is.character(col)) "left" else "right"
            format(colname, justify=justify, width=colW)
        }, character(1))
    fcols <- lapply(seq_along(extra_cols),
        function(j)
            format(extra_cols[[j]], width=nchar(fcolnames[[j]])))
    setNames(fcols, fcolnames)
}

.make_J_allele_header_line <-
    function(labelW, cdr3W, cdr3fwr4_sep, fwr4W, extra_colnames)
{
    stopifnot(isSingleInteger(labelW),
              isSingleInteger(cdr3W),
              isSingleString(cdr3fwr4_sep),
              isSingleInteger(fwr4W),
              is.character(extra_colnames))
    margin <- strrep(" ", labelW)
    if (cdr3W >= 5L) {
        cdr3 <- paste0(strrep(" ", cdr3W-5L), "CDR3>")
    } else if (cdr3W == 4L) {
        cdr3 <- "CDR3"
    } else {
        cdr3 <- strrep(" ", cdr3W)
    }
    if (fwr4W >= 5L) {
        fwr4 <- paste0("<FWR4", strrep(" ", fwr4W-5L))
    } else if (fwr4W == 4L) {
        fwr4 <- "FWR4"
    } else {
        fwr4 <- strrep(" ", fwr4W)
    }
    xcolnames <- paste(extra_colnames, collapse=" ")
    paste0(margin, cdr3, cdr3fwr4_sep, fwr4, xcolnames)
}

.add_colors <- function(x)
{
    stopifnot(is(x, "DNAString") || is(x, "AAString"))
    s <- as.character(x)
    class(s) <- c(seqtype(x), class(s))
    if (is(x, "DNAString"))
        return(Biostrings:::add_colors.DNA(s))
    Biostrings:::add_colors.AA(s)
}

.make_J_allele_line <- function(i, labels, Lfillers, cdr3,
                                   cdr3fwr4_sep, fwr4, Rfillers,
                                   full_seq, extra_cols)
{
    stopifnot(isSingleInteger(i),
              is.character(labels),
              is.character(Lfillers),
              is(cdr3, "XStringSet"),
              isSingleString(cdr3fwr4_sep),
              is(fwr4, "XStringSet"),
              is.character(Rfillers),
              length(labels) == length(Lfillers),
              length(labels) == length(cdr3),
              length(labels) == length(fwr4),
              length(labels) == length(Rfillers),
              is.list(extra_cols))
    seq1 <- .add_colors(cdr3[[i]])
    seq2 <- .add_colors(fwr4[[i]])
    fseq <- full_seq[[i]]
    ## Both 'Ltag' and 'Rtag' must be single characters. If that
    ## needs to change then the call to .make_J_allele_header_line()
    ## in .print_cdr3_fwr4_parts() below will need to be adjusted.
    if (length(fseq) == 0L) {
        Ltag <- Rtag <- " "
        stuff <- c(Lfillers[[i]], seq1, cdr3fwr4_sep, seq2, Rfillers[[i]])
    } else {
        Ltag <- Rtag <- "?"
        Linnertag <- ">"
        Rinnertag <- " <"  # note the space before the <
        LinnertagW <- nchar(Lfillers[[i]]) + nchar(cdr3fwr4_sep) +
                      nchar(Rfillers[[i]]) - nchar(fseq) - nchar(Rinnertag)
        Linnertag <- format(Linnertag, width=LinnertagW)
        stuff <- c(Linnertag, .add_colors(fseq), Rinnertag)
    }
    stuff <- paste(stuff, collapse="")
    xcols <- vapply(extra_cols, function(col) col[[i]], character(1))
    paste0(labels[[i]], Ltag, stuff, Rtag, paste(xcols, collapse=" "))
}

.print_cdr3_fwr4_parts <- function(cdr3_fwr4_parts,
                                   filler=".", cdr3fwr4_sep=" ")
{
    stopifnot(is(cdr3_fwr4_parts, "DataFrame"))
    allele_names <- cdr3_fwr4_parts[ , "allele_name"]
    cdr3         <- cdr3_fwr4_parts[ , "cdr3"]
    fwr4         <- cdr3_fwr4_parts[ , "fwr4"]
    full_seq     <- cdr3_fwr4_parts[ , "full_seq"]
    splitidx <- which(width(cdr3) != 0L | width(fwr4) != 0L)
    full_seq[splitidx] <- ""
    core_colnames <- c("allele_name", "cdr3", "fwr4", "full_seq")
    core_col_idx <- match(core_colnames, colnames(cdr3_fwr4_parts))
    extra_cols <- cdr3_fwr4_parts[-core_col_idx]

    ## Format all columns.
    labels <- paste0(format(seq_len(nrow(cdr3_fwr4_parts))), ". ",
                     format(paste0(allele_names, ": ")))
    labelW <- nchar(labels[[1L]])
    stopifnot(all(nchar(labels) == labelW))
    cdr3_maxwidth <- max(width(cdr3))
    Lfillers <- strrep(filler, cdr3_maxwidth - width(cdr3))
    fwr4_maxwidth <- max(width(fwr4))
    Rfillers <- strrep(filler, fwr4_maxwidth - width(fwr4))
    extra_cols <- .format_extra_cols(extra_cols)

    ## Print everything.

    locus <- substr(allele_names, 1L, 3L)
    group_lens <- runLength(Rle(locus))
    if (length(group_lens) == length(unique(locus))) {
        groupend_idx <- cumsum(head(group_lens, n=-1L))
    } else {
        groupend_idx <- integer(0)
    }

    ## We add 1 to 'cdr3_maxwidth' and 'fwr4_maxwidth' to account for
    ## the size of the 'Ltag' and 'Rtag' used in .make_J_allele_line().
    ## See .make_J_allele_line() above in this file.
    line <- .make_J_allele_header_line(labelW, cdr3_maxwidth+1L,
                                       cdr3fwr4_sep, fwr4_maxwidth+1L,
                                       names(extra_cols))
    stopifnot(is.character(line), length(line) == 1L)
    message(line)
    for (i in seq_len(nrow(cdr3_fwr4_parts))) {
        line <- .make_J_allele_line(i, labels, Lfillers, cdr3,
                                       cdr3fwr4_sep, fwr4, Rfillers,
                                       full_seq, extra_cols)
        stopifnot(is.character(line), length(line) == 1L)
        message(line)
        if (i %in% groupend_idx)
            message("")
    }
}

### Displays the J alleles sequences justified with respect to their
### CDR3/FWR4 junction.
print_J_alleles <- function(J_alleles, auxdata, translate=FALSE,
                            igblast_organism=NA, filler=".", cdr3fwr4_sep=" ")
{
    if (!is(J_alleles, "DNAStringSet")) {
        if (!isSingleNonWhiteString(J_alleles))
            stop(wmsg("'J_alleles' must be a DNAStringSet object ",
                      "containing germline J gene allele sequences, ",
                      "or a single string that is the name of a cached ",
                      "germline db"))
        db_name <- J_alleles
        J_alleles <- load_germline_sequences(db_name, region_types="J")
        if (!missing(auxdata))
            stop(wmsg("'auxdata' should not be supplied when 'J_alleles' ",
                      "is the name of a cached germline db"))
        auxdata <- load_auxdata(db_name)
        if (identical(igblast_organism, NA)) {
            igblast_organism <- infer_igblast_organism_from_db_name(db_name)
            if (is.na(igblast_organism))
                igblast_organism <- NA  # replace NA_character_ with NA
        }
    }
    if (!isSingleString(filler) || nchar(filler) != 1L)
        stop(wmsg("'filler' must be a single character"))
    if (!isSingleString(cdr3fwr4_sep))
        stop(wmsg("'cdr3fwr4_sep' must be a single string"))
    cdr3_fwr4_parts <- .extract_cdr3_fwr4_parts_from_J_alleles(
                                     J_alleles, auxdata, translate=translate,
                                     igblast_organism=igblast_organism)
    .print_cdr3_fwr4_parts(cdr3_fwr4_parts,
                           filler=filler, cdr3fwr4_sep=cdr3fwr4_sep)
    invisible(cdr3_fwr4_parts)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### translate_J_alleles()
###

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### J_allele_has_stop_codon()
###

### Only needs access to the "coding_frame_start" column in 'auxdata'.
### Returns a named logical vector that is parallel to 'J_alleles' and has
### the allele names on it.
### The returned vector will contain an NA for any allele that is not
### annotated in 'auxdata' or for which 'auxdata$coding_frame_start' has an NA.
J_allele_has_stop_codon <- function(J_alleles, auxdata)
{
    aa <- translate_J_alleles(J_alleles, auxdata)
    ans <- setNames(grepl("*", aa, fixed=TRUE), names(aa))
    ans[is.na(aa)] <- NA
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### translate_fwr4()
###

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

