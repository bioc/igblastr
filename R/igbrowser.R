### =========================================================================
### igbrowser()
### -------------------------------------------------------------------------


.NUCLEOTIDE_COLORS <- c(A="rgb(240, 128, 128)",
                        C="rgb(128, 192, 100)",
                        G="rgb(128, 176, 225)",
                        T="rgb(192, 176, 96)")

.make_style_attrib <- function(attribs)
{
    stopifnot(is.character(attribs))
    paste0("style=\"", paste(attribs, collapse="; "), "\"")
}

.set_background_color <- function(region, color)
{
    stopifnot(is.character(region), isSingleString(color))
    attribs <- c(sprintf("background: %s", color), "font-weight: bold")
    span_tag <- sprintf("<span %s>", .make_style_attrib(attribs))
    paste0(span_tag, region, "</span>")
}

.compute_aa_offset_in_v_sequence <- function(AIRR_df)
{
    stopifnot(is.data.frame(AIRR_df))
    v_sequence_start <- AIRR_df$v_sequence_start
    cdr1_start <- AIRR_df$cdr1_start
    stopifnot(is.integer(v_sequence_start), is.integer(cdr1_start))
    (cdr1_start - v_sequence_start) %% 3L
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_header_line1()
### .make_header_line2()
###

.make_header_line1 <- function(AIRR_df)
{
    stopifnot(is.data.frame(AIRR_df))
    sequence_id <- AIRR_df$sequence_id
    stopifnot(is.character(sequence_id))
    attribs <- c("padding-left: 15px", "padding-right: 8px")
    span_tag <- sprintf("<span %s>", .make_style_attrib(attribs))
    seq_nb <- paste0(span_tag, "<b>", seq_along(sequence_id), ".</b></span>")
    paste0(seq_nb, "<i>sequence_id:</i>&nbsp;<b>", sequence_id, "</b>")
}

.HEADER_LINE2_SPAN_ATTRIBS <- c("padding: 2px",
                                "padding-left: 8px",
                                "padding-right: 8px",
                                "border-radius: 8px")

.make_VDJheadbox <- function(prefix, vdj_call, vdj_identity, bg_color)
{
    stopifnot(isSingleNonWhiteString(prefix),
              is.vector(vdj_call),
              is.numeric(vdj_identity),
              isSingleNonWhiteString(bg_color))
    vdj_identity <- format(vdj_identity, digits=4L)
    attribs <- c(.HEADER_LINE2_SPAN_ATTRIBS,
                 sprintf("background: %s", bg_color))
    span_tag <- sprintf("<span %s>", .make_style_attrib(attribs))
    paste0(span_tag,
           "<i>", prefix, "_call:</i>&nbsp;<b>", vdj_call, "</b>",
           "&nbsp;&nbsp;&nbsp;",
           "<i>", prefix, "_identity:</i>&nbsp;<b>", vdj_identity, "</b>",
           "</span>")
}

.make_header_line2 <- function(AIRR_df,
                               Vcolor="#FFDDD2", Dcolor="#CFC", Jcolor="#CEF",
                               Ccolor="#EEC")
{
    stopifnot(is.data.frame(AIRR_df),
              isSingleNonWhiteString(Vcolor),
              isSingleNonWhiteString(Dcolor),
              isSingleNonWhiteString(Jcolor),
              isSingleNonWhiteString(Ccolor))
    locus <- AIRR_df$locus
    v_call <- AIRR_df$v_call
    v_identity <- AIRR_df$v_identity
    d_call <- AIRR_df$d_call
    d_identity <- AIRR_df$d_identity
    j_call <- AIRR_df$j_call
    j_identity <- AIRR_df$j_identity
    stopifnot(is.character(locus),
              is.character(v_call), is.numeric(v_identity),
              is.vector(d_call), is.numeric(d_identity),
              is.character(j_call), is.numeric(j_identity))
    attribs <- c(.HEADER_LINE2_SPAN_ATTRIBS, "background: #DDD")
    span_tag <- sprintf("<span %s>", .make_style_attrib(attribs))
    locus <- paste0(span_tag, "<i>locus:</i>&nbsp;",
                    "<b>", locus, "</b></span>")
    ans <- paste0(locus,
                  "&nbsp;", .make_VDJheadbox("v", v_call, v_identity, Vcolor),
                  "&nbsp;", .make_VDJheadbox("d", d_call, d_identity, Dcolor),
                  "&nbsp;", .make_VDJheadbox("j", j_call, j_identity, Jcolor))
    if ("c_call" %in% colnames(AIRR_df)) {
        c_call <- AIRR_df$c_call
        c_identity <- AIRR_df$c_identity
        stopifnot(is.character(c_call), is.numeric(c_identity))
        ans <- paste0(ans, "&nbsp;",
                      .make_VDJheadbox("c", c_call, c_identity, Ccolor))
    }
    paste0("<p>", ans, "</p>")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_FWR_CDR_line()
###

.make_double_arrow_ascii <- function(width, text="", color=NULL)
{
    stopifnot(is.integer(width), isSingleString(text))
    text <- rep.int(text, length(width))
    text[width < nchar(text) + 2L] <- ""
    ndashes <- width - 2L - nchar(text)
    Ldashes <- ndashes %/% 2L
    Rdashes <- (ndashes + 1L) %/% 2L
    ans <- paste0("<", strrep("-", pmax(Ldashes, 0L)), text,
                       strrep("-", pmax(Rdashes, 0L)), ">")
    ans[width == 1L] <- "."
    ans[width == 0L] <- ""
    stopifnot(identical(nchar(ans), width))
    if (!is.null(color)) {
        span_tag <- sprintf("<span style=\"background: %s\">", color)
        ans <- paste0(span_tag, ans, "</span>")
    }
    ans
}

.make_fwr1_margin <- function(AIRR_df)
{
    stopifnot(is.data.frame(AIRR_df))
    v_sequence_start <- AIRR_df$v_sequence_start
    fwr1_start <- AIRR_df$fwr1_start
    stopifnot(is.integer(v_sequence_start), is.integer(fwr1_start),
              identical(v_sequence_start, fwr1_start))
    strrep(" ", fwr1_start - 1L)
}

.make_FWR_CDR_line <- function(AIRR_df, show.full.sequence,
                                        FWRcolor="#C9D",
                                        CDRcolor="#EE4")
{
    stopifnot(is.data.frame(AIRR_df))

    ## Extract columns of interest.
    fwr1_start <- AIRR_df$fwr1_start
    fwr1_end   <- AIRR_df$fwr1_end
    cdr1_start <- AIRR_df$cdr1_start
    cdr1_end   <- AIRR_df$cdr1_end
    fwr2_start <- AIRR_df$fwr2_start
    fwr2_end   <- AIRR_df$fwr2_end
    cdr2_start <- AIRR_df$cdr2_start
    cdr2_end   <- AIRR_df$cdr2_end
    fwr3_start <- AIRR_df$fwr3_start
    fwr3_end   <- AIRR_df$fwr3_end
    cdr3_start <- AIRR_df$cdr3_start
    cdr3_end   <- AIRR_df$cdr3_end
    fwr4_start <- AIRR_df$fwr4_start
    fwr4_end   <- AIRR_df$fwr4_end

    ## Sanity checks.
    stopifnot(is.integer(fwr1_start), is.integer(fwr1_end),
              is.integer(cdr1_start), is.integer(cdr1_end),
              is.integer(fwr2_start), is.integer(fwr2_end),
              is.integer(cdr2_start), is.integer(cdr2_end),
              is.integer(fwr3_start), is.integer(fwr3_end),
              is.integer(cdr3_start), is.integer(cdr3_end),
              is.integer(fwr4_start), is.integer(fwr4_end))
    stopifnot(!anyNA(fwr1_start), !anyNA(fwr1_end),
              !anyNA(cdr1_start), !anyNA(cdr1_end),
              !anyNA(fwr2_start), !anyNA(fwr2_end),
              !anyNA(cdr2_start), !anyNA(cdr2_end),
              !anyNA(fwr3_start), !anyNA(fwr3_end))
    stopifnot(all(fwr1_start < fwr1_end),
              all(fwr1_end + 1L == cdr1_start),
              all(cdr1_start < cdr1_end),
              all(cdr1_end + 1L == fwr2_start),
              all(fwr2_start < fwr2_end),
              all(fwr2_end + 1L == cdr2_start),
              all(cdr2_start < cdr2_end),
              all(cdr2_end + 1L == fwr3_start),
              all(fwr3_start < fwr3_end),
              all(fwr3_end + 1L == cdr3_start, na.rm=TRUE),
              all(cdr3_start < cdr3_end, na.rm=TRUE),
              all(cdr3_end + 1L == fwr4_start, na.rm=TRUE),
              all(fwr4_start < fwr4_end, na.rm=TRUE))

    fwr1_width <- fwr1_end - fwr1_start + 1L
    cdr1_width <- cdr1_end - cdr1_start + 1L
    fwr2_width <- fwr2_end - fwr2_start + 1L
    cdr2_width <- cdr2_end - cdr2_start + 1L
    fwr3_width <- fwr3_end - fwr3_start + 1L
    cdr3_width <- cdr3_end - cdr3_start + 1L
    fwr4_width <- fwr4_end - fwr4_start + 1L
    cdr3fwr4_missing <- is.na(cdr3_start)
    cdr3_width[cdr3fwr4_missing] <- fwr4_width[cdr3fwr4_missing] <- 0L

    html_fwr1 <- .make_double_arrow_ascii(fwr1_width, "FWR1", color=FWRcolor)
    html_cdr1 <- .make_double_arrow_ascii(cdr1_width, "CDR1", color=CDRcolor)
    html_fwr2 <- .make_double_arrow_ascii(fwr2_width, "FWR2", color=FWRcolor)
    html_cdr2 <- .make_double_arrow_ascii(cdr2_width, "CDR2", color=CDRcolor)
    html_fwr3 <- .make_double_arrow_ascii(fwr3_width, "FWR3", color=FWRcolor)
    html_cdr3 <- .make_double_arrow_ascii(cdr3_width, "CDR3", color=CDRcolor)
    html_fwr4 <- .make_double_arrow_ascii(fwr4_width, "FWR4", color=FWRcolor)
    html_cdr3[cdr3fwr4_missing] <-
        "<span style=\"color: red\"> [CDR3/FWR4 info not available!]</span>"

    ans <- paste0(html_fwr1, html_cdr1, html_fwr2,
                  html_cdr2, html_fwr3, html_cdr3, html_fwr4)
    if (show.full.sequence) {
        fwr1_margin <- .make_fwr1_margin(AIRR_df)
        ans <- paste0(fwr1_margin, ans)
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_sequence_aa_line()
###

.make_aa_margin_in_v_sequence <- function(AIRR_df, Vcolor="#FFDDD2")
{
    aa_offset_in_v_sequence <- .compute_aa_offset_in_v_sequence(AIRR_df)
    paste0(sprintf("<span style=\"background: %s\">", Vcolor),
           strrep(" ", aa_offset_in_v_sequence),
           "</span>")
}

.stretch_3x_sequence_aa <- function(aa)
{
    stopifnot(isSingleString(aa))
    exploded_aa <- safeExplode(aa)
    codon_colors <- c("#DDD", "#EEE")
    span_tags <- sprintf("<span style=\"background: %s\">", codon_colors)
    span_tags <- rep_len(span_tags, length.out=length(exploded_aa))
    paste0(paste0(span_tags, " ", exploded_aa, " </span>"), collapse="")
}

.compute_v_offset <- function(AIRR_df)
{
    stopifnot(is.data.frame(AIRR_df))
    v_sequence_start <- AIRR_df$v_sequence_start
    stopifnot(is.integer(v_sequence_start))
    v_sequence_start - 1L
}

.make_sequence_aa_line <- function(AIRR_df, show.full.sequence,
                                            Vcolor="#FFDDD2")
{
    stopifnot(is.data.frame(AIRR_df))
    sequence_aa <- AIRR_df$sequence_aa
    stopifnot(is.character(sequence_aa))
    margin <- .make_aa_margin_in_v_sequence(AIRR_df, Vcolor=Vcolor)
    stretched_aa <- vapply(sequence_aa, .stretch_3x_sequence_aa, character(1))
    ans <- paste0(margin, stretched_aa)
    if (show.full.sequence) {
        v_offset <- .compute_v_offset(AIRR_df)
        v_margin <- strrep(" ", v_offset)
        ans <- paste0(v_margin, ans)
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_dna_line()
###

### Extract and check v_, d_, j_, c_ columns of interest.
.extract_vdjc_COI <- function(AIRR_df)
{
    stopifnot(is.data.frame(AIRR_df))

    sequence <- AIRR_df$sequence
    v_sequence_start <- AIRR_df$v_sequence_start
    v_sequence_end   <- AIRR_df$v_sequence_end
    d_sequence_start <- AIRR_df$d_sequence_start
    d_sequence_end   <- AIRR_df$d_sequence_end
    j_sequence_start <- AIRR_df$j_sequence_start
    j_sequence_end   <- AIRR_df$j_sequence_end

    ## Sanity checks.
    stopifnot(is.character(sequence),
              is.integer(v_sequence_start), is.integer(v_sequence_end),
              is.integer(d_sequence_start), is.integer(d_sequence_end),
              is.integer(j_sequence_start), is.integer(j_sequence_end))
    stopifnot(!anyNA(v_sequence_start), !anyNA(v_sequence_end),
              identical(is.na(d_sequence_start), is.na(d_sequence_end)),
              !anyNA(j_sequence_start), !anyNA(j_sequence_end))
    stopifnot(all(v_sequence_start <= v_sequence_end),
              all(v_sequence_end < d_sequence_start, na.rm=TRUE),
              all(d_sequence_start <= d_sequence_end, na.rm=TRUE),
              all(d_sequence_end < j_sequence_start, na.rm=TRUE),
              all(j_sequence_start <= j_sequence_end),
              all(v_sequence_end < j_sequence_start))

    COI <- list(
        sequence=sequence,
        v_sequence_start=v_sequence_start,
        v_sequence_end=v_sequence_end,
        d_sequence_start=d_sequence_start,
        d_sequence_end=d_sequence_end,
        j_sequence_start=j_sequence_start,
        j_sequence_end=j_sequence_end,
        has_no_D=is.na(d_sequence_start)
    )

    if ("c_call" %in% colnames(AIRR_df)) {
        ## Strangely, for most sequences, 'j_sequence_end' and
        ## 'c_sequence_start' are the same! Are 'c_sequence_start' values
        ## zero-based? We correct this by adding 1 to 'c_sequence_start'.
        c_sequence_start <- AIRR_df$c_sequence_start
        c_sequence_end   <- AIRR_df$c_sequence_end
        stopifnot(is.integer(c_sequence_start), is.integer(c_sequence_end),
                  identical(is.na(c_sequence_start), is.na(c_sequence_end)),
                  all(j_sequence_end <= c_sequence_start, na.rm=TRUE),
                  all(c_sequence_start <= c_sequence_end, na.rm=TRUE))
        COI <- c(COI, list(
            c_sequence_start=c_sequence_start+1L,  # increase by one
            c_sequence_end=c_sequence_end,
            has_no_C=is.na(c_sequence_start)
        ))
    }

    COI
}

.toupper_with_lower_prefix <- function(x, prefix_width)
{
    stopifnot(is.character(x), is.integer(prefix_width),
              length(x) == length(prefix_width))
    lower <- tolower(substr(x, 1L, prefix_width))
    upper <- toupper(substr(x, prefix_width + 1L, nchar(x)))
    ans <- paste0(lower, upper)
    ans[is.na(x)] <- NA_character_
    stopifnot(identical(nchar(x), nchar(ans)))
    ans
}

.color_DNA_sequence <- function(dna, nuc_colors)
{
    stopifnot(isSingleString(dna))
    exploded_dna <- safeExplode(dna)
    colors <- nuc_colors[toupper(exploded_dna)]
    span_tags <- sprintf("<span style=\"color: %s\">", colors)
    colored <- paste0(
        ifelse(is.na(colors), "", span_tags),
        exploded_dna,
        ifelse(is.na(colors), "", "</span>"))
    paste0(colored, collapse="")
}

.color_DNA_sequences <- function(sequences, nuc_colors)
{
    stopifnot(is.character(sequences),
              is.character(nuc_colors),
              identical(names(nuc_colors), c("A", "C", "G", "T")))
    vapply(sequences, .color_DNA_sequence, character(1), nuc_colors=nuc_colors)
}

.make_dna_line <- function(AIRR_df, show.full.sequence, nuc_colors,
                           Vcolor="#FFDDD2", Dcolor="#CFC", Jcolor="#CEF",
                           Ccolor="#EEC")
{
    stopifnot(isSingleNonWhiteString(Vcolor),
              isSingleNonWhiteString(Dcolor),
              isSingleNonWhiteString(Jcolor),
              isSingleNonWhiteString(Ccolor))

    COI <- .extract_vdjc_COI(AIRR_df)
    has_no_D <- COI$has_no_D
    has_no_C <- COI$has_no_C  # possibly NULL!

    L_segment  <- substr(COI$sequence, 1L, COI$v_sequence_start - 1L)
    v_segment  <- substr(COI$sequence, COI$v_sequence_start,
                                       COI$v_sequence_end)
    vd_segment <- substr(COI$sequence, COI$v_sequence_end + 1L,
                                       COI$d_sequence_start - 1L)
    d_segment  <- substr(COI$sequence, COI$d_sequence_start,
                                       COI$d_sequence_end)
    dj_segment <- substr(COI$sequence, COI$d_sequence_end + 1L,
                                       COI$j_sequence_start - 1L)
    j_segment  <- substr(COI$sequence, COI$j_sequence_start,
                                       COI$j_sequence_end)
    R_segment  <- substr(COI$sequence, COI$j_sequence_end + 1L,
                                       nchar(COI$sequence))
    vd_segment[has_no_D] <- d_segment[has_no_D] <- dj_segment[has_no_D] <- ""
    if (!is.null(has_no_C)) {
        jc_segment <- substr(COI$sequence, COI$j_sequence_end + 1L,
                                           COI$c_sequence_start - 1L)
        c_segment  <- substr(COI$sequence, COI$c_sequence_start,
                                           COI$c_sequence_end)
        jc_segment[has_no_C] <- c_segment[has_no_C] <- ""
        R_segment[!has_no_C] <-
            substr(COI$sequence, COI$c_sequence_end + 1L,
                                 nchar(COI$sequence))[!has_no_C]
    }

    mid_segment1 <- paste0(vd_segment, d_segment, dj_segment)
    mid_segment2 <- substr(COI$sequence, COI$v_sequence_end + 1L,
                                         COI$j_sequence_start - 1L)
    mid_segment <- ifelse(has_no_D, mid_segment2, mid_segment1)
    sequence2 <- paste0(L_segment, v_segment, mid_segment, j_segment)
    if (!is.null(has_no_C)) {
        sequence2[!has_no_C] <-
            paste0(sequence2, jc_segment, c_segment)[!has_no_C]
    }
    sequence2 <- paste0(sequence2, R_segment)
    stopifnot(identical(COI$sequence, sequence2))

    ## Set case.
    aa_offset_in_v_sequence <- .compute_aa_offset_in_v_sequence(AIRR_df)
    v_segment <- .toupper_with_lower_prefix(v_segment, aa_offset_in_v_sequence)
    vd_segment <- toupper(vd_segment)
    d_segment <- toupper(d_segment)
    dj_segment <- toupper(dj_segment)
    j_segment <- toupper(j_segment)
    R_segment <- toupper(R_segment)
    mid_segment2 <- toupper(mid_segment2)
    if (!is.null(has_no_C)) {
        jc_segment <- toupper(jc_segment)
        c_segment <- toupper(c_segment)
    }

    if (!is.null(nuc_colors)) {
        ## Color DNA letters.
        v_segment  <- .color_DNA_sequences(v_segment, nuc_colors)
        vd_segment <- .color_DNA_sequences(vd_segment, nuc_colors)
        d_segment  <- .color_DNA_sequences(d_segment, nuc_colors)
        dj_segment <- .color_DNA_sequences(dj_segment, nuc_colors)
        j_segment  <- .color_DNA_sequences(j_segment, nuc_colors)
        R_segment  <- .color_DNA_sequences(R_segment, nuc_colors)
        mid_segment2  <- .color_DNA_sequences(mid_segment2, nuc_colors)
        if (!is.null(has_no_C)) {
            jc_segment <- .color_DNA_sequences(jc_segment, nuc_colors)
            c_segment  <- .color_DNA_sequences(c_segment, nuc_colors)
        }
    }

    ## Set germline segment background color.
    v_segment <- .set_background_color(v_segment, Vcolor)
    d_segment <- .set_background_color(d_segment, Dcolor)
    j_segment <- .set_background_color(j_segment, Jcolor)
    if (!is.null(has_no_C))
        c_segment <- .set_background_color(c_segment, Ccolor)

    mid_segment1 <- paste0(vd_segment, d_segment, dj_segment)
    mid_segment <- ifelse(has_no_D, mid_segment2, mid_segment1)
    ans <- paste0(v_segment, mid_segment, j_segment)
    if (show.full.sequence) {
        L_segment <- tolower(L_segment)
        L_segment  <- .color_DNA_sequences(L_segment, nuc_colors)
        ans <- paste0(L_segment, ans)
    }
    if (!is.null(has_no_C))
        ans <- paste0(ans, jc_segment, c_segment)
    paste0(ans, R_segment)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igbrowser()
###

.display_html_in_browser <- function(html)
{
    stopifnot(is.character(html))
    html_file <- tempfile(fileext = ".html")
    writeLines(html, html_file)
    temp_url <- paste0("file://", html_file)
    browseURL(temp_url)
}

igbrowser <- function(AIRR_df,
                      show.full.sequence=FALSE,
                      dna.coloring=TRUE,
                      Vcolor="#FFDDD2",
                      Dcolor="#CFC",
                      Jcolor="#CEF",
                      Ccolor="#EEC",
                      FWRcolor="#C9D",
                      CDRcolor="#EE4")
{
    stopifnot(is.data.frame(AIRR_df),
              isTRUEorFALSE(show.full.sequence),
              isTRUEorFALSE(dna.coloring))

    header_line1 <- .make_header_line1(AIRR_df)
    header_line2 <- .make_header_line2(AIRR_df, Vcolor=Vcolor,
                                                Dcolor=Dcolor,
                                                Jcolor=Jcolor,
                                                Ccolor=Ccolor)
    sequence_aa_line <- .make_sequence_aa_line(AIRR_df, show.full.sequence,
                                                        Vcolor=Vcolor)
    FWR_CDR_line <- .make_FWR_CDR_line(AIRR_df, show.full.sequence,
                                                FWRcolor=FWRcolor,
                                                CDRcolor=CDRcolor)
    nuc_colors <- if (dna.coloring) .NUCLEOTIDE_COLORS else NULL
    dna_line1 <- .make_dna_line(AIRR_df, show.full.sequence,
                                         nuc_colors=nuc_colors,
                                         Vcolor=Vcolor,
                                         Dcolor=Dcolor,
                                         Jcolor=Jcolor,
                                         Ccolor=Ccolor)
    nuc_colors <- setNames(rep.int("transparent", 4L), c("A", "C", "G", "T"))
    dna_line2 <- .make_dna_line(AIRR_df, show.full.sequence,
                                         nuc_colors=nuc_colors,
                                         Vcolor=Vcolor,
                                         Dcolor=Dcolor,
                                         Jcolor=Jcolor,
                                         Ccolor=Ccolor)
    attribs <- c("border-top: 1px solid #BBB",
                 "padding: 5px", "padding-left: 0px")
    td_tag <- sprintf("<td %s>", .make_style_attrib(attribs))
    html <- paste0(
        "<tr>", td_tag, "\n",
        header_line1, "\n",
        header_line2, "\n",
        "<pre style=\"margin: 0px; margin-bottom: 6px;\">",
        FWR_CDR_line, "\n",
        sequence_aa_line, "\n",
        dna_line1, "\n",
        dna_line2,
        "</pre>\n</td></tr>\n")
    html <- c("<html><head><title>igbrowser()</title></head>",
              "<body><table style=\"border-collapse: collapse\">",
              html,
              "</table></body></html>\n")
    .display_html_in_browser(html)
}

bcr_browser <- function(...)
{
    .Deprecated("igbrowser")
    igbrowser(...)
}

