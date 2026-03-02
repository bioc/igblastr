### =========================================================================
### Various low-level internet-related utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### Additional curl configuration options can be passed thru the ellipsis
### as **named** arguments e.g. websiteIsUp(url, connecttimeout=20).
### See ?httr::httr_options for all available options.
websiteIsUp <- function(url, ...)
{
    if (!has_internet())
        stop(wmsg("no internet"))
    config <- config(...)
    response <- try(HEAD(url, config, user_agent("igblastr")), silent=TRUE)
    !inherits(response, "try-error")
}

### Additional curl configuration options can be passed thru the ellipsis
### as **named** arguments e.g. urlExists(url, connecttimeout=20).
### See ?httr::httr_options for all available options.
urlExists <- function(url, ...)
{
    if (!has_internet())
        stop(wmsg("no internet"))
    config <- config(...)
    response <- try(HEAD(url, config, user_agent("igblastr")), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(wmsg(as.character(response)))
    response$status_code != 404L
}

### Additional curl configuration options can be passed thru the ellipsis
### as **named** arguments. See ?httr::httr_options for all available options.
getUrlContent <- function(url, query=list(), type=NULL, encoding=NULL, ...)
{
    stopifnot(is.list(query))
    if (length(query) != 0L)
        stopifnot(!is.null(names(query)))
    if (!has_internet())
        stop(wmsg("no internet"))
    config <- config(...)
    response <- try(GET(url, config, user_agent("igblastr"), query=query),
                    silent=TRUE)
    if (inherits(response, "try-error"))
        stop(wmsg(as.character(response)))
    if (response$status_code == 404L)
        stop(wmsg("Not Found (HTTP 404): ", url))
    stop_for_status(response)
    content(response, type=type, encoding=encoding)
}


### A thin wrapper around utils::download.file() that tries to be more
### cautious.
### Returns 'destfile'.
### NOTE: A major problem with download.file() is that, in case of user
### interrupt (CTRL+C) or connectivity issue, we can end up with a partial
### download and corrupted file! Can we achieve atomic behavior?
### TODO: Try to make the behavior atomic, under any circumstance.
download_file <- function(url, destfile, ...)
{
    stopifnot(isSingleNonWhiteString(url), isSingleNonWhiteString(destfile))
    if (!has_internet())
        stop(wmsg("no internet"))
    code <- try(suppressWarnings(download.file(url, destfile, ...)),
                silent=TRUE)
    if (inherits(code, "try-error") || code != 0L)
        stop(wmsg("failed to download ", url))
    destfile
}

### NOT atomic! (see above NOTE and TODO)
download_as_tempfile <- function(dir_url, filename, ...)
{
    stopifnot(isSingleNonWhiteString(dir_url), isSingleNonWhiteString(filename))
    if (!has_internet())
        stop(wmsg("no internet"))
    url <- paste0(dir_url, filename)
    destfile <- tempfile()
    code <- try(suppressWarnings(download.file(url, destfile, ...)),
                silent=TRUE)
    if (inherits(code, "try-error") || code != 0L)
        stop(wmsg("failed to download ", filename, " from ", dir_url))
    destfile
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### scrape_html_dir_index()
###

.make_df_from_matrix_of_tds <- function(mat, suffix=NULL)
{
    stopifnot(is.matrix(mat), !is.null(colnames(mat)))

    EXPECTED_COLNAMES <- c("Name", "Last modified", "Size")
    m <- match(tolower(EXPECTED_COLNAMES), tolower(colnames(mat)))
    stopifnot(!anyNA(mat))
    mat <- mat[ , m, drop=FALSE]
    colnames(mat) <- EXPECTED_COLNAMES

    ## Remove "Parent Directory" row.
    if (tolower(trimws(mat[1L, 1L])) == "parent directory")
        mat <- mat[-1L, , drop=FALSE]

    if (!is.null(suffix))
        mat <- mat[has_suffix(mat[ , 1L], suffix), , drop=FALSE]
    as.data.frame(mat)
}

### Use for IMGT-style HTML directory index.
.IMGT_html_dir_parser <- function(xml, suffix=NULL)
{
    stopifnot(inherits(xml, "xml_document"))
    pre <- html_text(html_elements(xml, "body pre"))
    stopifnot(length(pre) == 1L)
    pre <- sub("\\.\\./", "", pre)
    pre <- strsplit(pre, "\r?\n")[[1L]]
    pre <- pre[nchar(pre) != 0L]
    pre <- strsplit(pre, " +")
    EXPECTED_NCOL <- 4L
    stopifnot(all(lengths(pre) == EXPECTED_NCOL))
    mat <- matrix(unlist(pre), ncol=EXPECTED_NCOL, byrow=TRUE)
    colnames(mat) <- c("Name", "Last modified", "Time last modified", "Size")
    df <- .make_df_from_matrix_of_tds(mat, suffix=suffix)
    df[[2L]] <- as.Date(df[[2L]], "%d-%b-%Y")
    df
}

### Use for OAS-style HTML directory index.
.OAS_html_dir_parser <- function(xml, suffix=NULL)
{
    stopifnot(inherits(xml, "xml_document"))
    EXPECTED_NCOL <- 5L
    all_ths <- html_text(html_elements(xml, "body table tr th"))
    stopifnot(length(all_ths) != 0L)
    all_tds <- html_text(html_elements(xml, "body table tr td"))
    stopifnot(length(all_tds) %% EXPECTED_NCOL == 0L)
    mat <- matrix(all_tds, ncol=EXPECTED_NCOL, byrow=TRUE)
    colnames(mat) <- all_ths[seq_len(EXPECTED_NCOL)]
    df <- .make_df_from_matrix_of_tds(mat, suffix=suffix)
    df[[2L]] <- as.Date(df[[2L]])
    df
}

### Additional curl configuration options can be passed thru the ellipsis
### as **named** arguments. See ?httr::httr_options for all available options.
### Returns a data.frame with 3 columns: Name, Last modified, Size
scrape_html_dir_index <- function(url, style=c("IMGT", "OAS"), suffix=NULL,
                                  ...)
{
    stopifnot(isSingleNonWhiteString(url), isSingleNonWhiteString(style))
    style <- match.arg(style)
    html <- getUrlContent(url, type="text", encoding="UTF-8", ...)
    xml <- read_html(html)
    fun <- if (style == "IMGT") .IMGT_html_dir_parser else .OAS_html_dir_parser
    fun(xml, suffix=suffix)
}

