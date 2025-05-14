### =========================================================================
### Various general purpose low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### TODO: wmsg() was replaced with this in S4Vectors >= 0.45.1 so drop
### wmsg2() and use wmsg() instead.
wmsg2 <- function(..., margin=2)
{
    width <- getOption("width") - margin
    paste0(strwrap(paste0(c(...), collapse=""), width=width),
           collapse=paste0("\n", strrep(" ", margin)))
}

### TODO: Move this to S4Vectors (or BiocBaseUtils).
load_package_gracefully <- function(package, ...)
{
    if (!requireNamespace(package, quietly=TRUE))
        stop("Could not load package ", package, ". Is it installed?\n\n  ",
             wmsg("Note that ", ..., " requires the ", package, " package. ",
                  "Please install it with:"),
             "\n\n    BiocManager::install(\"", package, "\")")
}

### "\xc2\xa0" is some kind of weird white space that sometimes creeps
### in when scrapping dirty HTML documents found on the internet.
.WHITESPACES <- c(" ", "\t", "\r", "\n", "\xc2\xa0")

### Vectorized. Note that NAs do **not** get propagated. NA elements in 'x'
### produce FALSE elements in the output.
has_whitespace <- function(x)
{
    stopifnot(is.character(x))
    pattern <- paste0("[", paste(.WHITESPACES, collapse=""), "]")
    grepl(pattern, x, perl=TRUE)
}

### A simple wrapper to base::trimws() that starts by replacing **all**
### whitespaces in 'x' with regular spaces (" "), even non-leading
### and non-trailing whitespaces.
### Like 'base::trimws()', trimws2() is vectorized and propagates NAs.
trimws2 <- function(x)
{
    stopifnot(is.character(x))
    old <- paste(.WHITESPACES, collapse="")
    new <- strrep(" ", nchar(old))
    trimws(chartr(old, new, x), whitespace=" ")
}

### Vectorized. Note that NAs do **not** get propagated. NA elements in 'x'
### produce FALSE elements in the output.
is_white_str <- function(x) !nzchar(trimws2(x))

isSingleNonWhiteString <- function(x) isSingleString(x) && !is_white_str(x)

drop_heading_and_trailing_white_lines <- function(lines)
{
    stopifnot(is.character(lines))
    ok <- vapply(lines, is_white_str, logical(1), USE.NAMES=FALSE)
    nonwhite_idx <- which(!ok)
    if (length(nonwhite_idx) == 0L) {
        keep_idx <- integer(0)
    } else {
        keep_idx <- (nonwhite_idx[[1L]]):(nonwhite_idx[[length(nonwhite_idx)]])
    }
    lines[keep_idx]
}

### Vectorized.
has_prefix <- function(x, prefix)
{
    stopifnot(is.character(x), isSingleString(prefix))
    substr(x, 1L, nchar(prefix)) == prefix
}

### Vectorized.
has_suffix <- function(x, suffix)
{
    stopifnot(is.character(x), isSingleString(suffix))
    x_nc <- nchar(x)
    substr(x, x_nc - nchar(suffix) + 1L, x_nc) == suffix
}

### Not used at the moment.
strslice <- function(x, width)
{
    stopifnot(isSingleString(x), isSingleNumber(width))
    chunks <- breakInChunks(nchar(x), chunksize=width)
    vapply(seq_along(chunks),
        function(i) substr(x, start(chunks)[i], end(chunks)[i]),
        character(1), USE.NAMES=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### interweave()
###

### 'x' and 'y' must be vectors or vector-like objects of the same length.
### They must support c() and subsetting.
interweave <- function(x, y)
{
    N <- length(x)
    stopifnot(length(y) == N)
    c(x, y)[S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(N)]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### align_vectors_by_names()
###

align_vectors_by_names <- function(vectors)
{
    stopifnot(is.list(vectors), length(vectors) != 0L)
    all_names <- lapply(vectors,
        function(v) {
            nms <- names(v)
            if (is.null(nms))
                stop(wmsg("all vectors must be named"))
            if (anyDuplicated(nms))
                stop(wmsg("some vectors have duplicated names"))
            nms
        })
    unique_names <- unique(unlist(all_names))
    ans <- lapply(vectors,
        function(v) setNames(v[unique_names], unique_names))
    stopifnot(all(lengths(ans) == length(ans[[1L]])))
    ans
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
    df <- as.data.frame(mat)
    df[[2L]] <- as.Date(df[[2L]])
    df
}

### 'css' must be a single string specifying the CSS selector to the table
### containing the index e.g. "body" or "body section".
### Returns a data.frame with 3 columns: Name, Last modified, Size
scrape_html_dir_index <- function(url, css="body", suffix=NULL)
{
    stopifnot(isSingleNonWhiteString(url), isSingleNonWhiteString(css))
    html <- getUrlContent(url, type="text", encoding="UTF-8")
    xml <- read_html(html)
    all_ths <- html_text(html_elements(xml, paste0(css, " table tr th")))
    all_tds <- html_text(html_elements(xml, paste0(css, " table tr td")))
    EXPECTED_NCOL <- 5L
    mat <- matrix(all_tds, ncol=EXPECTED_NCOL, byrow=TRUE)
    colnames(mat) <- all_ths[seq_len(EXPECTED_NCOL)]
    .make_df_from_matrix_of_tds(mat, suffix=suffix)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous stuff
###

websiteIsUp <- function(url)
{
    if (!has_internet())
        stop("no internet")
    response <- try(HEAD(url, user_agent("igblastr")), silent=TRUE)
    !inherits(response, "try-error")
}

urlExists <- function(url)
{
    response <- try(HEAD(url, user_agent("igblastr")), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(as.character(response), "  Please check your internet connection.")
    response$status_code != 404L
}

getUrlContent <- function(url, query=list(), type=NULL, encoding=NULL)
{
    stopifnot(is.list(query))
    if (length(query) != 0L)
        stopifnot(!is.null(names(query)))
    response <- try(GET(url, user_agent("igblastr"), query=query), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(as.character(response), "  Please check your internet connection.")
    if (response$status_code == 404L)
        stop(wmsg("Not Found (HTTP 404): ", url))
    stop_for_status(response)
    content(response, type=type, encoding=encoding)
}

### Note that in case of user interrupt (CTRL+C) we can end up with a
### partial download and corrupted file! Can we achieve atomic behavior?
### TODO: Try to make the behavior atomic, under any circumstance.
download_as_tempfile <- function(dir_url, filename, ...)
{
    url <- paste0(dir_url, filename)
    destfile <- tempfile()
    code <- try(suppressWarnings(download.file(url, destfile, ...)),
                silent=TRUE)
    if (inherits(code, "try-error") || code != 0L)
        stop(wmsg("Failed to download ", filename, " ",
                  "from ", dir_url, ". ",
                  "Are you connected to the internet?"))
    destfile
}

### A thin wrapper to untar() with more user-friendly error handling.
### 'exdir' should be the path to an existing directory that is
### preferrably empty.
untar2 <- function(tarfile, original_tarball_name, exdir=".")
{
    stopifnot(isSingleNonWhiteString(tarfile),
              isSingleNonWhiteString(original_tarball_name),
              isSingleNonWhiteString(exdir),
              dir.exists(exdir))
    code <- suppressWarnings(untar(tarfile, exdir=exdir))
    if (code != 0L)
        stop(wmsg("Anomaly: something went wrong during ",
                  "extraction of '", tarfile, "' (the local copy ",
                  "of '", original_tarball_name, "') to '", exdir, "'."))
}

### Returns the OS (e.g. Linux, Windows, or Darwin) and arch (e.g. x86_64
### or arm64) in a character vector of length 2, with names "OS" and "arch".
### Note that if the OS or arch cannot be obtained with Sys.info() then they
### get replaced with an NA.
get_OS_arch <- function()
{
    sys_info <- Sys.info()
    sysname <- sys_info[["sysname"]]
    if (!isSingleNonWhiteString(sysname))
        sysname <- NA_character_
    machine <- sys_info[["machine"]]
    if (!isSingleNonWhiteString(machine))
        machine <- NA_character_
    c(OS=sysname, arch=machine)
}

add_exe_suffix_on_Windows <- function(files, OS=get_OS_arch()[["OS"]])
{
    stopifnot(is.character(files), isSingleStringOrNA(OS))
    if (length(files) == 0L || is.na(OS) || !grepl("^win", tolower(OS)))
        return(files)
    paste0(files, ".exe")
}

### To use on the result of 'try(system2(..., stdout=TRUE, stderr=TRUE))'.
system_command_worked <- function(out)
{
    if (inherits(out, "try-error"))
        return(FALSE)
    status <- attr(out, "status")
    is.null(status) || isTRUE(all.equal(status, 0L))
}

system_command_works <- function(command, args=character())
{
    out <- try(suppressWarnings(system2(command, args=args,
                                        stdout=TRUE, stderr=TRUE)),
               silent=TRUE)
    system_command_worked(out)
}

has_perl <- function() system_command_works("perl", args="-v")

system3 <- function(command, outfile, errfile, args=character())
{
    status <- system2(command, args=args, stdout=outfile, stderr=errfile)
    if (file.exists(errfile)) {
        errmsg <- readLines(errfile)
        if (length(errmsg) != 0L)
            stop(paste(errmsg, collapse="\n"))
        unlink(errfile)
    }
    if (status != 0) {
        cmd_in_1string <- paste(c(command, args), collapse=" ")
        stop(wmsg("command '", cmd_in_1string, "' failed"))
    }
}

display_local_file_in_browser <- function(file)
{
    top_html <- tempfile()
    writeLines("<PRE>", top_html)
    bottom_html <- tempfile()
    writeLines("</PRE>", bottom_html)
    temp_html <- tempfile(fileext=".html")
    concatenate_files(c(top_html, file, bottom_html), out=temp_html)
    temp_url <- paste0("file://", temp_html)
    browseURL(temp_url)
}

display_data_frame_in_browser <- function(df)
{
    temp_html <- tempfile(fileext=".html")
    print(xtable(df), type="html", file=temp_html)
    temp_url <- paste0("file://", temp_html)
    browseURL(temp_url)
}

