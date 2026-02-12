### =========================================================================
### Various general purpose low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### TODO: wmsg() was replaced with this in S4Vectors >= 0.45.1 so drop
### wmsg2() and use wmsg() instead in BioC >= 3.23.
wmsg2 <- function(..., margin=2)
{
    width <- getOption("width") - margin
    paste0(strwrap(paste0(c(...), collapse=""), width=width),
           collapse=paste0("\n", strrep(" ", margin)))
}

### TODO: load_package_gracefully() was added to S4Vectors 0.47.6 so drop
### this and use S4Vectors:::load_package_gracefully() instead in BioC >= 3.23.
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

### TODO: has_prefix() was added to S4Vectors 0.47.5 so drop this and
### use S4Vectors:::has_prefix() instead in BioC >= 3.23.
has_prefix <- function(x, prefix)
{
    stopifnot(is.character(x), isSingleString(prefix))
    substr(x, 1L, nchar(prefix)) == prefix
}

### TODO: has_suffix() was added to S4Vectors 0.47.5 so drop this and
### use S4Vectors:::has_suffix() instead in BioC >= 3.23.
has_suffix <- function(x, suffix)
{
    stopifnot(is.character(x), isSingleString(suffix))
    x_nc <- nchar(x)
    substr(x, x_nc - nchar(suffix) + 1L, x_nc) == suffix
}

### Not used at the moment.
.tabulate_strings_by_prefix <- function(x, prefixes)
{
    stopifnot(is.character(x), is.character(prefixes), length(prefixes) >= 1L)
    nc <- nchar(prefixes)
    stopifnot(all(nc == nc[[1L]]))
    x_prefixes <- substr(x, 1L, nc)
    m <- match(x_prefixes, prefixes)
    if (anyNA(m)) {
        in1string <- paste0(prefixes, collapse=", ")
        stop(wmsg("all strings in 'x' must start with one of ",
                  "the following prefixes: ", in1string))
    }
    setNames(tabulate(m, length(prefixes)), prefixes)
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

check_seqlens <- function(seqlens, varname)
{
    stopifnot(is.integer(seqlens))
    seqids <- names(seqlens)
    stopifnot(!is.null(seqids))
    empty_idx <- which(seqlens == 0L)
    if (length(empty_idx) != 0L) {
        empty_seqids <- trimws2(seqids[empty_idx])
        if (!all(nzchar(empty_seqids)))
            stop(wmsg("some sequences in '", varname, "' are empty"))
        in1string <- paste(empty_seqids, collapse=", ")
        stop(wmsg("the following sequences in '", varname, "' ",
                  "are empty (showing seq ids): ", in1string))
    }
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
### read_version_file()
###

read_version_file <- function(dirpath)
{
    stopifnot(isSingleNonWhiteString(dirpath), dir.exists(dirpath))
    version_path <- file.path(dirpath, "version")
    if (!file.exists(version_path))
        stop(wmsg("missing 'version' file in ", dirpath, "/"))
    version <- readLines(version_path)
    if (length(version) != 1L)
        stop(wmsg("file '", version_path, "' should contain exactly one line"))
    version <- trimws2(version)
    if (version == "")
        stop(wmsg("file '", version_path, "' contains only white spaces"))
    version
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### drop_repeated_sequences()
###

### "Repeated" sequences in 'x' are sequences with identical names **and**
### identical content. Note that:
### - we keep sequences with identical content but different names;
### - we also keep sequences with identical names but different content.
drop_repeated_sequences <- function(x, what="", from="")
{
    stopifnot(is.character(x) || is(x, "XStringSet"),
              isSingleString(from))
    x_names <- names(x)
    stopifnot(!is.null(x_names))
    split_by_name <- splitAsList(x, x_names)
    drop_idx <- which(unsplit(duplicated(split_by_name), x_names))
    drop_idx_len <- length(drop_idx)
    if (drop_idx_len == 0L)
        return(x)  # no-op
    if (what != "" && from != "") {
        in1string <- paste(x_names[drop_idx], collapse=", ")
        msg <- c("Deleted ", drop_idx_len, " repeated ", what, " ",
                 "(out of ", length(x), ") from ", from, ": ", in1string)
        message("  o ", wmsg(msg, margin=4L), "\n")
    }
    x[-drop_idx]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### count_bin_hits()
###

### The "bins" are the adjacent integer intervals implicitly defined by
### integer vector 'bin_widths'. The values in 'bin_widths' must be
### non-negative integers that specify the "length" of each interval.
### Note that what we call the "length" (or "width") of an integer interval
### is simply the number of (integer) values in it. All the intervals are
### considered **adjacent** intervals, with the first interval starting on 1.
### So for example, setting 'bin_widths' to 'c(10L, 6L)' is defining the
### following 2 bins: the 1st bin is interval [1, 10], and the 2nd bin is
### interval [11, 16].
### Returns an integer vector parallel to 'bin_widths' where each element
### is the number of values in 'x' that fall in the corresponding bin.
### NAs in 'x' are ignored, as are values in 'x' that don't fall in any bin.
count_bin_hits <- function(x, bin_widths)
{
    stopifnot(is.integer(x), is.integer(bin_widths),
              !anyNA(bin_widths), all(bin_widths >= 0L))
    bin <- findInterval(x, c(0L, cumsum(bin_widths)), left.open=TRUE)
    setNames(tabulate(bin, nbins=length(bin_widths)), names(bin_widths))
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

### 'css' must be a single string specifying the CSS selector to
### the table containing the index e.g. "body" or "body section".
### Additional curl configuration options can be passed thru the ellipsis
### as **named** arguments. See ?httr::httr_options for all available options.
### Returns a data.frame with 3 columns: Name, Last modified, Size
scrape_html_dir_index <- function(url, css="body", suffix=NULL, ...)
{
    stopifnot(isSingleNonWhiteString(url), isSingleNonWhiteString(css))
    html <- getUrlContent(url, type="text", encoding="UTF-8", ...)
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

### Additional curl configuration options can be passed thru the ellipsis
### as **named** arguments e.g. websiteIsUp(url, connecttimeout=20).
### See ?httr::httr_options for all available options.
websiteIsUp <- function(url, ...)
{
    if (!has_internet())
        stop("no internet")
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
        stop("no internet")
    config <- config(...)
    response <- try(HEAD(url, config, user_agent("igblastr")), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(as.character(response))
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
        stop("no internet")
    config <- config(...)
    response <- try(GET(url, config, user_agent("igblastr"), query=query),
                    silent=TRUE)
    if (inherits(response, "try-error"))
        stop(as.character(response))
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
    if (!has_internet())
        stop("no internet")
    url <- paste0(dir_url, filename)
    destfile <- tempfile()
    code <- try(suppressWarnings(download.file(url, destfile, ...)),
                silent=TRUE)
    if (inherits(code, "try-error") || code != 0L)
        stop(wmsg("failed to download ", filename, " from ", dir_url))
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

### Returns a character vector with one string per key/val.
named_list_as_pretty_keyvals <- function(x, sep="; ")
{
    stopifnot(is.list(x))
    x_names <- names(x)
    stopifnot(!is.null(x_names))
    keys <- format(paste0(x_names, ": "))
    margin <- max(nchar(keys))
    x <- lapply(x, function(x) wmsg2(paste(x, collapse=sep), margin=margin))
    paste0(keys, as.character(x))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### system2e(), system3(), and related utils
###

### Returns TRUE if environment variable BLAST_USAGE_REPORT is set to true,
### and FALSE in **any** other case (i.e. if BLAST_USAGE_REPORT is not set,
### or is set to false or gibberish).
get_igblastr_usage_report_from_BLAST_USAGE_REPORT <- function()
{
    usage_report <- Sys.getenv("BLAST_USAGE_REPORT")
    nzchar(usage_report) && isTRUE(as.logical(toupper(usage_report)))
}

system2e <- function(...)
{
    igblastr_usage_report <- getOption("igblastr_usage_report")
    old_BLAST_USAGE_REPORT <- Sys.getenv("BLAST_USAGE_REPORT")
    on.exit(Sys.setenv(BLAST_USAGE_REPORT=old_BLAST_USAGE_REPORT))
    if (isTRUE(igblastr_usage_report)) {
        Sys.setenv(BLAST_USAGE_REPORT="true")
    } else {
        Sys.setenv(BLAST_USAGE_REPORT="false")
    }
    system2(...)
}

### To use on the result of 'try(system2e(..., stdout=TRUE, stderr=TRUE))'.
system_command_worked <- function(out)
{
    if (inherits(out, "try-error"))
        return(FALSE)
    status <- attr(out, "status")
    is.null(status) || isTRUE(all.equal(status, 0L))
}

system_command_works <- function(command, args=character())
{
    out <- try(suppressWarnings(system2e(command, args=args,
                                         stdout=TRUE, stderr=TRUE)),
               silent=TRUE)
    system_command_worked(out)
}

system3 <- function(command, outfile, errfile, args=character())
{
    status <- system2e(command, args=args, stdout=outfile, stderr=errfile)
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

has_perl <- function() system_command_works("perl", args="-v")

