### =========================================================================
### Low-level manipulation of IgBLAST data
### -------------------------------------------------------------------------
###


.path_to_igdata_store <- function()
{
    system.file(package="igblastr", "extdata", "igdata_store", mustWork=TRUE)
}

.path_to_igdata <- function(which=c("live", "original"))
{
    which <- match.arg(which)
    if (which == "live") {
        igdata_dir <- igblastr_cache(LIVE_IGDATA)
    } else {
        igdata_dir <- get_igblast_root()
    }
    igdata_dir
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### install_igdata_in_igblast_root()
### check_igdata_store()
###

### Only call this to repair an IgBLAST installation that is missing
### the internal_data/ or optional_file/ folders. See README.txt in
### igblastr/inst/extdata/igdata_store/ for more information.
install_igdata_in_igblast_root <- function(igblast_root)
{
    if (!isSingleNonWhiteString(igblast_root))
        stop(wmsg("'igblast_root' must be a single (non-empty) string"))
    if (!dir.exists(igblast_root))
        stop(wmsg("directory '", igblast_root, "' does not exist"))
    igdata_store <- .path_to_igdata_store()
    subdirs <- list.dirs(igdata_store, recursive=FALSE)
    file.copy(subdirs, igblast_root, recursive=TRUE)
}

### Compares static data found in igdata_store with data found under the
### specified igblast_root, and raises an error if they differ.
### Only used in unit tests at the moment. Relies on Unix command 'diff'
### so not supported on Windows.
check_igdata_store <- function(igblast_root)
{
    if (!isSingleNonWhiteString(igblast_root))
        stop(wmsg("'igblast_root' must be a single (non-empty) string"))
    if (!dir.exists(igblast_root))
        stop(wmsg("directory '", igblast_root, "' does not exist"))
    igdata_store <- .path_to_igdata_store()
    subdirs <- list.dirs(igdata_store, recursive=FALSE)
    for (subdir in subdirs) {
        subdir2 <- file.path(igblast_root, basename(subdir))
        if (!dir.exists(subdir2))
            stop(wmsg("invalid IgBLAST installation at '", igblast_root, "': ",
                      "directory has no '", basename(subdir), "' subdirectory"))
        args <- c("-r", "--brief", subdir, subdir2)
        out <- suppressWarnings(system2("diff", args, stdout=TRUE))
        status <- attr(out, "status")
        if (!is.null(status) || length(out) != 0L)
            warning(wmsg("content of ", subdir, " and ", subdir2, " differ"))
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reset_live_igdata()
###

reset_live_igdata <-
    function(subdirs=c("all", "internal_data", "optional_file"))
{
    subdirs <- match.arg(subdirs)
    if (subdirs == "all")
        subdirs <- c("internal_data", "optional_file")
    igdata_store <- .path_to_igdata_store()
    live_igdata <- igblastr_cache(LIVE_IGDATA)
    if (!dir.exists(live_igdata))
        dir.create(live_igdata, recursive=TRUE)
    for (subdir in subdirs) {
        nuke_file(file.path(live_igdata, subdir))
        from <- file.path(igdata_store, subdir)
        file.copy(from, live_igdata, recursive=TRUE)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### update_live_igdata()
###

.LATEST_IGBLAST_DATA_FTP_DIR <-
    "ftp.ncbi.nih.gov/blast/executables/igblast/release/patch/"

.summarize_update <- function(new_files, updated_files, check.only=FALSE)
{
    stopifnot(is.character(new_files), is.character(updated_files),
              isTRUEorFALSE(check.only))
    make_subject_verb <- function(n) {
        subject <- if (n == 1L) "file" else "files"
        verb <- if (check.only) "can be" else if (n == 1L) "was" else "were"
        paste0(subject, " ", verb)
    }
    n <- length(new_files)
    if (n != 0L) {
        what <- make_subject_verb(n)
        in1string <- paste(new_files, collapse=", ")
        message("The following new ", what, " installed: ", in1string)
    }
    n <- length(updated_files)
    if (n != 0L) {
        what <- make_subject_verb(n)
        in1string <- paste(updated_files, collapse=", ")
        message("The following ", what, " updated: ", in1string)
    }
}

### Returns the list of files that got actually updated i.e. that got
### replaced with a different file.
### WARNING: Assumes that the internal_data/ and optional_file/ folders
### at .LATEST_IGBLAST_DATA_FTP_DIR contain only files. Will break if
### they contain subfolders. Note that we don't expect this to happen
### for the optional_file/ folder but it WILL happen when the NCBI folks
### start populating internal_data/ with stuff (this subfolder is empty
### at the moment). When this happens, calling
###     .download_latest_igblast_data("internal_data")
### will break!
### FIXME: Make .download_latest_igblast_data() work on a hierarchy of
### folders.
.download_latest_igblast_data <-
    function(subdir=c("internal_data", "optional_file"), check.only=FALSE, ...)
{
    subdir <- match.arg(subdir)
    if (!isTRUEorFALSE(check.only))
        stop(wmsg("'check.only' must be TRUE or FALSE"))
    message("Checking available updates for '", subdir, "' ... ",
            appendLF=FALSE)
    ftp_dir <- paste0(.LATEST_IGBLAST_DATA_FTP_DIR, subdir, "/")
    listing <- try(suppressWarnings(list_ftp_dir(ftp_dir)), silent=TRUE)
    if (inherits(listing, "try-error"))
        stop(wmsg("Cannot open URL '", ftp_dir, "'. ",
                  "Are you connected to the internet?"))

    igdata_subdir <- file.path(igblastr_cache(LIVE_IGDATA), subdir)
    stopifnot(dir.exists(igdata_subdir))

    new_files <- updated_files <- character(0)
    for (datafile in listing) {
        downloaded_file <- download_as_tempfile(ftp_dir, datafile, ...)
        current_file <- file.path(igdata_subdir, datafile)
        if (file.exists(current_file)) {
            if (md5sum(downloaded_file) == md5sum(current_file))
                next
            updated_files <- c(updated_files, datafile)
        } else {
            new_files <- c(new_files, datafile)
        }
        if (!check.only)
            rename_file(downloaded_file, current_file, replace=TRUE)
    }
    message("ok")

    current_time <- format(Sys.time(), "%Y-%m-%d %H:%M")
    path <- file.path(igdata_subdir, "LAST_CHECKED")
    writeLines(current_time, path)
    if (length(new_files) + length(updated_files) == 0L) {
        message("No new updates found.")
    } else {
        if (!check.only) {
            path <- file.path(igdata_subdir, "LAST_UPDATED")
            writeLines(current_time, path)
        }
        .summarize_update(new_files, updated_files, check.only=check.only)
    }

    list(new=new_files, updated=updated_files)
}

update_live_igdata <- function(check.only=FALSE)
{
    ## Skipping internal_data/ for now. The folder at
    ## .LATEST_IGBLAST_DATA_FTP_DIR is currently empty, and,
    ## the day the NCBI folks start populating it with stuff,
    ## .download_latest_igblast_data("internal_data") will break.
    ## See WARNING/FIXME for .download_latest_igblast_data() above.
    #.download_latest_igblast_data("internal_data", quiet=TRUE)
    .download_latest_igblast_data("optional_file", check.only=check.only,
                                  quiet=TRUE)
    invisible(NULL)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_aux_data_md5sum_df()
###

.compute_aux_files_md5sums <- function(which=c("live", "original"))
{
    aux_dir <- file.path(.path_to_igdata(which), "optional_file")
    aux_files <- list.files(aux_dir, pattern="\\.(aux|frame)$")
    setNames(md5sum(file.path(aux_dir, aux_files)), aux_files)
}

.make_aux_data_md5sum_df <- function()
{
    live_md5sums <- .compute_aux_files_md5sums("live")
    orig_md5sums <- .compute_aux_files_md5sums("original")
    files <- sort(union(names(live_md5sums), names(orig_md5sums)))
    ans <- data.frame(file=files,
                      live=unname(live_md5sums[files]),
                      original=unname(orig_md5sums[files]))
    class(ans) <- c("aux_data_md5sum_df", class(ans))
    ans
}

.aux_data_md5sum_df_as_character <- function(x)
{
    EXPECTED_COLNAMES <- c("file", "live", "original")
    stopifnot(is.data.frame(x), identical(colnames(x), EXPECTED_COLNAMES))
    note <- character(nrow(x))
    note[x$live != x$original] <- " (updated)"
    note[is.na(x$original)] <- " (new)"
    ans <- character(3L * nrow(x))
    ans[3L * seq_len(nrow(x)) - 2L] <- paste0("- file '", x$file, "'")
    ans[3L * seq_len(nrow(x)) - 1L] <- paste0("    live:     ", x$live, note)
    ans[3L * seq_len(nrow(x))     ] <- paste0("    original: ", x$original)
    ans
}

print.aux_data_md5sum_df <- function(x, margin="", ...)
{
    y <- .aux_data_md5sum_df_as_character(x)
    cat(paste0(margin, y), sep="\n")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### time_since_live_igdata_last_checked()
### igdata_info()
###

.get_last_checked <- function(igdata_subdir)
{
    stopifnot(isSingleNonWhiteString(igdata_subdir), dir.exists(igdata_subdir))
    path <- file.path(igdata_subdir, "LAST_CHECKED")
    if (file.exists(path)) readLines(path) else NA_character_
}

.get_last_updated <- function(igdata_subdir)
{
    stopifnot(isSingleNonWhiteString(igdata_subdir), dir.exists(igdata_subdir))
    path <- file.path(igdata_subdir, "LAST_UPDATED")
    if (file.exists(path)) readLines(path) else NA_character_
}

.how_old <- function(datetime, units="days")
{
    stopifnot(isSingleNonWhiteString(datetime))
    as.double(difftime(Sys.time(), as.POSIXct(datetime), units=units))
}

### Used in .onLoad() hook.
### By default, time is returned in number of days.
time_since_live_igdata_last_checked <- function(units="days")
{
    igdata_subdir <- file.path(igblastr_cache(LIVE_IGDATA), "optional_file")
    last_checked <- .get_last_checked(igdata_subdir)
    if (is.na(last_checked))
        return(Inf)
    .how_old(last_checked, units=units)
}

.annotate_how_old <- function(dt)
{
    msg <- paste0(sprintf("%.2f", dt), " days ago")
    if (dt > 30)
        msg <- paste0(msg, " --> time to run 'update_live_igdata()' again")
    msg
}

.annotate_last_checked <- function(last_checked)
{
    if (is.na(last_checked))
        return(paste0("Not checked yet --> ",
                      "run 'update_live_igdata()' to check ",
                      "for new updates available at NCBI."))
    dt <- .how_old(last_checked)
    paste0(last_checked, " (", .annotate_how_old(dt), ")")
}

igdata_info <- function()
{
    live_igdata <- .path_to_igdata("live")
    igdata_subdir <- file.path(live_igdata, "optional_file")
    last_checked <- .get_last_checked(igdata_subdir)
    last_updated <- .get_last_updated(igdata_subdir)
    if (is.na(last_updated))
        last_updated <- "Not updated yet."
    ans <- list(
        live_igdata=live_igdata,
        last_checked=.annotate_last_checked(last_checked),
        last_updated=last_updated,
        original_igdata=.path_to_igdata("original"),
        aux_data_md5sums=.make_aux_data_md5sum_df()
    )
    class(ans) <- "igdata_info"
    ans
}

print.igdata_info <- function(x, ...)
{
    stopifnot(is.list(x))
    md5sum_df <- x$aux_data_md5sums
    if (!is.null(md5sum_df))
        x$aux_data_md5sums <- NULL
    keyvals <- named_list_as_pretty_keyvals(x)
    if (!is.null(md5sum_df)) {
        y <- .aux_data_md5sum_df_as_character(md5sum_df)
        keyval <- paste0("aux_data_md5sums:\n", paste0("  ", y, collapse="\n"))
        keyvals <- c(keyvals, keyval)
    }
    cat(keyvals, sep="\n")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_igblast_auxiliary_data()
###

get_igblast_auxiliary_data <- function(organism, which=c("live", "original"))
{
    organism <- normalize_igblast_organism(organism)
    which <- match.arg(which)
    aux_dir <- file.path(.path_to_igdata(which), "optional_file")
    aux_file <- file.path(aux_dir, paste0(organism, "_gl.aux"))
    if (!file.exists(aux_file))
        stop(wmsg("no auxiliary data found in ", aux_dir, " for ", organism))
    aux_file
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_igblast_auxiliary_data()
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
        sseqid            =           m[ , 1L],
        coding_frame_start=as.integer(m[ , 2L]),
        chaintype         =           m[ , 3L],
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
load_igblast_auxiliary_data <- function(organism, which=c("live", "original"))
{
    which <- match.arg(which)
    filepath <- get_igblast_auxiliary_data(organism, which)
    m <- .read_jagged_table_as_character_matrix(filepath)
    .make_auxiliary_data_df_from_matrix(m, filepath)
}

