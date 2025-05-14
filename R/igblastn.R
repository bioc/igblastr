### =========================================================================
### igblastn()
### -------------------------------------------------------------------------


### Not used at the moment.
.get_igblastr_tempdir <- function()
{
    dirpath <- file.path(tempdir(), "igblastr")
    if (!dir.exists(dirpath)) {
        if (file.exists(dirpath))
            stop(wmsg("Anomaly: '", dirpath, "' already exists ",
                      "as a file, not as a directory."))
        dir.create(dirpath)
    }
    dirpath
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_query()
###

.check_query_seq_lengths <- function(seqlens)
{
    stopifnot(is.integer(seqlens))
    seqids <- names(seqlens)
    stopifnot(!is.null(seqids))
    empty_idx <- which(seqlens == 0L)
    if (length(empty_idx) != 0L) {
        in1string <- paste(seqids[empty_idx], collapse=", ")
        stop(wmsg("the following sequences in 'query' are empty ",
                  "(showing seq ids): ", in1string))
    }
}

### Returns **absolute** path to FASTA file containing all the query sequences.
### The path is returned in a string with possibly the 'safe_to_remove'
### attribute on it indicating whether the caller can safely remove the
### FASTA file once it's done using it.
.normalize_character_query <- function(query)
{
    stopifnot(is.character(query))
    if (length(query) == 0L)
        stop(wmsg("character vector 'query' cannot be empty"))
    if (anyNA(query))
        stop(wmsg("character vector 'query' cannot contain NAs"))
    paths <- character(length(query))
    was_gz <- logical(length(query))
    for (i in seq_along(query)) {
        path <- file_path_as_absolute(query[[i]])
        .check_query_seq_lengths(fasta.seqlengths(path))
        if (has_suffix(path, ".gz")) {
	    was_gz[[i]] <- TRUE
            destpath <- tempfile("igblastn_query_", fileext=".fasta")
            gunzip(path, destname=destpath, overwrite=TRUE, remove=FALSE)
            path <- destpath
        }
        paths[[i]] <- path
    }
    if (length(query) == 1L) {
        ## 'paths' and 'was_gz' have length 1.
        if (was_gz)
            attr(paths, "safe_to_remove") <- TRUE
        return(paths)
    }
    path <- tempfile("igblastn_query_", fileext=".fasta")
    concatenate_files(paths, out=path)
    if (any(was_gz))
        unlink(paths[was_gz])
    attr(path, "safe_to_remove") <- TRUE
    path
}

### Returns **absolute** path to FASTA file containing all the query sequences.
### See comment for .normalize_character_query() above for the meaning
### of the 'safe_to_remove' attribute.
.normarg_query <- function(query)
{
    if (is.character(query))
        return(.normalize_character_query(query))
    if (is(query, "DNAStringSet")) {
        if (is.null(names(query)))
            stop(wmsg("DNAStringSet object 'query' must have names"))
        .check_query_seq_lengths(setNames(width(query), names(query)))
        path <- tempfile("igblastn_query_", fileext=".fasta")
        writeXStringSet(query, path)
        attr(path, "safe_to_remove") <- TRUE
        return(path)
    }
    stop(wmsg("'query' must be a character vector that contains ",
              "the paths to the input files (FASTA), or a named ",
              "DNAStringSet object"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_outfmt()
###

.stop_on_invalid_outfmt <- function()
{
    msg1 <- c("'outfmt' must be one of \"AIRR\", 3, 4, 7, or 19 (\"AIRR\" ",
              "is just an alias for 19), or a string describing a customized ",
              "format 7.")
    msg2 <- c("The string describing a customized format 7 must start ",
              "with a 7 followed by the desired hit table fields ",
              "(a.k.a. format specifiers) separated with spaces. ",
              "For example: \"7 std qseq sseq btop\". Note that 'std' ",
              "stands for 'qseqid sseqid pident length mismatch gapopen ",
              "gaps qstart qend sstart send evalue bitscore', which is ",
              "the default.")
    msg3 <- c("Use list_outfmt7_specifiers() to list all supported format ",
              "specifiers.")
    stop(wmsg(msg1), "\n  ", wmsg(msg2), "\n  ", wmsg(msg3))
}

### 'outfmt' is assumed to be a single string.
.check_customized_format_7 <- function(outfmt)
{
    if (!has_prefix(outfmt, "7 "))
        .stop_on_invalid_outfmt()
    user_specifiers <- substr(outfmt, 3L, nchar(outfmt))
    user_specifiers <- strsplit(user_specifiers, " ", fixed=TRUE)[[1L]]
    user_specifiers <- user_specifiers[nzchar(user_specifiers)]
    supported_specifiers <- c("std", names(list_outfmt7_specifiers()))
    invalid_specifiers <- setdiff(user_specifiers, supported_specifiers)
    if (length(invalid_specifiers) != 0L) {
        in1string <- paste(invalid_specifiers, collapse=", ")
        stop(wmsg("invalid format specifier(s): ", in1string))
    }
}

### 'outfmt' can be 3, 4, 7, 19, "AIRR", or a single string describing a
### customized format 7 e.g. "7 std qseq sseq btop".
### See .stop_on_invalid_outfmt() above for more information.
### Returns a single string.
.normarg_outfmt <- function(outfmt="AIRR")
{
    if (isSingleNumber(outfmt)) {
        if (!(outfmt %in% c(3, 4, 7, 19)))
            .stop_on_invalid_outfmt()
        outfmt <- as.character(as.integer(outfmt))
    } else if (isSingleNonWhiteString(outfmt)) {
        outfmt <- trimws2(outfmt)
        if (toupper(outfmt) == "AIRR") {
            outfmt <- "19"
        } else {
            if (!(outfmt %in% c("3", "4", "7", "19")))
                .check_customized_format_7(outfmt)
        }
    } else {
        .stop_on_invalid_outfmt()
    }
    outfmt
}

### Returns 3L, 4L, 7L, or 19L.
.extract_outfmt_nb <- function(outfmt)
{
    stopifnot(isSingleNonWhiteString(outfmt))
    outfmt_nb <- strsplit(trimws2(outfmt), " ", fixed=TRUE)[[1L]][[1L]]
    stopifnot(outfmt_nb %in% c("3", "4", "7", "19"))
    as.integer(outfmt_nb)
}

print.igblastn_raw_output <- function(x, ...) cat(x, sep="\n")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normarg_out()
###

### Must return an absolute path.
.normarg_out <- function(out)
{
    if (is.null(out))
        return(tempfile("igblastn_out_", fileext=".txt"))
    if (!isSingleNonWhiteString(out))
        stop(wmsg("'out' must be NULL or a single (non-empty) string"))
    dirpath <- dirname(out)
    if (!dir.exists(dirpath))
        stop(wmsg(dirpath, ": no such directory"))
    dirpath <- file_path_as_absolute(dirpath)
    out <- file.path(dirpath, basename(out))
    if (file.exists(out))
        stop(wmsg(out, ": file exists"))
    out
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .parse_igblastn_out()
###

.fix_AIRR_logical_cols <- function(AIRR_df)
{
    stopifnot(is.data.frame(AIRR_df))
    LOGICAL_COLS <- c("stop_codon", "vj_in_frame", "v_frameshift",
                      "productive", "rev_comp", "complete_vdj", "d_frame")
    for (colname in LOGICAL_COLS)
        AIRR_df[ , colname] <- as.logical(AIRR_df[ , colname])
    AIRR_df
}

### TODO: Parse output format 3 and 4.
.parse_igblastn_out <- function(out, outfmt_nb)
{
    stopifnot(isSingleNonWhiteString(out),
              isSingleInteger(outfmt_nb), outfmt_nb %in% c(3L, 4L, 7L, 19L))
    if (outfmt_nb == 19L) {
        ## Make sure to use 'tryLogical=FALSE'. By default, i.e.
        ## when 'tryLogical=TRUE', there's the risk that read.table()
        ## will erroneously decide to interpret some columns as logical!
        ## For example this can happen to column "d_sequence_alignment_aa"
        ## if it contains only single-letter amino acid sequences T (Thr)
        ## or F (Phe).
        AIRR_df <- read.table(out, header=TRUE, sep="\t",
                              tryLogical=FALSE, na.strings="")
        return(.fix_AIRR_logical_cols(AIRR_df))
    }
    out_lines <- readLines(out)
    if (outfmt_nb == 7L)
        return(parse_outfmt7(out_lines))
    warning(wmsg("parsing of igblastn output format ",
                 outfmt_nb, " is not supported yet"))
    class(out_lines) <- "igblastn_raw_output"
    out_lines
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastn()
###

.show_igblastn_command <- function(igblast_root, exe_args,
                                   show.in.browser=FALSE)
{
    igblastn_exe <- get_igblast_exe("igblastn", igblast_root=igblast_root)
    cmd <- c(igblastn_exe, exe_args)
    cmd_in_1string <- paste(cmd, collapse=" ")
    outfile <- if (show.in.browser)
               tempfile("igblastn_command_", fileext=".txt") else ""
    cat(cmd_in_1string, "\n", file=outfile, sep="")
    if (show.in.browser)
        display_local_file_in_browser(outfile)
    cmd  # returns the command in a character vector
}

.parse_warnings_or_errors <- function(lines, pattern)
{
    stopifnot(is.character(lines), isSingleNonWhiteString(pattern))
    keep_idx <- grep(pattern, lines, ignore.case=TRUE)
    msgs <- trimws2(lines[keep_idx])
    msgs[nzchar(msgs)]
}

.parse_and_issue_warnings <- function(stderr_file)
{
    warn_msgs <- .parse_warnings_or_errors(readLines(stderr_file), "warning:")
    for (msg in warn_msgs)
        warning(wmsg(msg))
}

.stop_on_igblastn_exe_error <- function(stderr_file)
{
    err_msgs <- .parse_warnings_or_errors(readLines(stderr_file), "error:")
    if (length(err_msgs) == 0L)  # could this ever happen?
        stop(wmsg("'igblastn' returned an unknown error"))
    err_msgs <- vapply(err_msgs, wmsg, character(1), USE.NAMES=FALSE)
    stop(paste(err_msgs, collapse="\n  "))
}

### The function calls setwd() before invoking the 'igblastn' executable so
### make sure that any file path passed thru 'exe_args' (e.g. the '-out' file
### path) is an **absolute** path.
.run_igblastn_exe <- function(igblast_root, exe_args)
{
    stopifnot(is.character(exe_args))
    igblastn_exe <- get_igblast_exe("igblastn", igblast_root=igblast_root)
    oldwd <- getwd()
    setwd(igblast_root)
    on.exit(setwd(oldwd))

    stderr_file <- tempfile("igblastn_stderr_", fileext=".txt")
    status <- system2(igblastn_exe, args=exe_args, stderr=stderr_file)
    .parse_and_issue_warnings(stderr_file)
    if (status != 0)
        .stop_on_igblastn_exe_error(stderr_file)
}

igblastn <- function(query, outfmt="AIRR",
                     germline_db_V="auto", germline_db_V_seqidlist=NULL,
                     germline_db_D="auto", germline_db_D_seqidlist=NULL,
                     germline_db_J="auto", germline_db_J_seqidlist=NULL,
                     organism="auto", c_region_db="auto",
                     auxiliary_data="auto",
                     ...,
                     out=NULL, parse.out=TRUE,
                     show.in.browser=FALSE, show.command.only=FALSE)
{
    if (is.null(out))
        on.exit(unlink(out))
    out <- .normarg_out(out)
    if (!isTRUEorFALSE(parse.out))
        stop(wmsg("'parse.out' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(show.in.browser))
        stop(wmsg("'show.in.browser' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(show.command.only))
        stop(wmsg("'show.command.only' must be TRUE or FALSE"))

    igblast_root <- get_igblast_root()
    query <- .normarg_query(query)
    if (isTRUE(attr(query, "safe_to_remove")))
        on.exit(unlink(query), add=TRUE)
    outfmt <- .normarg_outfmt(outfmt)
    outfmt_nb <- .extract_outfmt_nb(outfmt)
    cmd_args <- make_igblastn_command_line_args(
                     query, outfmt=outfmt,
                     germline_db_V=germline_db_V,
                     germline_db_V_seqidlist=germline_db_V_seqidlist,
                     germline_db_D=germline_db_D,
                     germline_db_D_seqidlist=germline_db_D_seqidlist,
                     germline_db_J=germline_db_J,
                     germline_db_J_seqidlist=germline_db_J_seqidlist,
                     organism=organism,
                     c_region_db=c_region_db,
                     auxiliary_data=auxiliary_data,
                     ...)
    exe_args <- make_exe_args(c(cmd_args, out=out))

    if (show.command.only) {
        ans <- .show_igblastn_command(igblast_root, exe_args,
                                      show.in.browser=show.in.browser)
        return(invisible(ans))
    }

    ## Run the 'igblastn' standalone executable included in IgBLAST.
    .run_igblastn_exe(igblast_root, exe_args)

    if (outfmt_nb == 19L) {
        if (show.in.browser || parse.out)
            AIRR_df <- .parse_igblastn_out(out, outfmt_nb)
        if (show.in.browser)
            display_data_frame_in_browser(AIRR_df)
        if (parse.out) {
            ans  <- tibble(AIRR_df)
        } else {
            ans <- readLines(out)
            class(ans) <- "igblastn_raw_output"
        }
        return(ans)
    }

    if (show.in.browser)
        display_local_file_in_browser(out)
    if (parse.out) {
        ans <- .parse_igblastn_out(out, outfmt_nb)
    } else {
        ans <- readLines(out)
        class(ans) <- "igblastn_raw_output"
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastn_help()
###

igblastn_help <- function(long.help=FALSE, show.in.browser=FALSE)
{
    if (!isTRUEorFALSE(long.help))
        stop(wmsg("'long.help' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(show.in.browser))
        stop(wmsg("'show.in.browser' must be TRUE or FALSE"))

    igblast_root <- get_igblast_root()
    igblastn_exe <- get_igblast_exe("igblastn", igblast_root=igblast_root)
    exe_args <- if (long.help) "-help" else "-h"

    oldwd <- getwd()
    setwd(igblast_root)
    on.exit(setwd(oldwd))
    outfile <- if (show.in.browser)
               tempfile("igblastn_help_", fileext=".txt") else ""
    status <- system2(igblastn_exe, args=exe_args, stdout=outfile)
    if (status != 0)
        stop(wmsg("'igblastn' returned an error"))
    if (show.in.browser)
        display_local_file_in_browser(outfile)
}

