### =========================================================================
### Low-level file manipulation
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nuke_file()
###

### Use responsibly!
### NOTE: Contrary to what its man page says, 'unlink(path, recursive=TRUE)'
### can still silently fail and return 0 if it cannot delete the
### supplied 'path'. This happens for example on Linux if 'path' points to
### a file or directory located in a directory witout x access (e.g.
### chmod 400 or 600). The "workaround" below that consists in checking
### for the existence of 'path' does not help either because in that case
### file.exists() can't be trusted either: it silently returns FALSE even
### though 'path' exists.
### TODO: Report this issue to R core.
nuke_file <- function(path)
{
    stopifnot(isSingleNonWhiteString(path))
    ret <- unlink(path, recursive=TRUE, force=TRUE)
    ## file.exists() "workaround" does not help. See NOTE above.
    if (ret != 0L || file.exists(path))
        stop(wmsg("failed to delete '", path, "'"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rename_file()
###

### Unfortunately, when 'from' is a directory, 'file.copy(from, to)'
### insists that 'to' should also be an existing directory, so that it
### can create a copy of 'from' **inside** it. So the path to the copied
### folder will always be 'file.path(to, basename(from))'. This means that,
### unlike with Unix 'cp' command, there's no way to copy **and** rename.
### This .copy_dir() helper tries to remedy this by mimicking the semantic
### of Unix 'cp' command. It proceeds in two steps: copy and rename.
.copy_dir <- function(from, to, copy.mode=TRUE, copy.date=FALSE)
{
    stopifnot(isSingleNonWhiteString(from), isSingleNonWhiteString(to))
    if (!dir.exists(from))
        stop(wmsg(from, ": not a directory"))

    errmsg <- c("failed to copy '", from, "' to '", to, "'")

    ## copy_dir_to_dir() expected to always be used on existing directories.
    copy_dir_to_dir <- function(.from, .to) {
        stopifnot(dir.exists(.from), dir.exists(.to))
        ## We do not overwrite the destination folder.
        dest <- file.path(.to, basename(.from))
        if (file.exists(dest))
            stop(wmsg(dest, ": file exists"))
        ok <- file.copy(.from, .to, recursive=TRUE,
                        copy.mode=copy.mode, copy.date=copy.date)
        if (!ok) stop(wmsg(errmsg))
    }
    ## In the context of the .copy_dir() function, rename_dir() should never
    ## be used to rename a directory across file systems.
    rename_dir <- function(.from, .to) {
        stopifnot(dir.exists(.from))
        ok <- file.rename(.from, .to)
        if (!ok) stop(wmsg(errmsg))
    }

    if (dir.exists(to))
        return(copy_dir_to_dir(from, to))

    to_parent <- dirname(to)
    if (!dir.exists(to_parent))
        stop(wmsg(to_parent, ": directory does not exist"))

    if (basename(from) == basename(to)) ## We're lucky!
        return(copy_dir_to_dir(from, to_parent))

    middle <- file.path(to_parent, basename(from))
    if (!file.exists(middle)) {
        ## We're still somewhat lucky!
        on.exit(nuke_file(middle))
        copy_dir_to_dir(from, to_parent)
        ## rename_dir() unlikely to fail (if we were able to copy to 'to_parent'
        ## in the first place, then we should be able to rename a subfolder
        ## of 'to_parent').
        rename_dir(middle, to)
        return(invisible(NULL))
    }

    ## We ran out of luck. Last chance is to copy to a temporary hidden
    ## subfolder of 'to_parent' and then to rename.
    repeat {
        temp_to_basename <- paste0(".", basename(tempfile()))
        temp_to <- file.path(to_parent, temp_to_basename)
        if (!file.exists(temp_to)) break
    }
    on.exit(nuke_file(temp_to))
    ok <- dir.create(temp_to)
    if (!ok)
        stop(wmsg(errmsg))
    copy_dir_to_dir(from, temp_to)
    ## rename_dir() unlikely to fail (if we were able to copy to 'to_parent'
    ## in the first place, then we should be able to rename a subfolder
    ## of 'to_parent').
    middle <- file.path(temp_to, basename(from))
    rename_dir(middle, to)
}

### Brute force version of file.rename().
### Renames 'from' to 'to'.
### 'from' can point to a file or directory.
### Raises an error if 'to' exists, unless 'replace' is TRUE in which
### case 'to' will be nuked.
### If 'to' does not exist, then its parent directory will be created.
### Works across file systems.
rename_file <- function(from, to, replace=FALSE)
{
    stopifnot(isSingleNonWhiteString(from),
              isSingleNonWhiteString(to),
              isTRUEorFALSE(replace))
    if (!file.exists(from))
        stop(wmsg(from, ": no such file or directory"))
    if (file.exists(to)) {
        if (!replace)
            stop(wmsg("'", to, "' already exists"))
        nuke_file(to)
    } else {
        to_parent <- dirname(to)
        if (!dir.exists(to_parent)) {
            if (file.exists(to_parent))
                stop(wmsg("cannot create dir '", to_parent, "' ",
                          "(path exists as a non-dir file ",
                          "so I won't try to delete it)"))
            ok <- suppressWarnings(dir.create(to_parent, recursive=TRUE))
            if (!ok)
                stop(wmsg("failed to create dir '", to_parent, "'"))
        }
    }
    ## Can fail if 'from' or 'to' are not on the same file system, or if
    ## one of them is located in a non-writable directory.
    ok <- suppressWarnings(file.rename(from, to))
    if (!ok) {
        if (dir.exists(from)) {
            suppressWarnings(.copy_dir(from, to, copy.date=TRUE))
        } else {
            ## 'from' is a file that is not a directory.
            ## Can fail if 'from' is not readable or if 'to' is located in
            ## a non-writable directory.
            ok <- suppressWarnings(file.copy(from, to, copy.date=TRUE))
            if (!ok)
                stop(wmsg("failed to rename '", from, "' to '", to, "'"))
        }
        nuke_file(from)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### remove_empty_dir()
###

### remove_empty_dir(path, parents=TRUE) same as 'rmdir --parents <path>'
### except that it doesn't return an error if 'path' is not empty.
### Returns the first anscestor that is not empty.
remove_empty_dir <- function(path, parents=FALSE)
{
    stopifnot(isSingleNonWhiteString(path), isTRUEorFALSE(parents))
    if (!dir.exists(path)) {
        if (file.exists(path))
            stop(wmsg(path, ": not a directory"))
        stop(wmsg(path, ": no such directory"))
    }
    if (path == ".")
        stop(wmsg("'path' cannot be \".\""))
    empty <- length(list.files(path, all.files=TRUE, no..=TRUE)) == 0L
    if (empty) {
        ret <- unlink(path, recursive=TRUE)
        if (ret != 0L || file.exists(path))
            stop(wmsg("failed to remove empty directory ", path))
        path <- dirname(path)
        if (path != "." && parents)
            path <- remove_empty_dir(path, parents=TRUE)
    } else {
        if (!parents)
            stop(wmsg(path, ": directory not empty"))
    }
    path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### remove_hidden_files()
###

### Does not recursively search hidden files i.e. it only removes the
### hidden files (and directories if 'include.hidden.dirs=TRUE') that
### are located **directly** under 'path'.
remove_hidden_files <- function(path=".", include.hidden.dirs=FALSE)
{
    stopifnot(isSingleNonWhiteString(path), dir.exists(path),
              isTRUEorFALSE(include.hidden.dirs))
    hidden_files <- list.files(path, pattern="^\\.", all.files=TRUE,
                               full.names=TRUE, no..=TRUE)
    unlink(hidden_files, recursive=include.hidden.dirs, force=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### copy_files_to_dir()
###

### The semantic of 'file.copy(files, destdir, overwrite=FALSE)' is not
### what we want (it simply skips copying the files for which the
### destination file exists). What we want is a copy everyting or nothing
### semantic (i.e. atomicity), which is what copy_files_to_dir() achieves.
copy_files_to_dir <- function(files, destdir, overwrite=FALSE)
{
    stopifnot(is.character(files),
              isSingleNonWhiteString(destdir), dir.exists(destdir),
              isTRUEorFALSE(overwrite))
    if (!overwrite) {
        filenames <- basename(files)
        destfiles <- file.path(destdir, filenames)
        already_there <- file.exists(destfiles)
        if (any(already_there)) {
            in1string <- paste0(filenames[already_there], collapse=", ")
            stop(wmsg("The following files are already present in 'destdir': ",
                      in1string, "."),
                 "\n  ",
                 wmsg("Use 'overwrite=TRUE' to replace them."))
        }
    }
    ## Note that file.copy() has its own way of handling 'overwrite=FALSE'.
    file.copy(files, destdir, overwrite=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### concatenate_files()
###

concatenate_files <- function(files, out=stdout(), n=50000L)
{
    stopifnot(is.character(files))
    if (is.character(out)) {
        out <- file(out, "wb")
        on.exit(close(out))
    }
    for (f in files) {
        con <- file(f, "rb")
        while (TRUE) {
            bytes <- readBin(con, what=raw(), n=n)
            if (length(bytes) == 0L)
                break
            writeBin(bytes, out)
        }
        close(con)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### display_local_file_in_browser()
### display_data_frame_in_browser()
###

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
### list_fasta_files()
###

### Returns the full names by default! (unlike list.files())
list_fasta_files <- function(path, pattern="\\.fasta$", full.names=TRUE,
                             recursive=FALSE)
{
    stopifnot(isSingleNonWhiteString(path), dir.exists(path),
              isSingleNonWhiteString(pattern),
              isTRUEorFALSE(full.names),
              isTRUEorFALSE(recursive))
    list.files(path, pattern=pattern, full.names=full.names,
               recursive=recursive)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fasta_files_as_one_uncompressed_file()
###

### Returns **absolute** path to the **uncompressed** FASTA file
### containing the sequences from all the input FASTA files.
### The path is returned in a string with possibly the 'safe_to_remove'
### attribute on it indicating whether the caller can safely remove the
### FASTA file once it's done using it.
fasta_files_as_one_uncompressed_file <- function(fasta_files, varname)
{
    stopifnot(is.character(fasta_files))
    if (length(fasta_files) == 0L)
        stop(wmsg("character vector '", varname, "' cannot be empty"))
    if (anyNA(fasta_files))
        stop(wmsg("character vector '", varname, "' cannot contain NAs"))
    temp_file_pattern <- paste0(varname, "_")
    paths <- character(length(fasta_files))
    was_gz <- logical(length(fasta_files))
    for (i in seq_along(fasta_files)) {
        path <- file_path_as_absolute(fasta_files[[i]])
        check_seqlens(fasta.seqlengths(path), varname)
        if (has_suffix(path, ".gz")) {
            was_gz[[i]] <- TRUE
            destpath <- tempfile(temp_file_pattern, fileext=".fasta")
            gunzip(path, destname=destpath, overwrite=TRUE, remove=FALSE)
            path <- destpath
        }
        paths[[i]] <- path
    }
    if (length(fasta_files) == 1L) {
        ## 'paths' and 'was_gz' have length 1.
        if (was_gz)
            attr(paths, "safe_to_remove") <- TRUE
        return(paths)
    }
    path <- tempfile(temp_file_pattern, fileext=".fasta")
    concatenate_files(paths, out=path)
    if (any(was_gz))
        unlink(paths[was_gz])
    attr(path, "safe_to_remove") <- TRUE
    path
}

