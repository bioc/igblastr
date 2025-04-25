### =========================================================================
### Low-level utilities related to germline and C-region db manipulation
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


VDJ_REGION_TYPES <- c("V", "D", "J")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### disambiguate_fasta_seqids()
###

### Similar to base::make.unique() but mangles with suffixes made of
### lowercase letters.
.make_pool_of_suffixes <- function(min_pool_size)
{
    max_pool_size <- (length(letters)**8 - 1) / (length(letters) - 1) - 1
    if (min_pool_size > max_pool_size)
        stop(wmsg("too many duplicate seq ids"))
    ans <- character(0)
    for (i in 1:7) {
        ans <- c(ans, mkAllStrings(letters, i))
        if (length(ans) >= min_pool_size)
            return(ans)
    }
    ## Should never happen because we checked for this condition earlier (see
    ## above).
    stop(wmsg("too many duplicate seq ids"))
}

.make_unique_seqids <- function(seqids)
{
    stopifnot(is.character(seqids))
    if (length(seqids) <= 1L)
        return(seqids)
    oo <- order(seqids)
    seqids2 <- seqids[oo]
    ir <- IRanges(1L, runLength(Rle(seqids2)))
    pool_of_suffixes <- .make_pool_of_suffixes(max(width(ir)))
    suffixes <- extractList(pool_of_suffixes, ir)  # CharacterList
    suffixes[lengths(suffixes) == 1L] <- ""
    unlist(suffixes, use.names=FALSE)
    seqids2 <- paste0(seqids2, unlist(suffixes, use.names=FALSE))
    ans <- seqids2[S4Vectors:::reverseIntegerInjection(oo, length(oo))]
    setNames(ans, names(seqids))
}

### In-place replacement!
disambiguate_fasta_seqids <- function(filepath)
{
    stopifnot(isSingleNonWhiteString(filepath))
    fasta_lines <- readLines(filepath)
    header_idx <- grep("^>", fasta_lines)
    header_lines <- fasta_lines[header_idx]
    if (anyDuplicated(header_lines)) {
        fasta_lines[header_idx] <- .make_unique_seqids(header_lines)
        writeLines(fasta_lines, filepath)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_version_file()
###

read_version_file <- function(dirpath)
{
    stopifnot(isSingleNonWhiteString(dirpath))
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
### sort_db_names()
###

### The sorting is guaranteed to be the same everywhere. In particular it's
### guaranteed to put the db names prefixed with an underscore first.
sort_db_names <- function(db_names, decreasing=FALSE)
{
    stopifnot(is.character(db_names))
    ok <- has_prefix(db_names, "_")
    ## We set LC_COLLATE to C so:
    ## 1. sort() gives the same output whatever the platform or country;
    ## 2. sort() will behave the same way when called in the context
    ##    of 'R CMD build' or 'R CMD check' (both set 'R CMD check'
    ##    LC_COLLATE to C when building the vignette or running the tests)
    ##    vs when called in the context of an interactive session;
    ## 3. sort() is about 4x faster vs when LC_COLLATE is set to en_US.UTF-8.
    prev_locale <- Sys.getlocale("LC_COLLATE")
    Sys.setlocale("LC_COLLATE", "C")
    on.exit(Sys.setlocale("LC_COLLATE", prev_locale))
    ans1 <- sort(db_names[ok], decreasing=decreasing)
    ans2 <- sort(db_names[!ok], decreasing=decreasing)
    c(ans1, ans2)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_db_fasta_file()
###

### Note that the returned path is NOT guaranteed to exist.
get_db_fasta_file <- function(db_path, region_type=c(VDJ_REGION_TYPES, "C"))
{
    stopifnot(isSingleNonWhiteString(db_path))
    region_type <- match.arg(region_type)
    file.path(db_path, paste0(region_type, ".fasta"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .tabulate_dbs_by_region_type()
###

### Returns a named integer vector with 'region_types' as names.
.tabulate_db_by_region_type <- function(db_path, region_types)
{
    vapply(region_types,
        function(region_type) {
            fasta_file <- get_db_fasta_file(db_path, region_type)
            length(fasta.seqlengths(fasta_file))
        }, integer(1))
}

### Returns an integer matrix with 1 row per db and 1 col per region type.
.tabulate_dbs_by_region_type <- function(dbs_path, db_names, region_types)
{
    all_counts <- lapply(db_names,
        function(db_name) {
            db_path <- file.path(dbs_path, db_name)
            .tabulate_db_by_region_type(db_path, region_types)
        })
    data <- unlist(all_counts, use.names=FALSE)
    if (is.null(data))
        data <- integer(0)
    matrix(data, ncol=length(region_types), byrow=TRUE,
                 dimnames=list(NULL, region_types))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .tabulate_db_by_group()
### .tabulate_c_region_db_by_locus()
###

### All prefixes must have the same length.
.tabulate_gene_names_by_prefix <- function(gene_names, prefixes)
{
    stopifnot(is.character(gene_names),
              is.character(prefixes), length(prefixes) >= 1L)
    nc <- nchar(prefixes)
    stopifnot(all(nc == nc[[1L]]))
    gene_prefixes <- substr(gene_names, 1L, nc)
    m <- match(gene_prefixes, prefixes)
    if (anyNA(m)) {
        in1string <- paste0(prefixes, collapse=", ")
        stop(wmsg("not all gene names start with one of the ",
                  "following prefixes: ", in1string))
    }
    setNames(tabulate(m, length(prefixes)), prefixes)
}

### Returns an integer matrix, typically with some NAs.
.tabulate_db_by_group <- function(db_path, groups)
{
    stopifnot(isSingleNonWhiteString(db_path),
              is.character(groups), all(nchar(groups, type="width") == 4L))
    prefixes <- substr(groups, 1L, 3L)
    suffixes <- substr(groups, 4L, 4L)
    loci <- unique(prefixes)
    region_types <- unique(suffixes)
    ans <- matrix(NA_integer_, nrow=length(loci),
                               ncol=length(region_types),
                               dimnames=list(loci, region_types))
    Mindex <- cbind(prefixes, suffixes)
    ans[Mindex] <- 0L
    for (region_type in region_types) {
        fasta_file <- get_db_fasta_file(db_path, region_type)
        gene_names <- names(fasta.seqlengths(fasta_file))
        counts <- .tabulate_gene_names_by_prefix(gene_names, groups)
        stopifnot(identical(names(counts), groups))
        ok <- has_suffix(groups, region_type) | counts == 0L
        if (!all(ok)) {
            in1string <- paste(groups[!ok], collapse=", ")
            warning(wmsg("some gene names in '" , fasta_file, "' ",
                         "unexpectedly belong to groups ", in1string))
        }
        ans[Mindex] <- ans[Mindex] + counts
    }
    ans
}

### Returns a named integer vector with 'loci' as names.
.tabulate_c_region_db_by_locus <- function(db_path, loci)
{
    stopifnot(isSingleNonWhiteString(db_path), is.character(loci))
    fasta_file <- get_db_fasta_file(db_path, "C")
    gene_names <- names(fasta.seqlengths(fasta_file))
    .tabulate_gene_names_by_prefix(gene_names, loci)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_dbs()
###

.IG_LOCI <- paste0("IG", c("H", "K", "L"))
.TR_LOCI <- paste0("TR", c("A", "B", "G", "D"))

### Group names are formed by concatenating a locus name (e.g. IGH or TRB)
### and a region type (e.g. V).

.VJ_REGION_TYPES <- VDJ_REGION_TYPES[-2L]

.IG_GERMLINE_GROUPS <- c(paste0("IGH", VDJ_REGION_TYPES),
                         paste0("IGK", .VJ_REGION_TYPES),
                         paste0("IGL", .VJ_REGION_TYPES))

.TR_GERMLINE_GROUPS <- c(paste0("TRA", .VJ_REGION_TYPES),
                         paste0("TRB", VDJ_REGION_TYPES),
                         paste0("TRG", .VJ_REGION_TYPES),
                         paste0("TRD", VDJ_REGION_TYPES))

### Returns a named list with 1 list element per db.
.make_long_listing_for_germline_dbs <- function(dbs_path, db_names, groups)
{
    lapply(setNames(db_names, db_names),
        function(db_name) {
            db_path <- file.path(dbs_path, db_name)
            .tabulate_db_by_group(db_path, groups)
        })
}

### Returns a named list with 1 list element per db.
.make_long_listing_for_c_region_dbs <- function(dbs_path, db_names, loci)
{
    lapply(setNames(db_names, db_names),
        function(db_name) {
            db_path <- file.path(dbs_path, db_name)
            .tabulate_c_region_db_by_locus(db_path, loci)
        })
}

### 'long.listing' is ignored when 'names.only' is TRUE.
list_dbs <- function(dbs_path, what=c("germline", "C-region"),
                     builtin.only=FALSE,
                     names.only=FALSE, long.listing=FALSE)
{
    stopifnot(isSingleNonWhiteString(dbs_path), dir.exists(dbs_path))
    what <- match.arg(what)
    if (!isTRUEorFALSE(builtin.only))
        stop(wmsg("'builtin.only' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(names.only))
        stop(wmsg("'names.only' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(long.listing))
        stop(wmsg("'long.listing' must be TRUE or FALSE"))
    ## Excluding the 'USING' file for backward compatibility reasons.
    ## See NOTE preceding '.DB_IN_USE_cache' below in this file.
    db_names <- setdiff(list.files(dbs_path), "USING")
    if (builtin.only)
        db_names <- db_names[has_prefix(db_names, "_")]
    db_names <- sort_db_names(db_names)
    if (names.only)
        return(db_names)
    if (!long.listing) {
        region_types <- if (what == "germline") VDJ_REGION_TYPES else "C"
        basic_stats <- .tabulate_dbs_by_region_type(dbs_path, db_names,
                                                    region_types)
        return(data.frame(db_name=db_names, basic_stats))
    }
    if (what == "germline") {
        .make_long_listing_for_germline_dbs(dbs_path, db_names,
                                            .IG_GERMLINE_GROUPS)
    } else {
        .make_long_listing_for_c_region_dbs(dbs_path, db_names, .IG_LOCI)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_db_in_use()
### set_db_in_use()
###

### NOTE: Early versions of igblastr (prior to submission to Bioconductor)
### were recording the selection of the germline and C-region dbs in a
### persistent manner in a file located in the 'dbs_path' folder called 'USING'.
### This was obviously a very bad idea because it meant that:
### - If the user was using igblastr in more than one R session, then all
###   sessions were forced to use the same dbs.
### - Changing the selection in one session would change it for all other
###   sessions, which was hugely problematic!
### So, starting with igblastr 0.99.0, the db selection is recorded in memory
### (in the '.DB_IN_USE_cache' environment below) rather than in a file. This
### means that it is now remembered for the duration of the session only, and
### does not persist across sessions. It also means that the user now needs to
### specify the selection at the beginning of each session, which is actually
### a good thing!
.DB_IN_USE_cache <- new.env(parent=emptyenv())

### Returns "" if no db is currently in use.
get_db_in_use <- function(dbs_path, what=c("germline", "C-region"))
{
    stopifnot(isSingleNonWhiteString(dbs_path), dir.exists(dbs_path))
    what <- match.arg(what)
    db_name <- .DB_IN_USE_cache[[what]]
    if (is.null(db_name) || db_name == "")
        return("")  # no db is currently in use

    db_path <- file.path(dbs_path, db_name)
    if (!dir.exists(db_path)) {
        if (what == "germline") {
            fun <- "use_germline_db"
        } else {
            fun <- "use_c_region_db"
        }
        repair_with <- paste0("Try to repair with ", fun, "(\"<db_name>\").")
        see <- paste0("See '?", fun, "' for more information.")
        stop(wmsg("Anomaly: \"", db_name, "\" is not the name ",
                  "of a cached ", what, " db. ",
                  repair_with, " ", see))
    }
    db_path
}

set_db_in_use <- function(what=c("germline", "C-region"), db_name="",
                          verbose=FALSE)
{
    what <- match.arg(what)
    stopifnot(isSingleString(db_name), isTRUEorFALSE(verbose))
    print_ok <- FALSE
    if (verbose) {
        if (db_name == "") {
            old_db_name <- .DB_IN_USE_cache[[what]]
            if (!(is.null(old_db_name) || old_db_name == "")) {
                message("Cancelling the current ", what, " selection ... ",
                        appendLF=FALSE)
                print_ok <- TRUE
            }
        } else {
            message("Selecting ", what, " db ", db_name, " for use ",
                    "with igblastn() ... ", appendLF=FALSE)
            print_ok <- TRUE
        }
    }
    .DB_IN_USE_cache[[what]] <- db_name
    if (print_ok)
        message("ok")
    invisible(db_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### print_dbs_df()
###

### Used by print.germline_dbs_df() and print.c_region_dbs_df().
print_dbs_df <- function(dbs_df, dbs_path, what=c("germline", "C-region"))
{
    stopifnot(is.data.frame(dbs_df),
              isSingleNonWhiteString(dbs_path), dir.exists(dbs_path))
    what <- match.arg(what)
    dbs_df <- as.data.frame(dbs_df)
    db_names <- dbs_df[ , "db_name"]
    db_path <- get_db_in_use(dbs_path, what=what)
    if (db_path != "") {
        ## Mark db in use with an asterisk in extra white column.
        used <- character(length(db_names))
        used[db_names %in% basename(db_path)] <- "*"
        dbs_df <- cbind(dbs_df, ` `=used)
    }
    ## Left-justify the "db_name" column (1st column).
    col1 <- format(c("db_name", db_names), justify="left")
    dbs_df[ , "db_name"] <- col1[-1L]
    colnames(dbs_df)[[1L]] <- col1[[1L]]
    ## Do not print the row names.
    print(dbs_df, row.names=FALSE)
}

