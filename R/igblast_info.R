### =========================================================================
### igblast_info() and related
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### has_igblast()
###

has_igblast <- function()
{
    igblast_root <- try(get_igblast_root(), silent=TRUE)
    !inherits(igblast_root, "try-error")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list_igblast_organisms()
###

list_igblast_organisms <- function(igblast_root=get_igblast_root())
{
    internal_data <- get_igblast_root_subdir(igblast_root, "internal_data")
    if (!dir.exists(internal_data))
        return(character(0))
    list.files(internal_data)
}

### Not exported!
normalize_igblast_organism <- function(organism)
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    all_organisms <- list_igblast_organisms()
    organism <- tolower(organism)
    if (!(organism %in% list_igblast_organisms())) {
        all_in_1string <- paste0("\"", all_organisms, "\"", collapse=", ")
        stop(wmsg("'organism' must be one of ", all_in_1string))
    }
    organism
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastn_version()
### makeblastdb_version()
### igblast_build()
### igblast_info()
###

.extract_version_string_from_raw_version <- function(raw_version)
{
    sub("[^:]*: *", "", raw_version[[1L]])
}

.extract_build_from_raw_version <- function(raw_version)
{
    sub("^ *Package: ", "", raw_version[[2L]])
}

.IGBLAST_VERSION_cache <- new.env(parent=emptyenv())

.get_igblast_exe_version <- function(cmdname, igblast_root, raw.version)
{
    if (!isTRUEorFALSE(raw.version))
        stop(wmsg("'raw.version' must be TRUE or FALSE"))
    ## 'cmdpath' is guaranteed to be an absolute path.
    cmdpath <- get_igblast_exe(cmdname, igblast_root=igblast_root, check=FALSE)
    raw_version <- .IGBLAST_VERSION_cache[[cmdpath]]
    if (is.null(raw_version)) {
        raw_version <- try(suppressWarnings(system2e(cmdpath, "-version",
                                                     stdout=TRUE, stderr=TRUE)),
                           silent=TRUE)
        if (!system_command_worked(raw_version))
            stop(wmsg("command '", cmdpath, " -version' failed"))
	.IGBLAST_VERSION_cache[[cmdpath]] <- raw_version
    }
    if (raw.version)
        return(raw_version)
    .extract_version_string_from_raw_version(raw_version)
}

igblastn_version <- function(igblast_root=get_igblast_root(),
                             raw.version=FALSE)
{
    .get_igblast_exe_version("igblastn", igblast_root, raw.version)
}

makeblastdb_version <- function(igblast_root=get_igblast_root(),
                                raw.version=FALSE)
{
    .get_igblast_exe_version("makeblastdb", igblast_root, raw.version)
}

### Returns a string of the form "igblast 1.22.0, build Oct 11 2023 16:06:20".
igblast_build <- function(igblast_root=get_igblast_root())
{
    raw_version <- igblastn_version(igblast_root, raw.version=TRUE)
    .extract_build_from_raw_version(raw_version)
}

print.igblast_info <- function(x, ...)
{
    x <- lapply(x, function(x) paste(x, collapse="; "))
    x <- paste0(names(x), ": ", as.character(x))
    cat(x, sep="\n")
}

igblast_info <- function(igblast_root=get_igblast_root())
{
    raw_version <- igblastn_version(igblast_root, raw.version=TRUE)
    igblast_bui <- .extract_build_from_raw_version(raw_version)
    igblastn_ver <- .extract_version_string_from_raw_version(raw_version)
    makeblastdb_ver <- makeblastdb_version(igblast_root)
    OS_arch <- get_OS_arch()
    all_organisms <- list_igblast_organisms(igblast_root)
    if (length(all_organisms) == 0L) {
        organisms <- "none!"
    } else {
        organisms <- paste0(all_organisms, collapse=", ")
    }

    ans <- list(
        igblast_root=igblast_root,
        igblast_build=igblast_bui,
        igblastn_version=igblastn_ver,
        makeblastdb_version=makeblastdb_ver,
        `OS/arch`=paste(OS_arch, collapse="/"),
        #igblastn_exe=igblastn_exe,
        #igblastn_version=raw_version[[1L]],
        #igblastp_exe=igblastp_exe,
        #igblastp_version=igblastp_ver[[1L]],
        organisms=organisms
    )
    class(ans) <- "igblast_info"
    ans
}

