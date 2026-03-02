### =========================================================================
### Low-level cache utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###
### The igblastr cache is divided into various compartments that are
### completely independent one from each other. These compartments are:
###
###   compartment     path
###   --------------  ------------------------------------
###   IGBLAST_ROOTS   <igblastr-cache>/igblast_roots
###   LIVE_IGDATA     <igblastr-cache>/live_igdata
###   GERMLINE_DBS    <igblastr-cache>/germline_dbs
###   C_REGION_DBS    <igblastr-cache>/c_region_dbs
###   IMGT_STORE      <igblastr-cache>/store/IMGT-releases
###   OGRDB_STORE     <igblastr-cache>/store/OGRDB
###
### Always use the igblastr_cache() function implemented in this file to
### obtain these paths.

IGBLAST_ROOTS <- "IGBLAST_ROOTS"
LIVE_IGDATA   <- "LIVE_IGDATA"
GERMLINE_DBS  <- "GERMLINE_DBS"
C_REGION_DBS  <- "C_REGION_DBS"
IMGT_STORE    <- "IMGT_STORE"
OGRDB_STORE   <- "OGRDB_STORE"

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastr_cache()
###

.IGBLASTR_CACHES <- c(IGBLAST_ROOTS, LIVE_IGDATA,
                      GERMLINE_DBS, C_REGION_DBS,
                      IMGT_STORE, OGRDB_STORE)

### Returns **absolute** path to igblastr's persistent cache.
### Note that the returned path is guaranteed to be a directory but is not
### guaranteed to be a **writable** directory. If it's not, then bad things
### will happen downstream.
.get_cache_path <- function()
{
    cache_path <- getOption("igblastr_cache")
    if (is.null(cache_path)) {
        cache_path <- R_user_dir("igblastr", "cache")
    } else {
        if (!isSingleNonWhiteString(cache_path))
            stop(wmsg("global option \"igblastr_cache\" must ",
                      "be set to a single (non-empty) string"))
    }
    if (!dir.exists(cache_path)) {
        ok <- dir.create(cache_path, recursive=TRUE)
        if (!ok) {
            msg1 <- c("Failed to create igblastr's persistent cache ",
                      "at: ", cache_path)
            msg2 <- c("Please use 'options()' to set global option ",
                      "\"igblastr_cache\" to the path of a directory ",
                      "that you are able to create.")
            stop(wmsg(msg1), "\n  ", wmsg(msg2))
        }
    }
    file_path_as_absolute(cache_path)
}

### Returns **absolute** path to the cache compartment specified via 'which'.
igblastr_cache <- function(which=NULL)
{
    cache_path <- .get_cache_path()
    if (is.null(which))
        return(cache_path)
    stopifnot(isSingleNonWhiteString(which))
    switch(which,
        IGBLAST_ROOTS=file.path(cache_path, "igblast_roots"),
        LIVE_IGDATA  =file.path(cache_path, "live_igdata"),
        GERMLINE_DBS =file.path(cache_path, "germline_dbs"),
        C_REGION_DBS =file.path(cache_path, "c_region_dbs"),
        IMGT_STORE   =file.path(cache_path, "store", "IMGT-releases"),
        OGRDB_STORE  =file.path(cache_path, "store", "OGRDB"),
        stop(wmsg("'which' must be one of ",
                  paste0("\"", .IGBLASTR_CACHES, "\"", collapse=", ")))
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_internal_igblast_roots()
###

get_internal_igblast_roots <- function() igblastr_cache(IGBLAST_ROOTS)

