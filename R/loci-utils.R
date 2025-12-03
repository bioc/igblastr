### =========================================================================
### Low-level loci utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some fundamental global constants
###

VDJ_REGION_TYPES <- c("V", "D", "J")
.VJ_REGION_TYPES <- VDJ_REGION_TYPES[-2L]

### Group names are formed by concatenating a locus name (e.g. IGH or TRB)
### and a region type (e.g. V).

.IG_LOCI_2_REGION_TYPES <- list(IGH=VDJ_REGION_TYPES,
                                IGK=.VJ_REGION_TYPES,
                                IGL=.VJ_REGION_TYPES)

.TR_LOCI_2_REGION_TYPES <- list(TRA=.VJ_REGION_TYPES,
                                TRB=VDJ_REGION_TYPES,
                                TRG=.VJ_REGION_TYPES,
                                TRD=VDJ_REGION_TYPES)

IG_LOCI <- names(.IG_LOCI_2_REGION_TYPES)
TR_LOCI <- names(.TR_LOCI_2_REGION_TYPES)

.revmap <- function(loci2regiontypes)
{
    loci <- rep.int(names(loci2regiontypes), lengths(loci2regiontypes))
    f <- factor(unlist(loci2regiontypes, use.names=FALSE),
                levels=VDJ_REGION_TYPES)
    split(loci, f)
}

IG_REGION_TYPES_2_LOCI <- .revmap(.IG_LOCI_2_REGION_TYPES)
TR_REGION_TYPES_2_LOCI <- .revmap(.TR_LOCI_2_REGION_TYPES)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### stop_if_malformed_loci_vector()
###

stop_if_malformed_loci_vector <- function(loci)
{
    if (!is.character(loci) || length(loci) == 0L)
        stop(wmsg("'loci' must be a non-empty character vector"))
    if (anyNA(loci) || anyDuplicated(loci))
        stop(wmsg("'loci' cannot contain NAs or duplicates"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_loci_prefix()
###

### Also checks that the supplied loci are valid loci.
### Returns "IG" or "TR".
extract_loci_prefix <- function(loci)
{
    stop_if_malformed_loci_vector(loci)
    if (all(loci %in% IG_LOCI))
        return("IG")
    if (all(loci %in% TR_LOCI))
        return("TR")
    IG_loci_in1string <- paste0("\"", IG_LOCI, "\"", collapse=", ")
    TR_loci_in1string <- paste0("\"", TR_LOCI, "\"", collapse=", ")
    stop(wmsg("'loci' must be a subset of 'c(", IG_loci_in1string, ")' ",
              "or a subset of 'c(", TR_loci_in1string, ")'"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normalize_user_supplied_loci()
###

### Raise an error if the supplied set of loci as a whole does not have
### regions of all 3 types (V/D/J). Note that all loci have at least V/J
### regions but some loci (e.g. IGK or TRA) don't have D regions.
.check_loci_for_missing_regions <- function(loci)
{
    loci_prefix <- extract_loci_prefix(loci)
    loci_2_region_types <- switch(loci_prefix,
                                  IG=.IG_LOCI_2_REGION_TYPES,
                                  TR=.TR_LOCI_2_REGION_TYPES,
                                  stop("unknown loci prefix: ", loci_prefix))
    region_types <- unique(unlist(loci_2_region_types[loci], use.names=FALSE))
    missing_regions <- setdiff(VDJ_REGION_TYPES, region_types)
    if (length(missing_regions) != 0L) {
        ## Note that 'missing_regions' can only be the single string "D",
        ## but we intentionally keep our code general.
        in1string1 <- paste0(missing_regions, collapse=", ")
        stop(wmsg("The selected subset of loci must have V/D/J regions. ",
                  "However, the current selection has no ", in1string1, " ",
                  "regions."))
    }
}

### If 'loci' is "auto" then 'stop.if.missing.regions' is ignored.
### Otherwise 'tcr.db' is ignored.
### Returns loci in canonical order.
normalize_user_supplied_loci <- function(loci="auto", tcr.db=FALSE,
                                         stop.if.missing.regions=FALSE)
{
    stop_if_malformed_loci_vector(loci)
    if (length(loci) == 1L) {
        if (loci == "auto") {
            if (!isTRUEorFALSE(tcr.db))
                stop(wmsg("'tcr.db' must be TRUE or FALSE"))
            return(if (tcr.db) TR_LOCI else IG_LOCI)
        }
        if (is_white_str(loci))
            stop(wmsg("'loci' cannot be a white string"))
        loci <- trimws2(strsplit(loci, "+", fixed=TRUE)[[1L]])
    }
    if (!isTRUEorFALSE(stop.if.missing.regions))
        stop(wmsg("'stop.if.missing.regions' must be TRUE or FALSE"))
    if (stop.if.missing.regions)
        .check_loci_for_missing_regions(loci)
    loci_prefix <- extract_loci_prefix(loci)
    if (!identical(tcr.db, FALSE))
        stop(wmsg("'tcr.db' should not be used when 'loci' is supplied"))
    valid_loci <- if (loci_prefix == "IG") IG_LOCI else TR_LOCI
    valid_loci[valid_loci %in% loci]  # return loci in canonical order
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### map_loci_to_region_types()
###

map_loci_to_region_types <- function(loci)
{
    loci_prefix <- extract_loci_prefix(loci)
    if (loci_prefix == "IG") {
        loci2regiontypes <- .IG_LOCI_2_REGION_TYPES
    } else {
        loci2regiontypes <- .TR_LOCI_2_REGION_TYPES
    }
    loci2regiontypes[loci]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_region_type_loci()
###

.get_all_region_type_loci <- function(region_type, loci_prefix=c("IG", "TR"))
{
    stopifnot(isSingleNonWhiteString(region_type))
    region_type <- match.arg(region_type, VDJ_REGION_TYPES)
    stopifnot(isSingleNonWhiteString(loci_prefix))
    loci_prefix <- match.arg(loci_prefix)
    if (loci_prefix == "IG") {
        region_types_2_loci <- IG_REGION_TYPES_2_LOCI
    } else {
        region_types_2_loci <- TR_REGION_TYPES_2_LOCI
    }
    region_types_2_loci[[region_type]]
}

get_region_type_loci <- function(region_type, selected_loci)
{
    loci_prefix <- extract_loci_prefix(selected_loci)
    all_loci <- .get_all_region_type_loci(region_type, loci_prefix)
    intersect(all_loci, selected_loci)
}

