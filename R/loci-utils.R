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
### normalize_loci()
###

### Returns loci in canonical order.
normalize_loci <- function(loci="auto", tcr.db=FALSE)
{
    if (!is.character(loci) || length(loci) == 0L)
        stop(wmsg("'loci' must be a non-empty character vector"))
    if (anyNA(loci))
        stop(wmsg("'loci' cannot contain NAs"))
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
    if (anyDuplicated(loci))
        stop(wmsg("'loci' cannot contain duplicates"))
    ## For now we don't allow mixing IG loci and TR loci. Note that it
    ## would not be too hard to support this if we wanted to but is there
    ## a reasonable use case for it?
    if (all(loci %in% IG_LOCI)) {
        loci_2_region_types <- .IG_LOCI_2_REGION_TYPES
        valid_loci <- IG_LOCI
    } else if (all(loci %in% TR_LOCI)) {
        loci_2_region_types <- .TR_LOCI_2_REGION_TYPES
        valid_loci <- TR_LOCI
    } else {
        IG_loci_in1string <- paste0("\"", IG_LOCI, "\"", collapse=", ")
        TR_loci_in1string <- paste0("\"", TR_LOCI, "\"", collapse=", ")
        stop(wmsg("'loci' must be a subset of 'c(", IG_loci_in1string, ")' ",
                  "or 'c(", TR_loci_in1string, ")'"))
    }
    ## Let's check that the user-selected subset of loci has regions of all
    ## 3 types (V/D/J). Note that all loci have at least V/J regions but
    ## some loci (e.g. IGK or TRA) don't have D regions.
    region_types <- unique(unlist(loci_2_region_types[loci], use.names=FALSE))
    missing_regions <- setdiff(VDJ_REGION_TYPES, region_types)
    if (length(missing_regions) != 0L) {
        ## Note that 'missing_regions' can only be the single string "D", but
        ## but we intentionally keep our code general.
        in1string1 <- paste0(missing_regions, collapse=", ")
        stop(wmsg("The selected subset of loci must have V/D/J regions. ",
                  "However, the current selection has no ", in1string1, " ",
                  "regions."))
    }
    if (!identical(tcr.db, FALSE))
        stop(wmsg("'tcr.db' is ignored when 'loci' is supplied"))
    valid_loci[valid_loci %in% loci]  # return loci in canonical order
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### stop_if_bad_loci()
###

stop_if_bad_loci <- function(loci)
{
    if (!is.character(loci) || length(loci) == 0L)
        stop(wmsg("'loci' must be a non-empty character vector"))
    if (anyNA(loci) || anyDuplicated(loci))
        stop(wmsg("'loci' cannot contain NAs or duplicates"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### map_loci_to_region_types()
###

map_loci_to_region_types <- function(loci)
{
    stop_if_bad_loci(loci)
    if (all(loci %in% IG_LOCI)) {
        loci2regiontypes <- .IG_LOCI_2_REGION_TYPES
    } else if (all(loci %in% TR_LOCI)) {
        loci2regiontypes <- .TR_LOCI_2_REGION_TYPES
    } else {
        what <- if (length(loci) != 1L) "set of " else ""
        in1string <- paste(loci, collapse=", ")
        stop(wmsg("invalid ", what, "loci: ", in1string))
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
    stop_if_bad_loci(selected_loci)
    loci_prefix <- unique(substr(selected_loci, 1L, 2L))
    stopifnot(length(loci_prefix) == 1L)
    all_loci <- .get_all_region_type_loci(region_type, loci_prefix)
    intersect(all_loci, selected_loci)
}

