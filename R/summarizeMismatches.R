### =========================================================================
### Tabulate mismatches and indels between query and germline sequences
### -------------------------------------------------------------------------

### The functionality implemented in this file relies on the
### GenomicAlignments package for the low-level CIGAR utilities and
### the sequenceLayer() function.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_ROI_ranges()
###

### TODO: Maybe add the C region to the list of Regions Of Interest. Should
### be on user request, not by default. Also, we can only tabulate mismatches
### and indels that fall in that region if a C-region db was specified.
### The C region corresponds to the (c_sequence_start, c_sequence_end) interval
### and is located after the (fwr4_start, fwr4_end) interval. In principle,
### it should be adjacent to the FWR4 region. However, one small gotcha is
### that, for some queries, fwr4_end is equal to c_sequence_start in the
### output produced by igblastn. This means that there's sometimes a
### 1-nucleotide overlap between the FWR4 region and the C region!
### So in order to keep our Regions Of Interest adjacent and non-overlapping,
### let's use the (fwr4_start, c_sequence_end) interval for this region. And
### let's report the counts of mismatches and/or indels that fall in this
### interval in the "c" column of the returned matrix.

.FWRCDR_JUNC_NAMES <- paste0(c("", FWRCDR_NAMES), ".", c(FWRCDR_NAMES, ""))

.insert_junc_names_between_fwrcdr_names <- function()
{
    ans <- character(length(.FWRCDR_JUNC_NAMES) + length(FWRCDR_NAMES))
    ans[c(TRUE, FALSE)] <- .FWRCDR_JUNC_NAMES
    ans[c(FALSE, TRUE)] <- FWRCDR_NAMES
    ans
}

### TODO: Don't use $ to access a column because that works even if the
### column does not exist (in which case it returns NULL). Instead, we
### want to fail with an informative error message if the column does
### not exist.
.extract_ROI_ranges <- function(AIRR_df, with.junctions=FALSE)
{
    stopifnot(is.data.frame(AIRR_df), isTRUEorFALSE(with.junctions))
    all_ranges <- list(
        fwr1=IRanges(AIRR_df$fwr1_start, AIRR_df$fwr1_end),
        cdr1=IRanges(AIRR_df$cdr1_start, AIRR_df$cdr1_end),
        fwr2=IRanges(AIRR_df$fwr2_start, AIRR_df$fwr2_end),
        cdr2=IRanges(AIRR_df$cdr2_start, AIRR_df$cdr2_end),
        fwr3=IRanges(AIRR_df$fwr3_start, AIRR_df$fwr3_end),
        cdr3=IRanges(AIRR_df$cdr3_start, AIRR_df$cdr3_end),
        fwr4=IRanges(AIRR_df$fwr4_start, AIRR_df$fwr4_end)
    )
    if (with.junctions) {
        ## Note that junction ranges will be used in the context of finding
        ## overlaps with deletion ranges which are 0-width ranges. However,
        ## a 0-width range only has a hit with another range if it falls
        ## strictly within the latter. This is why our junction ranges span
        ## two adjacent base positions.
        junction_ranges <- list(
            .fwr1    =IRanges(AIRR_df$fwr1_start-1L, AIRR_df$fwr1_start),
            fwr1.cdr1=IRanges(AIRR_df$fwr1_end, AIRR_df$cdr1_start),
            cdr1.fwr2=IRanges(AIRR_df$cdr1_end, AIRR_df$fwr2_start),
            fwr2.cdr2=IRanges(AIRR_df$fwr2_end, AIRR_df$cdr2_start),
            cdr2.fwr3=IRanges(AIRR_df$cdr2_end, AIRR_df$fwr3_start),
            fwr3.cdr3=IRanges(AIRR_df$fwr3_end, AIRR_df$cdr3_start),
            cdr3.fwr4=IRanges(AIRR_df$cdr3_end, AIRR_df$fwr4_start),
            fwr4.    =IRanges(AIRR_df$fwr4_end, AIRR_df$fwr4_end+1L)
        )
        all_ranges <- c(all_ranges, junction_ranges)
        stopifnot(identical(names(all_ranges),
                            c(FWRCDR_NAMES, .FWRCDR_JUNC_NAMES)))
        all_ranges_names <- .insert_junc_names_between_fwrcdr_names()
        all_ranges <- all_ranges[all_ranges_names]
        stopifnot(identical(names(all_ranges), all_ranges_names))
    }
    ## Regroup ranges by query sequence.
    all_ranges <- IRangesList(all_ranges)
    f <- rep.int(seq_len(nrow(AIRR_df)), length(all_ranges))
    ans <- unname(split(unlist(all_ranges), f))
    ## Sanity checks.
    stopifnot(
        identical(lengths(ans), rep.int(length(all_ranges), nrow(AIRR_df))),
        identical(names(unlist(ans, use.names=FALSE)),
                  rep.int(names(all_ranges), nrow(AIRR_df)))
    )
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .count_hits_per_ROI()
###

.check_ROI_ranges <- function(ROI_ranges, N, ROI_names)
{
    stopifnot(is(ROI_ranges, "CompressedIRangesList"),
              identical(lengths(ROI_ranges), rep.int(length(ROI_names), N)),
              identical(names(unlist(ROI_ranges, use.names=FALSE)),
                        rep.int(ROI_names, N)))
}

### Returns 'length(x)' x 'length(colnames)' matrix.
.list2matrix <- function(x, colnames)
{
    stopifnot(is.list(x), is.character(colnames),
              all(lengths(x) == length(colnames)))
    matrix(unlist(x, use.names=FALSE), ncol=length(colnames),
           byrow=TRUE, dimnames=list(NULL, colnames))
}

### 'query' and 'ROI_ranges' must be two **parallel** IntegerRangesList
### derivatives of length N.
### The following strong assumption is made: for all 1 <= i <= N, every
### range in 'query[[i]]' is expected to overlap with **exactly** one range
### in 'ROI_ranges[[i]]'. The function will return an error if this is not
### the case.
### Returns an integer matrix with N rows (one row per list element in 'query'
### or in 'ROI_ranges') and one column per Region Of Interest.
.count_hits_per_ROI <- function(query, ROI_ranges, ROI_names,
                                use.query.weights=FALSE)
{
    stopifnot(is(query, "CompressedIntegerRangesList"),
              is.character(ROI_names),
              isTRUEorFALSE(use.query.weights))
    N <- length(query)
    .check_ROI_ranges(ROI_ranges, N, ROI_names)
    if (use.query.weights)
        stopifnot("weight" %in% colnames(mcols(unlist(query, use.names=FALSE))))

    all_hits <- findOverlaps(query, ROI_ranges)

    ## Because of the "exactly one hit per range in 'query'" assumption,
    ## 'all_hits' is expected to have the same shape as 'query'.
    stopifnot(identical(lengths(all_hits), lengths(query)))

    ## Extract nb of hits per Region Of Interest.
    ## Note that we use an lapply() loop for this at the moment, which is
    ## not very efficient. TODO: Can we avoid the lapply() loop?
    all_counts <- lapply(seq_along(all_hits),
        function(i) {
            hits <- all_hits[[i]]
            q <- query[[i]]
            stopifnot(identical(queryHits(hits), seq_along(q)))
            if (use.query.weights) {
                qweights <- mcols(q)$weight
                hit_counts <- S4Vectors:::tabulate2(subjectHits(hits),
                                                    nbins=length(ROI_names),
                                                    weight=qweights)
            } else {
                hit_counts <- tabulate(subjectHits(hits),
                                       nbins=length(ROI_names))
            }
            setNames(hit_counts, ROI_names)
        })

    ## Turn list of integer vectors into a matrix with one row per list
    ## element in 'all_counts'.
    .list2matrix(all_counts, ROI_names)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .compute_mm_pos()
### .compute_ins_pos()
### .compute_del_ranges()
###

.extract_cigar <- function(AIRR_df, region_type)
{
    stopifnot(is.data.frame(AIRR_df),
              isSingleNonWhiteString(region_type), nchar(region_type) == 1L)
    cigar_colname <- paste0(tolower(region_type), "_cigar")
    cigar <- AIRR_df[[cigar_colname]]
    if (is.null(cigar))
        stop(wmsg("'AIRR_df' has no \"", cigar_colname, "\" column"))
    cigar[is.na(cigar)] <- ""
    cigar
}

.extract_call <- function(AIRR_df, region_type)
{
    stopifnot(is.data.frame(AIRR_df),
              isSingleNonWhiteString(region_type), nchar(region_type) == 1L)
    call_colname <- paste0(tolower(region_type), "_call")
    call <- AIRR_df[[call_colname]]
    if (is.null(call))
        stop(wmsg("'AIRR_df' has no \"", call_colname, "\" column"))
    call
}

.extract_rseqs <- function(germline_db, call)
{
    stopifnot(is(germline_db, "DNAStringSet"), !is.null(names(germline_db)))
    NOT_A_VALID_ALLELE_NAME <- "_not_a_valid_allele_name_"
    call[is.na(call)] <- NOT_A_VALID_ALLELE_NAME
    c(germline_db, DNAStringSet(setNames("", NOT_A_VALID_ALLELE_NAME)))[call]
}

### Returns the positions of the nucleotide mismatches in an IPosList object
### that has one list element per row in 'AIRR_df'.
### Only the positions that fall in the (fwr1_start, fwr4_end) interval
### are returned.
### TODO: Look into obtaining 'rseqs' (see below) directly from
### the 'germline_alignment' column (or one of the '[vdj]_germline_alignment'
### columns). If that works, then we can get rid of the 'germline_db' argument
### which is inconvenient to use because it requires that the user still has
### access to the germline or C-region db that was used by igblastn() to
### produce 'AIRR_df'. In other words, .compute_mm_pos() would be able
### to get everything it needs from 'AIRR_df' and nothing else, like
### .compute_ins_pos() and .compute_del_pos().
.compute_mm_pos <- function(AIRR_df, region_type, germline_db)
{
    cigar <- .extract_cigar(AIRR_df, region_type)
    call <- .extract_call(AIRR_df, region_type)
    nona_idx <- which(!is.na(call))
    rseqs <- .extract_rseqs(germline_db, call)
    region_type <- tolower(region_type)
    sequence_start_colname <- paste0(region_type, "_sequence_start")
    sequence_end_colname <- paste0(region_type, "_sequence_end")
    sequence_start <- AIRR_df[[sequence_start_colname]]
    sequence_end <- AIRR_df[[sequence_end_colname]]

    projected_rseqs <-
        GenomicAlignments::sequenceLayer(rseqs, cigar,
                                         from="reference", to="query")
    qseqs <- DNAStringSet(AIRR_df$sequence)
    stopifnot(width(projected_rseqs)[nona_idx] == width(qseqs)[nona_idx])

    trimmed_qseqs <- subseq(qseqs, sequence_start, sequence_end)
    trimmed_rseqs <- subseq(projected_rseqs, sequence_start, sequence_end)

    ## Sanity checks.
    expected_trimmed_lens <-
        GenomicAlignments::cigarWidthAlongQuerySpace(cigar,
                                                     after.soft.clipping=TRUE)
    expected_trimmed_lens <- expected_trimmed_lens[nona_idx]
    stopifnot(
        identical(width(trimmed_qseqs)[nona_idx], expected_trimmed_lens),
        identical(width(trimmed_rseqs)[nona_idx], expected_trimmed_lens)
    )

    ## Compute 'mm_pos' as an IPosList object.
    ## Unfortunately we use an lapply() loop for this, which is not very
    ## efficient. TODO: Can we avoid the lapply() loop?
    mm_pos <- lapply(seq_len(nrow(AIRR_df)),
        function(i) {
            if (is.na(call[[i]]))
                return(integer(0))
            qbytes <- as.raw(trimmed_qseqs[[i]])
            rbytes <- as.raw(trimmed_rseqs[[i]])
            is_mm <- qbytes != rbytes & rbytes != as.raw(DNAString("-"))
            mm_pos <- which(is_mm) + sequence_start[[i]] - 1L
            ## Keep only positions that fall in the (fwr1_start, fwr4_end)
            ## interval.
            keep_idx <- mm_pos >= AIRR_df$fwr1_start[[i]] &
                        mm_pos <= AIRR_df$fwr4_end[[i]]
            mm_pos[keep_idx]
        })

    relist(IPos(unlist(mm_pos, use.names=FALSE)), mm_pos)
}

### Returns the positions of the inserted nucleotides in an IPosList object
### that has one list element per row in 'AIRR_df'.
### Only the positions that fall in the (fwr1_start, fwr4_end) interval
### are returned.
.compute_ins_pos <- function(AIRR_df, region_type)
{
    cigar <- .extract_cigar(AIRR_df, region_type)

    ## Compute 'ins_ranges'.
    ins_ranges <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar, ops="I")

    ## Restrict ranges in 'ins_ranges' to the (fwr1_start, fwr4_end)
    ## interval.
    ins_ranges <- restrict(ins_ranges, start=AIRR_df$fwr1_start,
                                       end=AIRR_df$fwr4_end)

    ## Turn IRangesList object 'ins_ranges' into IPosList object.
    unlisted <- IPos(unlist(ins_ranges, use.names=FALSE))
    partitioning <- PartitioningByWidth(sum(width(ins_ranges)))
    relist(unlisted, partitioning)
}

### Returns the ranges of the deletions (0-width ranges) in an IRangesList
### object that has one list element per row in 'AIRR_df'.
### Only the ranges that are within (or adjacent to) the (fwr1_start, fwr4_end)
### interval are returned.
.compute_del_ranges <- function(AIRR_df, region_type)
{
    cigar <- .extract_cigar(AIRR_df, region_type)

    ## Compute 'del_ranges' and 'ref_del_ranges'.
    ## The ranges in 'del_ranges' are the deletion ranges w.r.t. the
    ## query sequences, so they're all expected to be 0-width ranges.
    ## The ranges in 'ref_del_ranges' are the deletion ranges w.r.t. the
    ## germline sequences. These should all have a width >= 1.
    ## Note that 'del_ranges' and 'ref_del_ranges' are both
    ## IRangesList objects that are expected to have the same shape, that is:
    ## - They should have the same length N, where N is the number of
    ##   query sequences or 'nrow(AIRR_df)'. In other words, they're both
    ##   expected to have one list element per row in 'AIRR_df'.
    ## - The number of ranges in 'del_ranges[[i]]'
    ##   and 'ref_del_ranges[[i]]' must be the same for all 1 <= i <= N.
    del_ranges <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar, ops="D")
    stopifnot(all(width(unlist(del_ranges, use.names=FALSE)) == 0L))
    ref_del_ranges <-
        GenomicAlignments::cigarRangesAlongReferenceSpace(cigar, ops="D")
    stopifnot(identical(lengths(del_ranges), lengths(ref_del_ranges)))

    ## Compute the "weights" of the deletions.
    ## We call them "weights" but they're actually the lengths of the
    ## deletions i.e. the number of nucleotides that is deleted for each
    ## deletion.
    del_weights <- width(ref_del_ranges)  # IntegerList
    ## Not proud of this (there's a cleaner way to do this).
    mcols(del_ranges@unlistData)$weight <- del_weights@unlistData

    ## Keep only ranges in 'del_ranges' that are within (or adjacent to)
    ## the (fwr1_start, fwr4_end) interval.
    ## Because restrict() will possibly drop some 0-width ranges,
    ## we will no longer be able to assume that 'del_ranges'
    ## and 'ref_del_ranges' (or 'del_weights') have the same shape.
    ## So it's important that we put the weights on 'del_ranges'
    ## **before** we call restrict().
    restrict(del_ranges, start=AIRR_df$fwr1_start, end=AIRR_df$fwr4_end)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .tabulate_mismatches_for_region_type()
### .tabulate_insertions_for_region_type()
### .tabulate_deletions_for_region_type()
###

### Infer nb of nucleotide mismatches in each FWR/CDR region from
### a given CIGAR column (one of "[vdjc]_cigar").
.tabulate_mismatches_for_region_type <-
    function(AIRR_df, region_type, ROI_ranges, germline_db)
{
    mm_pos <- .compute_mm_pos(AIRR_df, region_type, germline_db)
    .count_hits_per_ROI(mm_pos, ROI_ranges, FWRCDR_NAMES)
}

### Infer nb of single nucleotide insertions in each FWR/CDR region from
### a given CIGAR column (one of "[vdjc]_cigar").
.tabulate_insertions_for_region_type <-
    function(AIRR_df, region_type, ROI_ranges)
{
    ins_pos <- .compute_ins_pos(AIRR_df, region_type)
    .count_hits_per_ROI(ins_pos, ROI_ranges, FWRCDR_NAMES)
}

### Infer nb of single nucleotide deletions in each FWR/CDR region and
### junction from a given CIGAR column (one of "[vdjc]_cigar").
.tabulate_deletions_for_region_type <-
    function(AIRR_df, region_type, ROI_ranges)
{
    del_ranges <- .compute_del_ranges(AIRR_df, region_type)
    ROI_names <- .insert_junc_names_between_fwrcdr_names()
    .count_hits_per_ROI(del_ranges, ROI_ranges, ROI_names,
                        use.query.weights=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### tabulate_mismatches()
### tabulate_insertions()
### tabulate_deletions()
###

.collapse_counts <- function(all_counts)
{
    a_dimnames <- list(NULL, colnames(all_counts[[1L]]), names(all_counts))
    a <- array(unlist(all_counts, use.names=FALSE),
               dim=c(dim(all_counts[[1L]]), length(all_counts)),
               dimnames=a_dimnames)
    rowSums(a, dims=2L)
}

tabulate_mismatches <- function(AIRR_df, germline_db_name, c_region_db_name="")
{
    ## TODO: Try to load GenomicAlignments namespace and fail graciously
    ## if the package is not installed.
    ROI_ranges <- .extract_ROI_ranges(AIRR_df)
    region_types <- VDJ_REGION_TYPES
    all_counts <- lapply(setNames(region_types, region_types),
        function(region_type) {
            germline_db <- load_germline_sequences(germline_db_name,
                                                   region_type)
            .tabulate_mismatches_for_region_type(AIRR_df, region_type,
                                                 ROI_ranges, germline_db)
        })
    if ("c_cigar" %in% colnames(AIRR_df)) {
        if (c_region_db_name == "")
            stop(wmsg("please specify the name of the C-region db ",
                      "that was used to produce 'AIRR_df'"))
        c_region_db <- load_c_region_sequences(c_region_db_name)
        c_counts <- .tabulate_mismatches_for_region_type(AIRR_df, "C",
                                                 ROI_ranges, c_region_db)
        all_counts <- c(all_counts, list(C=c_counts))
    }
    .collapse_counts(all_counts)
}

tabulate_insertions <- function(AIRR_df)
{
    ## TODO: Try to load GenomicAlignments namespace and fail graciously
    ## if the package is not installed.
    ROI_ranges <- .extract_ROI_ranges(AIRR_df)
    region_types <- VDJC_REGION_TYPES
    all_counts <- lapply(setNames(region_types, region_types),
        function(region_type) {
            .tabulate_insertions_for_region_type(AIRR_df, region_type,
                                                 ROI_ranges)
        })
    .collapse_counts(all_counts)
}

tabulate_deletions <- function(AIRR_df)
{
    ## TODO: Try to load GenomicAlignments namespace and fail graciously
    ## if the package is not installed.
    ROI_ranges <- .extract_ROI_ranges(AIRR_df, with.junctions=TRUE)
    region_types <- VDJC_REGION_TYPES
    all_counts <- lapply(setNames(region_types, region_types),
        function(region_type) {
            .tabulate_deletions_for_region_type(AIRR_df, region_type,
                                                ROI_ranges)
        })
    .collapse_counts(all_counts)
}

