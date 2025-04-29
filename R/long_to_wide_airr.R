.sorry_we_ve_moved <- function(fun)
    c("The ", fun, " function has moved to the bcrutils package available ",
      "at https://github.com/HyrienLab/bcrutils and will be removed ",
      "from igblastr soon. Please install bcrutils and use ",
      "bcrutils::", fun, " instead.")

### Assumes that the input data.frame has at least columns "sequence_id"
### and "locus". The pairing of the rows will be inferred from these
### two columns.
long_to_wide_airr <- function(df)
{
    .Deprecated(msg=wmsg(.sorry_we_ve_moved("long_to_wide_airr()")))
    stopifnot(is.data.frame(df),
              "sequence_id" %in% colnames(df),
              "locus" %in% colnames(df))

    if (is_tibble(df)) {
        if (!requireNamespace("dplyr", quietly=TRUE))
            stop("dplyr is required when the input is a tibble")
        cbind <- dplyr::bind_cols
    }

    ## Infer the chain from the locus.
    locus <- sub("^IG", "", df$locus)
    stopifnot(locus %in% c("H", "K", "L"))
    chain <- ifelse(locus == "H", "heavy", "light")

    ## Order rows first by sequence id, then by chain.
    oo <- order(df$sequence_id, chain)
    df <- df[oo, , drop=FALSE]

    ## Check that all rows are properly paired.
    rle <- S4Vectors::Rle(df$sequence_id)
    if (!all(S4Vectors::runLength(rle) == 2L))
        stop("not all rows are properly paired")

    ## Split the rows in two groups: 'heavy_df' and 'light_df'.
    ans_nrow <- nrow(df) %/% 2L
    light_rowidx <- seq_len(ans_nrow) * 2L
    heavy_rowidx <- light_rowidx - 1L
    heavy_df <- df[heavy_rowidx, , drop=FALSE]  # extract odd rows
    light_df <- df[light_rowidx, , drop=FALSE]  # extract even rows
    rownames(heavy_df) <- rownames(light_df) <- NULL

    ## Sanity checks.
    heavy_locus <- sub("^IG", "", heavy_df$locus)
    if (!all(heavy_locus == "H"))
        stop("some pairs don't have a heavy member")
    light_locus <- sub("^IG", "", light_df$locus)
    if (!all(light_locus %in% c("K", "L")))
        stop("some pairs don't have a light member")
    stopifnot(all(heavy_df$sequence_id == light_df$sequence_id))

    ## Remove "sequence_id" column from 'heavy_df' and 'light_df'.
    sequence_id <- heavy_df$sequence_id
    m <- match("sequence_id", colnames(heavy_df))
    heavy_df <- heavy_df[ , -m, drop=FALSE]
    m <- match("sequence_id", colnames(light_df))
    light_df <- light_df[ , -m, drop=FALSE]

    ## Make and return the wide data.frame.
    colnames(heavy_df) <- paste0(colnames(heavy_df), "_heavy")
    colnames(light_df) <- paste0(colnames(light_df), "_light")
    cbind(sequence_id=sequence_id, heavy_df, light_df)
}

wide_to_long_airr <- function(df)
{
    .Deprecated(msg=wmsg(.sorry_we_ve_moved("wide_to_long_airr()")))
    stopifnot(is.data.frame(df))

    if (is_tibble(df)) {
        if (!requireNamespace("dplyr", quietly=TRUE))
            stop("dplyr is required when the input is a tibble")
        cbind <- dplyr::bind_cols
    }

    ## Identify _heavy and _light columns.
    is_heavy_col <- has_suffix(colnames(df), "_heavy")
    heavy_col_count <- sum(is_heavy_col)
    if (heavy_col_count == 0L)
        stop("input has no _heavy columns")
    is_light_col <- has_suffix(colnames(df), "_light")
    if (sum(is_light_col) != heavy_col_count)
        stop("set of _heavy columns doesn't match set of _light columns")

    ## Split the columns in two groups: 'heavy_df' and 'light_df'.
    heavy_df <- df[ , is_heavy_col, drop=FALSE]
    colnames(heavy_df) <- sub("_heavy$", "", colnames(heavy_df))
    light_df <- df[ , is_light_col, drop=FALSE]
    colnames(light_df) <- sub("_light$", "", colnames(light_df))
    if (!setequal(colnames(heavy_df), colnames(light_df)))
        stop("set of _heavy columns doesn't match set of _light columns")

    ## Make the long data.frame.
    light_df <- light_df[ , colnames(heavy_df), drop=FALSE]
    other_cols <- df[ , !(is_heavy_col | is_light_col), drop=FALSE]
    heavy_df <- cbind(other_cols, heavy_df)
    light_df <- cbind(other_cols, light_df)
    long_df <- rbind(heavy_df, light_df)

    ## Reorder rows of long data.frame so pairs are made of adjacent rows.
    oo <- integer(nrow(long_df))
    oo[c(TRUE, FALSE)] <- seq_len(nrow(df))
    oo[c(FALSE, TRUE)] <- seq_len(nrow(df)) + nrow(df)
    long_df <- long_df[oo, , drop=FALSE]
    rownames(long_df) <- NULL
    long_df
}

