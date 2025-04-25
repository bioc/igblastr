
### Assumes that the input data.frame has at least columns "sequence_id"
### and "locus". The pairing of the rows will be inferred from these
### two columns.
long_to_wide_airr <- function(df)
{
    stopifnot(is.data.frame(df),
              "sequence_id" %in% colnames(df),
              "locus" %in% colnames(df))

    ## Infer the chain from the locus.
    locus <- sub("^IG", "", df$locus)
    stopifnot(locus %in% c("H", "K", "L"))
    chain <- ifelse(locus == "H", "heavy", "light")

    ## Order rows first by sequence id, then by chain.
    oo <- order(df$sequence_id, chain)
    df <- df[oo, ]

    ## Check that all rows are properly paired.
    rle <- S4Vectors::Rle(df$sequence_id)
    if (!all(S4Vectors::runLength(rle) == 2L))
        stop("not all rows are properly paired")

    heavy_df <- df[c(TRUE, FALSE), ]  # extract odd rows
    light_df <- df[c(FALSE, TRUE), ]  # extract even rows

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
    heavy_df <- heavy_df[ , -match("sequence_id", colnames(heavy_df))]
    light_df <- light_df[ , -match("sequence_id", colnames(light_df))]

    ## Make and return the wide data.frame.
    colnames(heavy_df) <- paste0(colnames(heavy_df), "_heavy")
    colnames(light_df) <- paste0(colnames(light_df), "_light")
    cbind(sequence_id=sequence_id, heavy_df, light_df)
}

