
rows_with_same_key_are_identical <- function(df, key)
{
    stopifnot(is.data.frame(df), igblastr:::isSingleNonWhiteString(key))
    keys <- df[ , key]
    m <- match(keys, keys)
    df2 <- df[m, ]
    rownames(df2) <- NULL
    identical(df, df2)
}

