### =========================================================================
### allele2gene()
### -------------------------------------------------------------------------
###

allele2gene <- function(allele_names)
{
    if (!is.character(allele_names))
        stop(wmsg("'allele_names' must be a character vector"))
    sub("\\*.*$", "", allele_names)
}

