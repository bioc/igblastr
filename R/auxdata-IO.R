### =========================================================================
### Read/write ".aux" files
### -------------------------------------------------------------------------
###


### Not exported!
### Not the true colnames used in IgBLAST ".aux" files: ours are shorter,
### all lowercase, and contain underscores instead of spaces.
AUXDATA_COL2CLASS <- c(
    allele_name="character",
    coding_frame_start="integer",
    chain_type="character",
    cdr3_end="integer",
    extra_bps="integer"
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### check_auxdata_col2class()
###

### Not exported!
check_auxdata_col2class <- function(auxdata, what="'auxdata'")
{
    if (!is.data.frame(auxdata))
        stop(wmsg(what, " must be a data.frame"))
    expected_colnames <- names(AUXDATA_COL2CLASS)
    if (!identical(colnames(auxdata), expected_colnames)) {
        in1string <- paste(expected_colnames, collapse=", ")
        stop(wmsg(what, " must have the following columns ",
                  "(in this order): ", in1string))
    }
    col2class <- vapply(auxdata, function(x) class(x)[[1L]], character(1))
    if (!identical(col2class, AUXDATA_COL2CLASS)) {
        in1string <- paste0("    ", names(AUXDATA_COL2CLASS), " -> ",
                            AUXDATA_COL2CLASS, collapse="\n")
        stop(wmsg(what, " must have the following column types:"), "\n",
             in1string)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_auxdata()
### write_auxdata()
###

read_auxdata <- function(filepath)
{
    read_broken_table(filepath, AUXDATA_COL2CLASS)
}

write_auxdata <- function(auxdata, file="")
{
    check_auxdata_col2class(auxdata)
    header <- paste0("#", paste(colnames(auxdata), collapse=", "))
    cat(header, "\n", sep="", file=file)
    write.table(auxdata, file, append=TRUE, quote=FALSE,
                sep="\t", row.names=FALSE, col.names=FALSE)
}

