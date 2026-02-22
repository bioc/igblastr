### Typical usage:
###
###     expect_error2(some op , "not supported")
###
### This won't break like expect_error(some op, "not supported") does if
### the error message happens to contain other white spaces like \n's or
### \t's instead of a single space between "not" and "supported". These
### can be introduced in an unpredictable way by 'stop(wmsg(...))'.
expect_error2 <- function(object, expected_string, ...)
{
    regexp <- gsub("[\\s]+", "[\\\\s]+", expected_string, perl=TRUE)
    expect_error(object, regexp=regexp, perl=TRUE, ...)
}

### Fix human aux data on-the-fly.
### We know that NCBI originally messed up with the 'extra_bps' value
### for alleles IGHJ6*02 and IGHJ6*03 in the original human_gl.aux. They
### corrected this later in the updated human_gl.aux that they released
### in April 2025. We do our own correction here.
load_and_fix_human_auxdata <- function()
{
    auxdata <- load_auxdata("human", which="original")
    fixme <- auxdata[ , "allele_name"] %in% c("IGHJ6*02", "IGHJ6*03")
    auxdata[fixme, "extra_bps"] <- 1L  # replace 0L with 1L
    auxdata
}

