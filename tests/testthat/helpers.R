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

