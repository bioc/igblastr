### Run this code **in this folder** to generate all the FASTA files in the
### various subfolders of this folder.

library(igblastr)

RHESUS_MONKEY_GERMLINE_SETS <- list(
    `202601`=c(IGH_VDJ=1, IGK_VJ=1, IGL_VJ=1),
    `202602`=c(IGH_VDJ=2, IGK_VJ=2, IGL_VJ=2)
)

download_rhesus_monkey_germline_sequences <- function(overwrite=FALSE)
{
    expected_files <- paste0(igblastr:::IG_GROUPS, ".fasta")
    for (i in seq_along(RHESUS_MONKEY_GERMLINE_SETS)) {
        destdir <- names(RHESUS_MONKEY_GERMLINE_SETS)[[i]]
        germline_sets <- RHESUS_MONKEY_GERMLINE_SETS[[i]]
        message("Creating FASTA files in ", destdir, " ... ",
                appendLF=FALSE)
        filenames <- download_OGRDB_germline_sequences("Macaca mulatta",
                                    germline_sets,
                                    destdir=destdir, overwrite=overwrite)
        stopifnot(length(filenames) == length(expected_files),
                  setequal(filenames, expected_files))
        message("ok")
    }
}

download_rhesus_monkey_germline_sequences()


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### validate_rhesus_monkey_intdata()
###
### Validates the intdata associated with the rhesus monkey germline sets
### listed in RHESUS_MONKEY_GERMLINE_SETS by comparing the two methods of
### acquisition:
###   1. Infer intdata from the gaps in the V allele sequences.
###   2. Extract intdata from OGRDB json file.
###
### Note that all rhesus monkey germline sets have consistent internal data.

validate_rhesus_monkey_intdata <- function()
{
    organism <- "Macaca mulatta"
    for (i in seq_along(RHESUS_MONKEY_GERMLINE_SETS)) {
        germline_sets <- RHESUS_MONKEY_GERMLINE_SETS[[i]]
        for (j in seq_along(germline_sets)) {
            germline_set <- germline_sets[j]
            what <- c(organism, " germline set ",
                      "\"", names(germline_set), "\" ",
                      "(version ", germline_set, ")")
            message("Validating intdata for ", what, " ... ",
                    appendLF=FALSE)
            ok <- igblastr:::validate_OGRDB_intdata(organism, germline_set)
            msg <- if (ok) "ok" else "DATA IS INCONSISTENT!"
            message(msg)
        }
    }
}

