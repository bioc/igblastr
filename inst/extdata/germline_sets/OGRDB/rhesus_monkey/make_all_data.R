### Run this code **in this folder** to generate all the FASTA files in the
### various subfolders of this folder.

library(igblastr)

RHESUS_MONKEY_GERMLINE_SETS <- list(
    `202601`=c(IGH_VDJ=1, IGK_VJ=1, IGL_VJ=1),
    `202602`=c(IGH_VDJ=2, IGK_VJ=2, IGL_VJ=2)
)

download_rhesus_monkey_germline_sequences <- function(overwrite=FALSE)
{
    expected_fasta_files <- paste0(igblastr:::IG_GROUPS, ".fasta")
    for (i in seq_along(RHESUS_MONKEY_GERMLINE_SETS)) {
        destdir <- names(RHESUS_MONKEY_GERMLINE_SETS)[[i]]
        germline_sets <- RHESUS_MONKEY_GERMLINE_SETS[[i]]
        message("Creating FASTA files in ", destdir, " ... ",
                appendLF=FALSE)
        filenames <- download_OGRDB_germline_sequences("Macaca mulatta",
                                    germline_sets,
                                    destdir=destdir, overwrite=overwrite)
        stopifnot(length(filenames) == length(expected_fasta_files),
                  setequal(filenames, expected_fasta_files))
        message("ok")
    }
}

download_rhesus_monkey_germline_sequences()

make_rhesus_monkey_auxdata <- function(version, germline_sets, overwrite=FALSE)
{
    dir.create(json_dir <- tempfile())
    message("Creating IG[HKL]J_gl.aux files in ", version, " ... ",
            appendLF=FALSE)
    json_files <- download_OGRDB_germline_json("Macaca mulatta",
                                 germline_sets,
                                 destdir=json_dir, overwrite=TRUE)
    stopifnot(identical(names(json_files), names(germline_sets)))
    alleles_with_neg_cdr3_end <- list(IGH="IGHJ0-ZXTW*01")
    auxdata_files <-
        make_auxdata_files_from_ogrdb_jsons(json_dir, destdir=version,
                           alleles_with_neg_cdr3_end=alleles_with_neg_cdr3_end,
                           overwrite=overwrite)
    expected_files <- sprintf("IG%sJ_gl.aux", c("H", "K", "L"))
    stopifnot(identical(auxdata_files, expected_files))
    message("ok")
}

### Extract auxdata from the OGRDB json files.
for (i in seq_along(RHESUS_MONKEY_GERMLINE_SETS)) {
    version <- names(RHESUS_MONKEY_GERMLINE_SETS)[[i]]
    germline_sets <- RHESUS_MONKEY_GERMLINE_SETS[[i]]
    make_rhesus_monkey_auxdata(version, germline_sets)
}


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

validate_rhesus_monkey_intdata()

