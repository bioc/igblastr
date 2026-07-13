### Run this code **in this folder** to generate all the FASTA files in the
### various subfolders of this folder.

library(igblastr)

RAINBOW_TROUT_GERMLINE_SETS <- list(
    `202606`=c(IGH_VDJ=1)
)

download_rainbow_trout_germline_sequences <- function(overwrite=FALSE)
{
    #expected_fasta_files <- paste0(igblastr:::IG_GROUPS, ".fasta")
    expected_fasta_files <- paste0("IGH", c("V", "D", "J"), ".fasta")
    for (i in seq_along(RAINBOW_TROUT_GERMLINE_SETS)) {
        destdir <- names(RAINBOW_TROUT_GERMLINE_SETS)[[i]]
        germline_sets <- RAINBOW_TROUT_GERMLINE_SETS[[i]]
        message("Creating FASTA files in ", destdir, " ... ",
                appendLF=FALSE)
        filenames <- download_OGRDB_germline_sequences("Oncorhynchus mykiss",
                                    germline_sets,
                                    destdir=destdir, overwrite=overwrite)
        stopifnot(length(filenames) == length(expected_fasta_files),
                  setequal(filenames, expected_fasta_files))
        message("ok")
    }
}

make_rainbow_trout_auxdata <- function(version, germline_sets)
{
    dir.create(jsondir <- tempfile())
    message("Creating IG[HKL]J_gl.aux files in ", version, " ... ",
            appendLF=FALSE)
    filenames <- download_OGRDB_germline_json("Oncorhynchus mykiss",
                                germline_sets,
                                destdir=jsondir, overwrite=TRUE)
    stopifnot(identical(names(filenames), names(germline_sets)))
    for (filename in filenames) {
        locus <- substr(filename, 1L, 3L)
        json_path <- file.path(jsondir, filename)
        auxdata <- extract_auxdata_from_ogrdb_json(json_path)
        destfile <- file.path(version, paste0(locus, "J_gl.aux"))
        write_auxdata(auxdata, destfile)
    }
    message("ok")
}

download_rainbow_trout_germline_sequences()

### Extract auxdata from the OGRDB json files.
for (i in seq_along(RAINBOW_TROUT_GERMLINE_SETS)) {
    version <- names(RAINBOW_TROUT_GERMLINE_SETS)[[i]]
    germline_sets <- RAINBOW_TROUT_GERMLINE_SETS[[i]]
    make_rainbow_trout_auxdata(version, germline_sets)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### validate_rainbow_trout_intdata()
###
### Validates the intdata associated with the rainbow trout germline sets
### listed in RAINBOW_TROUT_GERMLINE_SETS by comparing the two methods of
### acquisition:
###   1. Infer intdata from the gaps in the V allele sequences.
###   2. Extract intdata from OGRDB json file.
###
### Note that all rainbow trout germline sets have consistent internal data.

validate_rainbow_trout_intdata <- function()
{
    organism <- "Oncorhynchus mykiss"
    for (i in seq_along(RAINBOW_TROUT_GERMLINE_SETS)) {
        germline_sets <- RAINBOW_TROUT_GERMLINE_SETS[[i]]
        for (j in seq_along(germline_sets)) {
            germline_set <- germline_sets[j]
            what <- c(organism, " germline set ",
                      "\"", names(germline_set), "\" ",
                      "(version ", germline_set, ")")
            message("Validating intdata for ", what, " ... ",
                    appendLF=FALSE)
            ok <- igblastr:::validate_OGRDB_intdata(organism, germline_set,
                               fwrcdr_ends=igblastr:::RAINBOW_TROUT_FWRCDR_ENDS)
            msg <- if (ok) "ok" else "DATA IS INCONSISTENT!"
            message(msg)
        }
    }
}

validate_rainbow_trout_intdata()

