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

make_rhesus_monkey_auxdata <- function(version, germline_sets)
{
    dir.create(jsondir <- tempfile())
    message("Creating IG[HKL]J_gl.aux files in ", version, " ... ",
            appendLF=FALSE)
    filenames <- download_OGRDB_germline_json("Macaca mulatta",
                                germline_sets,
                                destdir=jsondir, overwrite=TRUE)
    stopifnot(identical(names(filenames), names(germline_sets)))
    for (filename in filenames) {
        locus <- substr(filename, 1L, 3L)
        json_path <- file.path(jsondir, filename)
        if (locus == "IGH") {
            ## cdr3_end is -1 for rhesus monkey allele IGHJ0-ZXTW*01.
            ## This triggers a warning that we suppress. Also write_auxdata()
            ## does not allow negative values so we replace it with NA.
            auxdata <- suppressWarnings(
                extract_auxdata_from_ogrdb_json(json_path)
            )
            bad_idx <- which(auxdata[ , "cdr3_end"] < 0L)
            bad_alleles <- auxdata[bad_idx, "allele_name"]
            stopifnot(identical(bad_alleles, "IGHJ0-ZXTW*01"))
            auxdata[bad_idx, "cdr3_end"] <- NA_integer_
        } else {
            auxdata <- extract_auxdata_from_ogrdb_json(json_path)
        }
        destfile <- file.path(version, paste0(locus, "J_gl.aux"))
        write_auxdata(auxdata, destfile)
    }
    message("ok")
}

download_rhesus_monkey_germline_sequences()

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

