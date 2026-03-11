### Run this code **in this folder** to generate all the FASTA files in the
### strain-specific subfolders.

library(igblastr)

STRAIN_SPECIFIC_GERMLINE_SETS <- list(
    CAST_EiJ  =c(IGH=1, IGKV=1, IGLV=1),
    LEWES_EiJ =c(IGH=1, IGKV=1, IGLV=1),
    MSM_MsJ   =c(IGH=1, IGKV=1, IGLV=1),
    NOD_ShiLtJ=c(IGH=1, IGKV=1, IGLV=1),
    PWD_PhJ   =c(IGH=2, IGKV=1, IGLV=1)
)

download_mouse_germline_sequences <- function(overwrite=FALSE)
{
    shared_germline_sets <- c(`IGKJ (all strains)`=1, `IGLJ (all strains)`=1)
    expected_files <- paste0(igblastr:::IG_GROUPS, ".fasta")
    for (i in seq_along(STRAIN_SPECIFIC_GERMLINE_SETS)) {
        destdir <- names(STRAIN_SPECIFIC_GERMLINE_SETS)[[i]]
        strain <- chartr("_", "/", destdir)
        message("Generating FASTA files for ", strain, " ... ", appendLF=FALSE)
        germline_sets <- STRAIN_SPECIFIC_GERMLINE_SETS[[i]]
        names(germline_sets) <- paste(strain, names(germline_sets))
        germline_sets <- c(germline_sets, shared_germline_sets)
        filenames <- download_OGRDB_germline_sequences("Mus musculus",
                                    germline_sets,
                                    destdir=destdir, overwrite=overwrite)
        stopifnot(length(filenames) == length(expected_files),
                  setequal(filenames, expected_files))
        message("ok")
    }
}

download_mouse_germline_sequences()


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### validate_mouse_intdata()
###
### Validates the intdata associated with the mouse germline sets listed
### in STRAIN_SPECIFIC_GERMLINE_SETS by comparing the two methods of
### acquisition:
###   1. Infer intdata from the gaps in the V allele sequences.
###   2. Extract intdata from OGRDB json file.
###
### WARNING: It looks like the following mouse germline sets have inconsistent
### internal data:
###
###   o "CAST/EiJ IGH"  version 1
###     Alleles with differences:
###     - IGHV0-NTLL*00: fwr3_end is 285 when inferred from the gaps,
###                      but is 283 in the JSON file;
###
###   o "CAST/EiJ IGLV" version 1
###     Alleles with differences:
###     - IGLV0-DUHW*00, IGLV0-EYCQ*00, IGLV0-JEYS*00, IGLV0-RQWY*00,
###       IGLV0-ZQBR*00: fwr3_end is 273 when inferred from the gaps,
###                      but is 245 in the JSON file;
###
###   o "LEWES/EiJ IGH" version 1
###     Alleles with differences:
###     - IGHV0-MJAE*00: fwr3_end is 288 when inferred from the gaps,
###                      but is 284 in the JSON file;
###
###   o "MSM/MsJ IGLV"  version 1
###     Alleles with differences:
###     - IGLV0-37NY*00, IGLV0-COOR*00: fwr3_end is 273 when inferred from
###                      the gaps, but is 245 in the JSON file.

validate_mouse_intdata <- function()
{
    organism <- "Mus musculus"
    for (i in seq_along(STRAIN_SPECIFIC_GERMLINE_SETS)) {
        germline_sets <- STRAIN_SPECIFIC_GERMLINE_SETS[[i]]
        strain <- chartr("_", "/", names(STRAIN_SPECIFIC_GERMLINE_SETS)[[i]])
        names(germline_sets) <- paste(strain, names(germline_sets))
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

