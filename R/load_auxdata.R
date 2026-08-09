### =========================================================================
### load_auxdata() and related
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_germline_db_auxdata_path()
###

### Not exported!
make_germline_db_auxdata_path <- function(db_path)
{
    stopifnot(isSingleNonWhiteString(db_path), dir.exists(db_path))
    ## gl.aux: germline auxiliary file ("gl" in gl.aux stands for germline).
    file.path(db_path, "auxiliary_data", "gl.aux")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_auxdata_path()
###

.get_igblast_auxdata_path <- function(which, igblast_organism)
{
    igblast_organism <- normalize_igblast_organism(igblast_organism)
    auxdata_dir <- file.path(path_to_igdata(which), "optional_file")
    auxdata_filename <- paste0(igblast_organism, "_gl.aux")
    auxdata_path <- file.path(auxdata_dir, auxdata_filename)
    if (!file.exists(auxdata_path))
        stop(wmsg("no auxiliary data found in ",
                  auxdata_dir, " for ", igblast_organism))
    auxdata_path
}

get_auxdata_path <- function(igblast_organism, which=c("live", "original"))
{
    if (!isSingleNonWhiteString(igblast_organism))
        stop(wmsg("'igblast_organism' must be a single (non-empty) string"))
    which <- match.arg(which)

    if (!germline_db_exists(igblast_organism))
        return(.get_igblast_auxdata_path(which, igblast_organism))

    ## Treat 'igblast_organism' as the name of an existing germline db.
    db_path <- get_germline_db_path(igblast_organism)
    auxdata_path <- make_germline_db_auxdata_path(db_path)

    if (!dir.exists(dirname(auxdata_path)))
        stop(wmsg("no auxiliary data found in germline db ", igblast_organism))
    if (!file.exists(auxdata_path))
        stop(wmsg("auxiliary data file ", basename(auxdata_path), " ",
                  "not found in germline db ", igblast_organism))
    auxdata_path
}

get_igblast_auxiliary_data <- function(...)
{
    .Deprecated("get_auxdata_path")
    get_auxdata_path(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_auxdata()
###

### IMPORTANT NOTE: Unlike with the data.frame returned by load_intdata(),
### all the positions in the data.frame returned by load_auxdata() (that is,
### the positions reported in columns 'coding_frame_start' and 'cdr3_end')
### are 0-based!
load_auxdata <- function(igblast_organism, which=c("live", "original"))
{
    which <- match.arg(which)
    auxdata_path <- get_auxdata_path(igblast_organism, which=which)
    read_auxdata(auxdata_path)
}

load_igblast_auxiliary_data <- function(...)
{
    .Deprecated("load_auxdata")
    load_auxdata(...)
}

