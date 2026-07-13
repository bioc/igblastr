### =========================================================================
### load_intdata() and related
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_germline_db_intdata_path()
### write_ndm_data_to_db()
###

.make_intdata_file_suffix <- function(for.aa, domain_system)
{
    stopifnot(isTRUEorFALSE(for.aa), isSingleNonWhiteString(domain_system))
    paste0(".", if (for.aa) "pdm" else "ndm", ".", domain_system)
}

### Not exported!
make_germline_db_intdata_path <- function(db_path, for.aa, domain_system)
{
    stopifnot(isSingleNonWhiteString(db_path), dir.exists(db_path))
    file_suffix <- .make_intdata_file_suffix(for.aa, domain_system)
    intdata_filename <- paste0("V", file_suffix)
    file.path(db_path, "internal_data", intdata_filename)
}

### Not exported!
### Note that we don't handle custom "pdm" data at the moment. So the
### function does not need a 'for.aa' argument. Also the name of the
### function reflects the fact that it only handles "ndm" data.
write_ndm_data_to_db <- function(ndm_data, db_path,
                                 domain_system=c("imgt", "kabat"),
                                 check.and.reorder=FALSE)
{
    stopifnot(is.data.frame(ndm_data), isTRUEorFALSE(check.and.reorder))
    domain_system <- match.arg(domain_system)
    intdata_path <- make_germline_db_intdata_path(db_path, FALSE, domain_system)
    intdata_dir <- dirname(intdata_path)
    stopifnot(!dir.exists(intdata_dir))

    ## Even though write_ndm_data() will call check_ndm_data_col2class()
    ## internally, we prefer to fail **before** creating the 'intdata_dir'
    ## folder.
    check_ndm_data_col2class(ndm_data)
    if (check.and.reorder) {
        db_V_fasta_file <- get_db_fasta_file(db_path, "V")
        db_V_allele_names <- names(fasta.seqlengths(db_V_fasta_file))
        ndm_data <- check_and_reorder_igdata_rows(ndm_data, db_V_allele_names)
    }
    stopifnot(dir.create(intdata_dir))
    write_ndm_data(ndm_data, intdata_path)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_intdata_path()
###

.get_igblast_intdata_path <- function(which, igblast_organism,
                                      for.aa, domain_system)
{
    igblast_organism <- normalize_igblast_organism(igblast_organism)
    intdata_dir <- file.path(path_to_igdata(which), "internal_data",
                             igblast_organism)
    if (!dir.exists(intdata_dir))
        stop(wmsg("no internal data found for organism ", igblast_organism))
    file_suffix <- .make_intdata_file_suffix(for.aa, domain_system)
    intdata_filename <- paste0(igblast_organism, file_suffix)
    intdata_path <- file.path(intdata_dir, intdata_filename)
    if (!file.exists(intdata_path))
        stop(wmsg("internal data file ", intdata_filename, " ",
                  "not found in ", intdata_dir))
    intdata_path
}

get_intdata_path <- function(igblast_organism, for.aa=FALSE,
                             domain_system=c("imgt", "kabat"),
                             which=c("live", "original"))
{
    if (!isSingleNonWhiteString(igblast_organism))
        stop(wmsg("'igblast_organism' must be a single (non-empty) string"))
    if (!isTRUEorFALSE(for.aa))
        stop(wmsg("'for.aa' must be TRUE or FALSE"))
    domain_system <- match.arg(domain_system)
    which <- match.arg(which)

    if (!valid_germline_db_name(igblast_organism))
        return(.get_igblast_intdata_path(which, igblast_organism,
                                         for.aa, domain_system))

    ## Treat 'igblast_organism' as a valid germline db name.
    db_path <- get_germline_db_path(igblast_organism)
    intdata_path <- make_germline_db_intdata_path(db_path,
                                                  for.aa, domain_system)
    if (!dir.exists(dirname(intdata_path)))
        stop(wmsg("no internal data found in germline db ", igblast_organism))
    if (!file.exists(intdata_path))
        stop(wmsg("internal data file ", basename(intdata_path), " ",
                  "not found in germline db ", igblast_organism))
    intdata_path
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_intdata()
###

### IMPORTANT NOTE: The FWR/CDR positions in the returned data.frame are
### 1-based while the coding frame start positions are 0-based!
load_intdata <- function(igblast_organism, for.aa=FALSE,
                         domain_system=c("imgt", "kabat"),
                         which=c("live", "original"))
{
    domain_system <- match.arg(domain_system)
    which <- match.arg(which)
    intdata_path <- get_intdata_path(igblast_organism, for.aa=for.aa,
                                     domain_system=domain_system,
                                     which=which)
    read_ndm_data(intdata_path)
}

