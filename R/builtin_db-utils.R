### =========================================================================
### Creation of built-in germline and C-region dbs
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_builtin_AIRR_germline_db()
###

.form_builtin_AIRR_germline_db_name <- function(organism, fasta_dir)
{
    stopifnot(isSingleNonWhiteString(organism))
    loci <- get_loci_from_input_germline_fasta_set(fasta_dir)
    version <- read_version_file(fasta_dir)
    sprintf("_AIRR.%s.%s.%s", organism, paste(loci, collapse="+"), version)
}

create_builtin_AIRR_germline_db <- function(fasta_dir, organism, destdir)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir))
    db_name <- .form_builtin_AIRR_germline_db_name(organism, fasta_dir)
    db_path <- file.path(destdir, db_name)
    create_germline_db(fasta_dir, db_path)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### create_builtin_IMGT_c_region_db()
###

.form_builtin_IMGT_c_region_db_name <- function(organism, fasta_dir,
                                                tcr.db=FALSE)
{
    stopifnot(isSingleNonWhiteString(organism))
    loci <- get_loci_from_input_c_region_fasta_set(fasta_dir, tcr.db=tcr.db)
    version <- read_version_file(fasta_dir)
    sprintf("_IMGT.%s.%s.%s", organism, paste(loci, collapse="+"), version)
}

create_builtin_IMGT_c_region_db <- function(fasta_dir, organism, destdir,
                                            tcr.db=FALSE)
{
    stopifnot(isSingleNonWhiteString(destdir), dir.exists(destdir))
    db_name <- .form_builtin_IMGT_c_region_db_name(organism, fasta_dir,
                                                   tcr.db=tcr.db)
    db_path <- file.path(destdir, db_name)
    create_c_region_db(fasta_dir, db_path, tcr.db=tcr.db)
}

