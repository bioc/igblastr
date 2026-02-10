### =========================================================================
### use_c_region_db() and related
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### use_c_region_db()
###

### Returns "" if no db is currently in use.
.get_selected_c_region_db <- function(verbose=FALSE)
{
    c_region_dbs_home <- get_c_region_dbs_home(TRUE)  # guaranteed to exist
    db_path <- get_db_in_use(c_region_dbs_home, what="C-region")
    if (db_path == "")
        return(db_path)
    make_blastdbs(db_path, verbose=verbose)
    basename(db_path)
}

.how_to_suppress_use_c_region_db_msg <- function(db_name)
{
    msg1 <- "To suppress this message, use:"
    msg2 <- c("suppressMessages(use_c_region_db(\"", db_name, "\"))")
    c(wmsg(msg1), "\n    ", wmsg(msg2))
}

.note_on_selecting_IMGT_c_region_db <- function(db_name)
{
    is_imgt_db <- grepl("^_?IMGT\\.", db_name)
    if (!is_imgt_db)
        return()
    message("  ", wmsg(IMGT_TERMS_OF_USE), "\n\n  ",
            .how_to_suppress_use_c_region_db_msg(db_name))
}

.select_c_region_db <- function(db_name, verbose=FALSE)
{
    if (!isSingleString(db_name))
        stop(wmsg("'db_name' must be a single string"))

    if (db_name != "") {
        check_c_region_db_name(db_name)
        .note_on_selecting_IMGT_c_region_db(db_name)
        db_path <- get_c_region_db_path(db_name)
        make_blastdbs(db_path, verbose=verbose)
    }

    ## Returns 'db_name' invisibly.
    set_db_in_use("C-region", db_name, verbose=verbose)
}

### Passing 'db_name=""' will cancel the current selection.
use_c_region_db <- function(db_name=NULL, verbose=FALSE)
{
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))
    if (is.null(db_name))
        return(.get_selected_c_region_db(verbose=verbose))
    .select_c_region_db(db_name, verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### clean_c_region_blastdbs()
###

### Not used at the moment and not exported!
clean_c_region_blastdbs <- function()
{
    c_region_dbs_home <- get_c_region_dbs_home()  # NOT guaranteed to exist
    if (dir.exists(c_region_dbs_home)) {
        all_db_names <- list_c_region_dbs(names.only=TRUE)
        for (db_name in all_db_names) {
            db_path <- get_c_region_db_path(db_name)
            clean_blastdbs(db_path)
        }
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_c_region_db()
###

### Returns the C regions in a DNAStringSet object.
load_c_region_db <- function(db_name)
{
    check_c_region_db_name(db_name)
    db_path <- get_c_region_db_path(db_name)
    db_fasta_file <- get_db_fasta_file(db_path, "C")
    readDNAStringSet(db_fasta_file)
}

