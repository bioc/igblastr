### =========================================================================
### R wrapper to the edit_imgt_file.pl script included in IgBLAST
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_edit_imgt_file_Perl_script()
###

### Requires Perl!
### Checks that Perl script edit_imgt_file.pl is available and that Perl
### is functioning.
get_edit_imgt_file_Perl_script <- function()
{
    igblast_root <- get_igblast_root()
    bin_dir <- get_igblast_root_subdir(igblast_root, "bin")
    script <- file.path(bin_dir, "edit_imgt_file.pl")
    if (!file.exists(script)) {
        details <- c("Perl script 'edit_imgt_file.pl' (needed by ",
                     "internal helper edit_imgt_file()) not found ",
                     "in 'bin' subdirectory.")
        stop_on_invalid_igblast_root(igblast_root, details)
    }
    if (!has_perl())
        stop(wmsg("Setup error: Perl not found."))
    script
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### edit_imgt_file()
###

### Requires Perl!
edit_imgt_file <- function(infasta, outfasta, errfile=NULL, Perl_script=NULL)
{
    if (!isSingleNonWhiteString(infasta))
        stop(wmsg("'infasta' must be a single (non-empty) string"))
    if (!isSingleNonWhiteString(outfasta))
        stop(wmsg("'outfasta' must be a single (non-empty) string"))
    if (is.null(errfile)) {
        errfile <- tempfile("edit_imgt_file_errors", fileext=".txt")
        on.exit(unlink(errfile))
    } else if (!isSingleNonWhiteString(errfile)) {
        stop(wmsg("'errfile' must be NULL or a single (non-empty) string"))
    }
    if (is.null(Perl_script)) {
        Perl_script <- get_edit_imgt_file_Perl_script()
    } else if (!isSingleNonWhiteString(Perl_script)) {
        stop(wmsg("'Perl_script' must be NULL or a single (non-empty) string"))
    }

    ## This does not work on Windows!
    #system3(Perl_script, outfasta, errfile, args=infasta)

    ## Note that running the Perl script with 'script ...' works on Linux
    ## and Mac but not on Windows. So we run it with 'perl script ...'
    ## instead. This seems to run everywhere.
    system3("perl", outfasta, errfile, args=c(Perl_script, infasta))
}

