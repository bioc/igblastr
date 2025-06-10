.onLoad <- function(libname, pkgname)
{
    igblastr_usage_report <- getOption("igblastr_usage_report")
    if (is.null(igblastr_usage_report)) {
        igblastr_usage_report <-
            get_igblastr_usage_report_from_BLAST_USAGE_REPORT()
        options(igblastr_usage_report=igblastr_usage_report)
    }

    ## Removing the blast dbs from all the local germline and C-region dbs
    ## at load-time is actually a terrible idea because it will pull the rug
    ## out from under any other R session currently using them!
    #clean_germline_blastdbs()  # a VERY BAD idea!
    #clean_c_region_blastdbs()  # a VERY BAD idea!
}

