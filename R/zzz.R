.onLoad <- function(libname, pkgname)
{
    igblastr_usage_report <- getOption("igblastr_usage_report")
    if (is.null(igblastr_usage_report)) {
        igblastr_usage_report <-
            get_igblastr_usage_report_from_BLAST_USAGE_REPORT()
        options(igblastr_usage_report=igblastr_usage_report)
    }

    ## https://www.imgt.org can sometimes be really slow to answer so
    ## we'll set curl configuration option 'connecttimeout' to 20s when
    ## we query it with httr::HEAD() or httr::GET() (looks like these
    ## functions give up after 10s by default).
    IMGT_connecttimeout <- getOption("IMGT_connecttimeout")
    if (is.null(IMGT_connecttimeout))
        options(IMGT_connecttimeout=20)  # 20 seconds

    ## Removing the blast dbs from all the cached germline and C-region dbs
    ## at load-time is actually a terrible idea because it will pull the rug
    ## out from under any other R session currently using them!
    #clean_germline_blastdbs()  # a VERY BAD idea!
    #clean_c_region_blastdbs()  # a VERY BAD idea!
}

