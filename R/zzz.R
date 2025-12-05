.print_msg_if_live_igdata_needs_check <- function(dt)
{
    stopifnot(isSingleNumber(dt))
    if (is.infinite(dt)) {
        msg <- c("You can run 'update_live_igdata()' to check for new ",
                 "IgBLAST auxiliary or internal data files available at ",
                 "NCBI. See '?update_live_igdata' for more information.")
    } else if (dt > 30) {
        msg <- c("More than 30 days have passed since the last time ",
                 "you ran 'update_live_igdata()'. Time to run it again!")
    } else {
        return(invisible(NULL))
    }
    packageStartupMessage(wmsg("igblastr tip: ", msg))
}

.onLoad <- function(libname, pkgname)
{
    if (!dir.exists(igblastr_cache(LIVE_IGDATA)))
        reset_live_igdata()

    dt <- time_since_live_igdata_last_checked()
    .print_msg_if_live_igdata_needs_check(dt)

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

