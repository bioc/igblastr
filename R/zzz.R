.onLoad <- function(libname, pkgname)
{
    ## Removing the blast dbs from all the local germline and C-region dbs
    ## at load-time is actually a terrible idea because it will pull the rug
    ## out from under any other R session currently using them!
    #clean_germline_blastdbs()  # a VERY BAD idea!
    #clean_c_region_blastdbs()  # a VERY BAD idea!
}

