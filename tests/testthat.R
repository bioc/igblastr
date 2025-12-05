library(testthat)

is_bioc_build_machine <- isTRUE(as.logical(Sys.getenv("IS_BIOC_BUILD_MACHINE")))
if (!is_bioc_build_machine) {
    ## We temporarily set the cache to a different (non-persistent)
    ## location to prevent the tests from messing up the real cache
    ## located at 'R_user_dir("igblastr", "cache")'. Additionally, this
    ## allows running the tests in different R sessions without the
    ## concurrent runs stepping on each other's toes.
    ## Make sure to do this BEFORE loading igblastr.
    tmp_cache <- tempfile("igblastr_tests_cache_")
    options(igblastr_cache=tmp_cache)
}

library(igblastr)

if (!has_igblast()) install_igblast()
test_check("igblastr")

if (!is_bioc_build_machine) {
    options(igblastr_cache=NULL)
}

