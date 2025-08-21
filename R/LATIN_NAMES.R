### A named character vector that represents a list of key/value pairs,
### where the keys are the short names used by IgBLAST for the 5 organisms
### that it officially supports. See function list_igblast_organisms() in
### file igblast_info.R.

LATIN_NAMES <- c(
    human="Homo sapiens",
    mouse="Mus musculus",
    rabbit="Oryctolagus cuniculus",
    rat="Rattus norvegicus",
    rhesus_monkey="Macaca mulatta"
)

### Treats 'organism' as a regular expression.
### Returns one of 'names(LATIN_NAMES)'.
find_organism_shortname <- function(organism)
{
    stopifnot(isSingleNonWhiteString(organism))
    ans <- grep(chartr(" ", "_", organism), names(LATIN_NAMES),
                ignore.case=TRUE, value=TRUE)
    if (length(ans) == 1L)
        return(ans)
    if (length(ans) >= 2L)
        stop(wmsg("ambigous organism abbreviation: ", organism))
    idx <- grep(chartr("_", " ", organism), LATIN_NAMES,
                ignore.case=TRUE)
    if (length(idx) == 0L)
        stop(wmsg("unrecognized organism: ", organism))
    if (length(idx) != 1L)
        stop(wmsg("ambigous organism abbreviation: ", organism))
    names(LATIN_NAMES)[idx]
}

