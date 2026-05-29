\name{use_germline_dbs}

\alias{use_germline_db}
\alias{load_germline_db}

\title{Select cached germline db to use with igblastn()}

\description{
  \code{use_germline_db()} allows the user to select the cached germline
  db to use with \code{\link{igblastn}()}. This choice will be remembered
  for the duration of the current R session but can be changed anytime.

  \code{load_germline_db()} allows the user to load the nucleotide sequences
  of the germline gene alleles stored in a cached germline db.
}

\usage{
use_germline_db(db_name=NULL, verbose=FALSE)

load_germline_db(db_name, region_types=NULL)
}

\arguments{
  \item{db_name}{
    For \code{use_germline_db()}:

    \code{NULL} or a single string specifying the name of the
    cached germline db to use. Use \code{\link{list_germline_dbs}()}
    to list all the cached germline dbs.

    If set to \code{NULL} (the default), then
    \code{use_germline_db()} returns the name of the cached
    germline db that is currently in use, if any. Otherwise it
    raises an error.

    For \code{load_germline_db()}:

    A single string specifying the name of the cached germline db
    from which to load the V, D, and/or J allele sequences. Use
    \code{\link{list_germline_dbs}()} to list all the cached germline dbs.
  }
  \item{verbose}{
    If set to \code{TRUE}, then \code{use_germline_db()} will display
    some information about its internal operations.
  }
  \item{region_types}{
    The types of regions (V, D, and/or J) to load from the database.
    Specified as a single string (e.g. \code{"DJ"}) or as a character
    vector of single-letter elements (e.g. \code{c("D", "J")}).
    By default (i.e. when \code{region_types} is \code{NULL}), all the
    regions are returned.
  }
}

\value{
  When called with no argument, \code{use_germline_db()} returns a single
  string containing the name of the cached germline db currently used by
  \code{\link{igblastn}()} if any, or it raises an error if no germline db
  has been selected yet.

  When called with the \code{db_name} argument, \code{use_germline_db(db_name)}
  returns \code{db_name} invisibly.

  \code{load_germline_db()} returns the allele sequences from the
  specified germline db in a named \link[Biostrings]{DNAStringSet} object.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \code{\link{use_c_region_db}} to select the cached C-region db
          to use with \code{igblastn()}.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

list_germline_dbs()  # all cached germline dbs

## Select germline db to use with igblastn():
db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"
use_germline_db(db_name)

## Get current selection:
use_germline_db()

## Load germline gene allele sequences:
load_germline_db(db_name)
load_germline_db(db_name, region_types="D")
load_germline_db(db_name, region_types="DJ")
}

\keyword{utilities}
