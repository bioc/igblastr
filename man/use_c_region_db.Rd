\name{use_c_region_db}

\alias{use_c_region_db}
\alias{load_c_region_sequences}
\alias{load_c_region_db}

\title{Select cached C-region db to use with igblastn()}

\description{
  \code{use_c_region_db()} allows the user to select the cached C-region
  db to use with \code{\link{igblastn}()}. This choice will be remembered
  for the duration of the current R session but can be changed anytime.

  \code{load_c_region_sequences()} allows the user to load the nucleotide
  sequences of the gene alleles stored in a cached C-region db.
}

\usage{
use_c_region_db(db_name=NULL, verbose=FALSE)

load_c_region_sequences(db_name)
}

\arguments{
  \item{db_name}{
    For \code{use_c_region_db()}:

    \code{NULL} or a single string specifying the name of the
    cached C-region db to use. Use \code{\link{list_c_region_dbs}()}
    get the list of all cached C-region dbs.

    If set to \code{NULL} (the default), then
    \code{use_c_region_db()} returns the name of the cached
    C-region db that is currently in use, if any. Otherwise it
    returns the empty string (\code{""}).

    Note that the current selection can be cancelled with
    \code{use_c_region_db("")}.

    For \code{load_c_region_sequences()}:

    A single string specifying the name of the cached C-region db
    from which to load the allele sequences. Use
    \code{\link{list_c_region_dbs}()} to get the list of all cached
    C-region dbs.
  }
  \item{verbose}{
    If set to \code{TRUE}, then \code{use_c_region_db()} will display
    some information about its internal operations.
  }
}

\value{
  When called with no argument, \code{use_c_region_db()} returns a single
  string containing the name of the cached C-region db currently used by
  \code{\link{igblastn}()} if any, or the empty string (\code{""}) if
  \code{\link{igblastn}()} is not using any C-region db.

  When called with the \code{db_name} argument, \code{use_c_region_db(db_name)}
  returns \code{db_name} invisibly.

  \code{load_c_region_sequences()} returns the allele sequences from the
  specified C-region db in a named \link[Biostrings]{DNAStringSet} object.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{list_c_region_dbs}} to list the cached C-region dbs.

    \item \code{\link{use_germline_db}} to select the cached germline db
          to use with \code{igblastn()}.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

list_c_region_dbs()  # all cached C-region dbs

## Select C-region db to use with igblastn():
db_name <- "_IMGT.human.IGH+IGK+IGL.202605"
use_c_region_db(db_name)

## Get current selection:
use_c_region_db()    # get current selection

## Cancel current selection:
use_c_region_db("")
use_c_region_db()

## Load C-region sequences:
load_c_region_sequences(db_name)
load_c_region_sequences("_IMGT.mouse.IGH+IGK+IGL.202605")
}

\keyword{utilities}
