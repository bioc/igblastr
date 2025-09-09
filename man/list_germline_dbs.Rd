\name{list_germline_dbs}

\alias{list_germline_dbs}
\alias{print.germline_dbs_df}
\alias{use_germline_db}
\alias{load_germline_db}

\title{List cached germline dbs and select one to use with igblastn()}

\description{
  Use \code{list_germline_dbs()} to list all the \emph{cached germline dbs},
  that is, all the germline databases currently installed in \pkg{igblastr}'s
  persistent cache.

  Use \code{use_germline_db()} to select the cached germline db to use with
  \code{\link{igblastn}()}. This choice will be remembered for the duration
  of the current R session but can be changed anytime.

  Use \code{load_germline_db()} to load the nucleotide sequences of the
  gene regions stored in a cached germline db.
}

\usage{
list_germline_dbs(builtin.only=FALSE, names.only=FALSE, long.listing=FALSE)

use_germline_db(db_name=NULL, verbose=FALSE)

load_germline_db(db_name, region_types=NULL)
}

\arguments{
  \item{builtin.only}{
    By default \code{list_germline_dbs()} returns the list of all cached
    germline dbs, including built-in germline dbs. Set \code{builtin.only}
    to \code{TRUE} to return only the list of built-in germline dbs.
    Note that built-in dbs are prefixed with an underscore (\code{_}).
  }
  \item{names.only}{
    By default \code{list_germline_dbs()} returns the list of cached
    germline dbs in a data.frame with one db per row. Set \code{names.only}
    to \code{TRUE} to return only the db names in a character vector.
  }
  \item{long.listing}{
    \code{TRUE} or \code{FALSE}. If set to \code{TRUE}, then
    \code{list_germline_dbs()} returns a named list with one
    list element per germline db. Each list element is an integer
    matrix that indicates the number of germline sequences per locus
    and region type.

    Ignored if \code{names.only} is set to \code{TRUE}.
  }
  \item{db_name}{
    For \code{use_germline_db()}:

    \code{NULL} or a single string specifying the name of the
    cached germline db to use. Use \code{list_germline_dbs()} to
    list all the cached germline dbs.

    If set to \code{NULL} (the default), then
    \code{use_germline_db()} returns the name of the cached
    germline db that is currently in use, if any. Otherwise it
    raises an error.

    For \code{load_germline_db()}:

    A single string specifying the name of the cached germline db
    from which to load the V, D, and/or J regions. Use
    \code{list_germline_dbs()} to list all the cached germline dbs.
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

\details{
  The \pkg{igblastr} package provides a small set of utilities to
  manage the cached germline and C-region databases to use with
  \code{\link{igblastn}()}.

  Terminology:
  \itemize{
    \item A \emph{cached germline db} contains the nucleotide sequences of
          the V, D, and J gene regions for a given organism.
    \item A \emph{cached C-region db} contains the nucleotide sequences of
          the C regions (i.e. constant gene regions) for a given organism.
  }

  This man page documents the basic utilities to operate on cached
  germline dbs: \code{list_germline_dbs()}, \code{use_germline_db()},
  and \code{load_germline_db()}.

  The basic utilities to operate on cached C-region dbs are documented
  in the man page for \code{\link{list_c_region_dbs}}.
}

\value{
  \code{list_germline_dbs()} returns the list of all cached germline dbs
  in a data.frame with one db per row (if \code{names.only} is \code{FALSE},
  which is the default), or in a character vector (if \code{names.only}
  is \code{TRUE}). Columns V, D, J in the data.frame indicate the number
  of germline sequences for each region in each db.

  Built-in dbs are prefixed with an underscore (\code{_}).
  Note that the germline built-in dbs from AIRR were obtained from
  \url{https://ogrdb.airr-community.org/germline_sets/Homo\%20sapiens}
  and \url{https://ogrdb.airr-community.org/germline_sets/Mus\%20musculus}
  and included in the \pkg{igblastr} package on the date indicated by the
  suffix of the db name.

  When called with no argument, \code{use_germline_db()} returns a single
  string containing the name of the cached germline db currently used by
  \code{\link{igblastn}()} if any, or it raises an error if no germline db
  has been selected yet.

  When called with the \code{db_name} argument, \code{use_germline_db(db_name)}
  returns \code{db_name} invisibly.

  \code{load_germline_db()} returns the nucleotide sequences from the
  specified germline db in a named \link[Biostrings]{DNAStringSet} object.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

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

## Get list of built-in germline dbs only.
list_germline_dbs(builtin.only=TRUE)
list_germline_dbs(builtin.only=TRUE, names.only=TRUE)  # db names only

## Long listing:
list_germline_dbs(long.listing=TRUE)

if (IMGT_is_up()) {
  ## Install Mouse germline db from IMGT (Perl required!):
  install_IMGT_germline_db("202518-3", "Mus_musculus", force=TRUE)

  list_germline_dbs()  # all germline dbs

  ## Select germline db to use with igblastn():
  db_name <- "IMGT-202518-3.Mus_musculus.IGH+IGK+IGL"
  use_germline_db(db_name)  # select germline db to use

  use_germline_db()  # get current selection

  ## Load germline sequences:
  load_germline_db(db_name)
  load_germline_db(db_name, region_types="D")
  load_germline_db(db_name, region_types="DJ")
}
}

\keyword{utilities}
