\name{list_c_region_dbs}

\alias{list_c_region_dbs}
\alias{print.c_region_dbs_df}
\alias{use_c_region_db}
\alias{load_c_region_db}

\title{List local C-region dbs and select one to use with igblastn()}

\description{
  Use \code{list_c_region_dbs()} to list all the \emph{local C-region dbs},
  that is, all the C-region databases currently installed in \pkg{igblastr}'s
  persistent cache.

  Use \code{use_c_region_db()} to select the local C-region db to use with
  \code{\link{igblastn}()}. This choice will be remembered for the duration
  of the current R session but can be changed anytime.

  Use \code{load_c_region_db()} to load the nucleotide sequences of the
  gene regions stored in a local C-region db.
}

\usage{
list_c_region_dbs(builtin.only=FALSE, names.only=FALSE, long.listing=FALSE)

use_c_region_db(db_name=NULL, verbose=FALSE)

load_c_region_db(db_name)
}

\arguments{
  \item{builtin.only}{
    By default \code{list_c_region_dbs()} returns the list of all local
    C-region dbs, including built-in C-region dbs. Set \code{builtin.only}
    to \code{TRUE} to return only the list of built-in C-region dbs.
    Note that built-in dbs are prefixed with an underscore (\code{_}).
  }
  \item{names.only}{
    By default \code{list_c_region_dbs()} returns the list of local
    C-region dbs in a data.frame with one db per row. Set \code{names.only}
    to \code{TRUE} to return only the db names in a character vector.
  }
  \item{long.listing}{
    \code{TRUE} or \code{FALSE}. If set to \code{TRUE}, then
    \code{list_c_region_dbs()} returns a named list with one
    list element per C-region db. Each list element is a named integer
    vector that indicates the number of C-region sequences per locus.

    Ignored if \code{names.only} is set to \code{TRUE}.
  }
  \item{db_name}{
    For \code{use_c_region_db()}:

    \code{NULL} or a single string specifying the name of the
    local C-region db to use. Use \code{list_c_region_dbs()} to
    list all the local C-region dbs.

    If set to \code{NULL} (the default), then
    \code{use_c_region_db()} returns the name of the local
    C-region db that is currently in use, if any. Otherwise it
    returns the empty string (\code{""}).

    Note that the current selection can be cancelled with
    \code{use_c_region_db("")}.

    For \code{load_c_region_db()}:

    A single string specifying the name of the local C-region db
    from which to load the gene regions. Use \code{list_c_region_dbs()}
    to list all the local C-region dbs.
  }
  \item{verbose}{
    If set to \code{TRUE}, then \code{use_c_region_db()} will display
    some information about its internal operations.
  }
}

\details{
  The \pkg{igblastr} package provides a small set of utilities to
  manage the local germline and C-region databases to use with
  \code{\link{igblastn}()}.

  Terminology:
  \itemize{
    \item A \emph{local germline db} contains the nucleotide sequences of
          the V, D, and J gene regions for a given organism.
    \item A \emph{local C-region db} contains the nucleotide sequences of
          the C regions (i.e. constant gene regions) for a given organism.
  }

  This man page documents the basic utilities to operate on local
  C-region dbs: \code{list_c_region_dbs()}, \code{use_c_region_db()},
  and \code{load_c_region_db()}.

  The basic utilities to operate on local germline dbs are documented
  in the man page for \code{\link{list_germline_dbs}}.
}

\value{
  \code{list_c_region_dbs()} returns the list of all local C-region dbs
  in a data.frame with one db per row (if \code{names.only} is \code{FALSE},
  which is the default), or in a character vector (if \code{names.only}
  is \code{TRUE}). Column C in the data.frame indicates the number of
  C-region sequences in each db.

  Built-in dbs are prefixed with an underscore (\code{_}).
  Note that the C-region built-in dbs from IMGT were downloaded from
  \url{https://www.imgt.org/vquest/refseqh.html#constant-sets} and
  included in the \pkg{igblastr} package on the date indicated by the
  suffix of the db name.

  When called with no argument, \code{use_c_region_db()} returns a single
  string containing the name of the local C-region db currently used by
  \code{\link{igblastn}()} if any, or the empty string (\code{""}) if
  \code{\link{igblastn}()} is not using any C-region db.

  When called with the \code{db_name} argument, \code{use_c_region_db(db_name)}
  returns \code{db_name} invisibly.

  \code{load_c_region_db()} returns the nucleotide sequences from the
  specified C-region db in a named \link[Biostrings]{DNAStringSet} object.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{use_germline_db}} to select the local germline db
          to use with \code{igblastn()}.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## 3 built-in C-region dbs prefixed with an underscore:
list_c_region_dbs()

list_c_region_dbs(names.only=TRUE)  # db names only

list_c_region_dbs(long.listing=TRUE)  # long listing

## Select C-region db to use with igblastn():
use_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")

use_c_region_db()  # get current selection

use_c_region_db("")  # cancel current selection

use_c_region_db()

## Load C-region sequences:
load_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")
load_c_region_db("_IMGT.mouse.IGH.202412")
}

\keyword{utilities}
