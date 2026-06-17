\name{list_c_region_dbs}

\alias{list_c_region_dbs}
\alias{print.c_region_dbs_df}
\alias{rm_c_region_db}

\title{List cached C-region dbs}

\description{
  List the \emph{cached C-region dbs}, that is, the IgBLAST-compatible
  C-region databases currently installed in \pkg{igblastr}'s persistent cache.
}

\usage{
## List the cached C-region dbs:
list_c_region_dbs(builtin.only=FALSE, names.only=FALSE, long.listing=FALSE)

## Remove a C-region db from igblastr's persistent cache:
rm_c_region_db(db_name)
}

\arguments{
  \item{builtin.only}{
    By default \code{list_c_region_dbs()} returns the list of all cached
    C-region dbs, including built-in C-region dbs. Set \code{builtin.only}
    to \code{TRUE} to return only the list of built-in C-region dbs.
    Note that built-in dbs are prefixed with an underscore (\code{_}).
  }
  \item{names.only}{
    By default \code{list_c_region_dbs()} returns the list of cached
    C-region dbs in a data.frame with one db per row. Set \code{names.only}
    to \code{TRUE} to return the db names only, in which case they are
    returned in a character vector.
  }
  \item{long.listing}{
    \code{TRUE} or \code{FALSE}. If set to \code{TRUE}, then
    \code{list_c_region_dbs()} returns a named list with one
    list element per C-region db. Each list element is a named integer
    vector that indicates the number of C-region sequences per locus.

    Ignored if \code{names.only} is set to \code{TRUE}.
  }
  \item{db_name}{
    A single string specifying the name of the C-region db to remove
    from \pkg{igblastr}'s persistent cache. This cannot be a built-in db.
  }
}

\details{
  \subsection{Cached germline dbs and C-region dbs}{
    The \pkg{igblastr} package provides utility functions to manage the
    \emph{cached germline dbs} and \emph{cached C-region dbs} to use with
    \code{\link{igblastn}()}.

    Terminology used across \pkg{igblastr} documentation:
    \itemize{
      \item A \emph{cached germline db} contains the nucleotide sequences of
            the germline V, D, and J gene alleles for a given organism.
      \item A \emph{cached C-region db} contains the nucleotide sequences of
            the constant gene regions for a given organism.
    }

    This man page documents utilities that manage cached C-region dbs.

    See \code{?\link{list_germline_dbs}} for utilities that manage
    cached germline dbs.
  }
  \subsection{Built-in dbs}{
    The \pkg{igblastr} package comes with preinstalled cached germline and
    C-region dbs that are called \emph{built-in dbs}.

    Built-in dbs have their name prefixed with an underscore (\code{_}).

    Note that the built-in C-region dbs from IMGT were downloaded from
    \url{https://www.imgt.org/vquest/refseqh.html#constant-sets} and
    included in the \pkg{igblastr} package on the date indicated by the
    suffix of the db name.
  }
}

\value{
  \code{list_c_region_dbs()} returns the list of all cached C-region dbs
  in a data.frame with one db per row (if \code{names.only} is \code{FALSE},
  which is the default), or in a character vector (if \code{names.only}
  is \code{TRUE}). Column C in the data.frame indicates the number of
  C-region sequences in each db.

  \code{rm_c_region_db()} does not return anything (invisible \code{NULL}).
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{use_c_region_db}} to select the cached C-region db
          to use with \code{igblastn()}.

    \item \code{\link{reset_c_region_dbs}} to reset the cached C-region dbs.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## 7 built-in C-region dbs (prefixed with an underscore):
list_c_region_dbs()
list_c_region_dbs(names.only=TRUE)    # db names only
list_c_region_dbs(long.listing=TRUE)  # long listing
}

\keyword{utilities}
