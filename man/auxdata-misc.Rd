\name{auxdata-misc}

\alias{auxdata-misc}
\alias{auxdata_misc}

\alias{compute_germline_db_auxdata}

\title{Miscellaneous auxiliary data utilities}

\description{
  \code{compute_germline_db_auxdata()} can be used to compute
  the \emph{auxiliary data} associated with the J alleles in a
  cached germline dbs managed by \pkg{igblastr}.
}

\usage{
compute_germline_db_auxdata(db_name, ...)
}

\arguments{
  \item{db_name}{
    A single string specifying the name of a cached germline db.
    Use \code{\link{list_germline_dbs}()} to get the list of all cached
    germline dbs.
  }
  \item{...}{
    Extra arguments to be passed to the internal call to
    \code{compute_auxdata()}. See \code{?\link{compute_auxdata}} for
    what arguments are supported.
  }
}

\value{
  \code{compute_germline_db_auxdata()} returns the computed
  \emph{auxiliary data} in a data.frame with 1 row per J allele in
  the germline db and the same columns as the data.frame returned
  by \code{compute_auxdata()} or \code{read_auxdata()}.
  See \code{?\link{compute_auxdata}} for more information.
}

\seealso{
  \itemize{
    \item \code{\link{load_auxdata}} to access IgBLAST \emph{auxiliary data}.

    \item \code{\link{compute_auxdata}} to annotate a set of germline J
          gene allele sequences.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
          for important information about IgBLAST \emph{auxiliary data}.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## ---------------------------------------------------------------------
## Compute the auxiliary data associated with the J alleles in a
## germline db
## ---------------------------------------------------------------------

## The auxiliary data included in germline db
## _OGRDB.rhesus_monkey.IGH+IGK+IGL.202602 can be computed
## from scratch with compute_germline_db_auxdata():
db_name <- "_OGRDB.rhesus_monkey.IGH+IGK+IGL.202602"
computed_auxdata <- compute_germline_db_auxdata(db_name)

## Replace the negative "cdr3_end" value with NA:
replace_idx <- which(computed_auxdata[ , "cdr3_end"] < 0L)
computed_auxdata[replace_idx, "cdr3_end"] <- NA

## Sanity check:
stopifnot(identical(computed_auxdata, load_auxdata(db_name)))
}

\keyword{utilities}
