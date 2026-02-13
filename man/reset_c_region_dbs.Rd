\name{reset_c_region_dbs}

\alias{reset_c_region_dbs}

\title{Reset the cached C-region dbs}

\description{
  \code{reset_c_region_dbs()} will reset the cache of C-region dbs to a
  pristine state by:
  \enumerate{
    \item Removing \emph{all} the C-region dbs currently in \pkg{igblastr}'s
          persistent cache, including the built-in C-region dbs.
    \item Recreating the predefined set of built-in C-region dbs.
  }
  This means that any non built-in C-region db present in the cache WILL
  BE LOST!

  Also note that \code{reset_c_region_dbs()} will cancel the current
  selection if any C-region db was previously selected with
  \code{\link{use_c_region_db}()}.
}

\usage{
reset_c_region_dbs(verbose=FALSE)
}

\arguments{
  \item{verbose}{
    Set to \code{TRUE} to have \code{reset_c_region_dbs()} display
    some details about its internal operations.
  }
}

\value{
  Nothing (i.e. an invisible \code{NULL}).
}

\seealso{
  \itemize{
    \item \code{\link{list_c_region_dbs}} to list the cached C-region dbs.

    \item \code{\link{use_c_region_db}} to select the cached C-region db
          to use with \code{igblastn()}.

    \item \code{\link{install_IMGT_c_region_db}} to install a C-region db
          from IMGT.
  }
}

\examples{
## DANGER ZONE!
\dontrun{
  reset_c_region_dbs(verbose=TRUE)
}
}

\keyword{utilities}
