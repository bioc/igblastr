\name{reset_germline_dbs}

\alias{reset_germline_dbs}

\title{Reset the cached germline dbs}

\description{
  \code{reset_germline_dbs()} will reset the cache of germline dbs to a
  pristine state by:
  \enumerate{
    \item Removing \emph{all} the germline dbs currently in \pkg{igblastr}'s
          persistent cache, including the built-in germline dbs.
    \item Recreating the predefined set of built-in germline dbs.
  }
  This means that any non built-in germline db present in the cache WILL
  BE LOST!
}

\usage{
reset_germline_dbs(verbose=FALSE)
}

\arguments{
  \item{verbose}{
    Set to \code{TRUE} to have \code{reset_germline_dbs()} display
    some details about its internal operations.
  }
}

\value{
  Nothing (i.e. an invisible \code{NULL}).
}

\seealso{
  \itemize{
    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.
  }
}

\examples{
## DANGER ZONE!
\dontrun{
  reset_germline_dbs(verbose=TRUE)
}
}

\keyword{utilities}
