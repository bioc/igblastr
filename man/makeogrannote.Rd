\name{makeogrannote}

\alias{makeogrannote}
\alias{check_V_ndm_data}
\alias{read_V_ndm_data}
\alias{write_V_ndm_data}

\title{Extract V gene FWR/CDR annotation from OGRDB JSON file}

\description{
  The \code{makeogrannote} function is a reimplementation in R
  of the Python script \code{makeogrannote.py} included in IgBLAST.
}

\usage{
makeogrannote(germline_file)

## Related utilities:
check_V_ndm_data(V_ndm_data, allow.dup.entries=FALSE)
read_V_ndm_data(filepath)
write_V_ndm_data(V_ndm_data, file="", check.data=FALSE)
}

\arguments{
  \item{germline_file}{
    COMING SOON...
  }
  \item{V_ndm_data}{
    COMING SOON...
  }
  \item{allow.dup.entries}{
    COMING SOON...
  }
  \item{filepath}{
    COMING SOON...
  }
  \item{file}{
    COMING SOON...
  }
  \item{check.data}{
    COMING SOON...
  }
}

\details{
  COMING SOON...
}

\value{
  COMING SOON...
}

\seealso{
  \itemize{
    \item \code{\link{compute_imgt_intdata}} to annotate a set of germline V
          gene allele sequences.

    \item \link{intdata_utils} to access and manipulate IgBLAST
          \emph{internal data}.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
## COMING SOON...
}

\keyword{utilities}
