\name{intdata-IO-utils}

\alias{intdata-IO-utils}
\alias{intdata_IO_utils}

\alias{check_V_ndm_data}
\alias{read_V_ndm_data}
\alias{write_V_ndm_data}

\title{Read/write IgBLAST internal data}

\description{
  Low-level functions to read/write data in \emph{IgBLAST internal data
  format}.
}

\usage{
check_V_ndm_data(V_ndm_data, allow.dup.entries=FALSE)
read_V_ndm_data(filepath)
write_V_ndm_data(V_ndm_data, file="", check.data=FALSE)
}

\arguments{
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

    \item \code{\link{makeogrannote}} to extract the germline V gene allele
          FWR/CDR annotations (a.k.a. V gene delineations) from an
          OGRDB JSON file.

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
