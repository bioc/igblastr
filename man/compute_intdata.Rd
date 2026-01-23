\name{compute_intdata}

\alias{compute_intdata}

\title{Compute IgBLAST internal data}

\description{
  A utility function to annotate a set of germline V gene allele sequences
  in a way similar to how they are annotated in IgBLAST \emph{internal data}.

  Note that the annotation produced by the function can be used by IgBLAST
  as a substitute to the \emph{internal data} shipped with the IgBLAST
  software (and typically included in a standard IgBLAST installation).

  See \code{?\link{load_intdata}} for more information about IgBLAST
  \emph{internal data}.
}

\usage{
compute_intdata(V_alleles)
}

\arguments{
  \item{V_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline V
    gene allele \emph{gapped} sequences.
  }
}

\details{
  COMING SOON...
}

\value{
  Returns the computed \emph{internal data} in a data.frame with 1 row
  per supplied germline V allele sequence and the same columns as the
  data.frame returned by \code{\link{load_intdata}()}.
}

\seealso{
  \itemize{
    \item \code{\link{load_intdata}} to access and manipulate IgBLAST
          \emph{internal data}.

    \item \code{\link{compute_auxdata}} to annotate a set of germline J
          gene allele sequences.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

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
