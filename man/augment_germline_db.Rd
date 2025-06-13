\name{augment_germline_db}

\alias{augment_germline_db}
\alias{augment_germline_db_V}
\alias{augment_germline_db_D}
\alias{augment_germline_db_J}

\title{Add novel gene alleles to a germline db}

\description{
  Three functions to add novel V, D, or J gene alleles to a germline db.
}

\usage{
augment_germline_db_V(db_name, novel_alleles, destdir=".", overwrite=FALSE)
augment_germline_db_D(db_name, novel_alleles, destdir=".", overwrite=FALSE)
augment_germline_db_J(db_name, novel_alleles, destdir=".", overwrite=FALSE)
}

\arguments{
  \item{db_name}{
    A single string that is the name of the cached germline db that
    contains the set of gene alleles to augment.
    Use \code{\link{list_germline_dbs}()} to list the cached germline dbs.

    Depending on the function used (i.e. \code{augment_germline_db_V()},
    \code{augment_germline_db_D()}, or \code{augment_germline_db_J()}),
    the set of alleles to augment is the alleles from the V, D, or J region.
  }
  \item{novel_alleles}{
    A single string that is the path to a FASTA file (possibly gz-compressed)
    where the novel alleles are stored.

    Alternatively, the novel alleles can be supplied as a \emph{named}
    \link[Biostrings]{DNAStringSet} object.
  }
  \item{destdir}{
    A single string that is the path to the "destination directory", that is,
    the directory where the augmented V-, D-, or J-region db is to be created.
    This directory will be created if it doesn't exist already. Note that, by
    default, the augmented region db will be created in the current directory.
  }
  \item{overwrite}{
    If the "destination directory" already contains a V-, D-, or J-region db,
    should it be overwritten?
  }
}

\value{
  An invisible NULL.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

}

\keyword{utilities}
