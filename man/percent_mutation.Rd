\name{percent_mutation}

\alias{percent_mutation}

\title{Computes percent mutation in V, D, J segments}

\description{
  Computes the percent mutation in the V, D, J segments of a set of BCR
  or TCR sequences, at the nucleotide or amino acid levels.
}

\usage{
percent_mutation(AIRR_df, for.aa=FALSE)
}

\arguments{
  \item{AIRR_df}{
    The AIRR-formatted data.frame or \link[tibble]{tibble} returned by
    \code{\link{igblastn}()}.
  }
  \item{for.aa}{
    By default, the percent mutation is computed at the nucleotide level.
    Set \code{for.aa} to \code{TRUE} to compute it at the amino acid level.
  }
}

\value{
  A data.frame or \link[tibble]{tibble} with 1 row per row in
  \code{AIRR_df}, and with the following columns:
  \enumerate{
    \item \code{sequence_id}: the \code{sequence_id} column taken as-is from
          \code{AIRR_df};
    \item \code{locus}: the \code{locus} column taken as-is from
          \code{AIRR_df};
    \item \code{v_perc_mut}: the percent mutation in the V segment of
          the BCR or TCR sequences;
    \item \code{d_perc_mut}: the percent mutation in the D segment of
          the BCR or TCR sequences;
    \item \code{j_perc_mut}: the percent mutation in the J segment of
          the BCR or TCR sequences.
  }
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.

    \item \link[tibble]{tibble} objects implemented in the
          \pkg{tibble} package.
  }
}

\examples{
if (!has_igblast()) install_igblast()

query <- system.file(package="igblastr", "extdata",
                     "BCR", "heavy_sequences.fasta")
use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")

AIRR_df <- igblastn(query)

percent_mutation(AIRR_df)

percent_mutation(AIRR_df, for.aa=TRUE)
}

\keyword{utilities}
