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

\details{
  The percent mutation in the V segment is computed as follow.

  For each sequence in \code{AIRR_df}:
  \enumerate{
    \item First the number of mutations is obtained by computing
          the Hamming distance d between the string reported in the
          \code{v_sequence_alignment} column and the string reported
          in the \code{v_germline_alignment} column. Note that the
          two strings are guaranteed to \emph{always} have the same
          length L.
    \item The percent mutation is 100 * d / L.
  }
  If \code{for.aa} is set to \code{TRUE}, then replace
  \code{v_sequence_alignment} and \code{v_germline_alignment} with
  \code{v_sequence_alignment_aa} and \code{v_germline_alignment_aa} above.

  The percent mutation in the D and J segments is computed in the same
  manner using the corresponding \code{d_*} and \code{j_*} columns.
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
