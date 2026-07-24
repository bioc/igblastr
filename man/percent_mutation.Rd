\name{percent_mutation}

\alias{percent_mutation}

\title{Computes percent mutation in V, D, J segments}

\description{
  Computes the percent mutation in the V, D, J segments of a set of BCR
  or TCR sequences, at the nucleotide or amino acid levels.
}

\usage{
percent_mutation(AIRR_df, for.aa=FALSE, as.matrix=FALSE)
}

\arguments{
  \item{AIRR_df}{
    An AIRR-formatted data.frame or \link[tibble]{tibble} representing
    a set of annotated BCR or TCR sequences. This is typically obtained
    with \code{\link{igblastn}()}.
  }
  \item{for.aa}{
    By default, the percent mutation is computed at the nucleotide level.
    Set \code{for.aa} to \code{TRUE} to compute it at the amino acid level.
  }
  \item{as.matrix}{
    By default, the function returns the results in a data.frame.
    Set \code{as.matrix} to \code{TRUE} to get the results in a matrix.
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
  By default, the function returns a data.frame or \link[tibble]{tibble} with
  1 row per row in \code{AIRR_df}, and with the following columns:
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
  If \code{for.aa} is set to \code{TRUE}, then columns 3 to 5 are named
  \code{v_perc_mut_aa}, \code{d_perc_mut_aa}, and \code{j_perc_mut_aa}.

  If \code{as.matrix} is set to \code{TRUE}, then the function returns
  a 3-column numeric matrix with 1 row per row in \code{AIRR_df}.
  The rownames on the matrix are the sequence ids from
  \code{AIRR_df$sequence_id}.
  The colnames on the matrix are \code{[vdj]_perc_mut} if \code{for.aa}
  is \code{FALSE}, or \code{[vdj]_perc_mut_aa} if it's \code{TRUE}.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item The Rearrangement Schema from the AIRR Standards Working Group
          is documented at
    \url{https://docs.airr-community.org/en/latest/datarep/rearrangements.html}.

    \item \code{\link{extract_sequence_region}} to extract a region of
          interest from a set of analyzed sequences.

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

## Percent mutation in V, D, J segments at the amino acid level:
percent_mutation(AIRR_df, for.aa=TRUE)

## Get the results in a 3-column matrix:
m1 <- percent_mutation(AIRR_df, for.aa=TRUE, as.matrix=TRUE)
m1[1:4, ]
colMeans(m1)

## Percent mutation in V, D, J segments at the nucleotide level:
percent_mutation(AIRR_df)
m2 <- percent_mutation(AIRR_df, as.matrix=TRUE)
m2[1:4, ]
colMeans(m2)

## Note that columns 'v_identity', 'd_identity', and 'j_identity'
## in 'AIRR_df' are equal to '1 - v_perc_mut', '1 - d_perc_mut',
## and '1- j_perc_mut', respectively:
ident <- AIRR_df[ , c("v_identity", "d_identity", "j_identity")]
ident

## Sanity check:
m <- as.matrix(ident)
stopifnot(max((m + m2 - 100)^2, na.rm=TRUE) < 1e-6)
}

\keyword{utilities}
