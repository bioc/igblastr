\name{extract_sequence_region}

\alias{extract_sequence_region}
\alias{extract_region_length}

\title{Extract region of interest from a set of analyzed sequences}

\description{
  Extract the V, D, J, FWR1-4, or CDR1-3 region of each sequence in a
  set of analyzed BCR or TCR sequences.

  \code{extract_sequence_region()} extract the region sequences (as
  subsequences of the original BCR or TCR sequences), while
  \code{extract_region_length()} extract their lengths only.
}

\usage{
extract_sequence_region(AIRR_df, region_type, as.aa=FALSE)

extract_region_length(AIRR_df, region_type)
}

\arguments{
  \item{AIRR_df}{
    An AIRR-formatted data.frame or \link[tibble]{tibble} representing
    a set of annotated BCR or TCR sequences. This is typically obtained
    with \code{\link{igblastn}()}.
  }
  \item{region_type}{
    The type of region (V, D, J, FWR1-4, or CDR1-3) to extract from
    \code{AIRR_df}.
  }
  \item{as.aa}{
    Whether the regions should be extracted as nucleotide sequences (the
    default) or amino acid sequences.
  }
}

\value{
  \code{extract_sequence_region()} returns a character vector (when
  \code{as.aa} is \code{FALSE}) or an \link[Biostrings]{AAStringSet}
  object (when \code{as.aa} is \code{TRUE}) with one element per row
  in \code{AIRR_df}.

  Note that:
  \itemize{
    \item If \code{region_type} is \code{"D"} or \code{"C"}, then the
          character vector returned when \code{as.aa} is \code{FALSE}
          can contain \code{NA}s. These correspond to analyzed sequences
          for which the requested region is not defined or could not
          be identified by IgBLAST. This is the case for example
          if \code{region_type} is set to \code{"D"} but the analyzed
          sequence is from the light chain, or if \code{region_type} is
          set to \code{"C"} but the analyzed sequence is not long enough
          to reach into the C region.
    \item For the same reason, if \code{region_type} is \code{"D"} or
          \code{"C"}, then the \link[Biostrings]{AAStringSet} object
          returned when \code{as.aa} is \code{TRUE} can contain empty
          sequences. In this case, empty sequences are used instead of
          \code{NA}s, with the same meaning as above.
  }

  \code{extract_region_length()} returns an integer vector with
  one element per row in \code{AIRR_df}.
  Note that it just does \code{nchar(extract_sequence_region(...))}.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item The Rearrangement Schema from the AIRR Standards Working Group
          is documented at
    \url{https://docs.airr-community.org/en/latest/datarep/rearrangements.html}.

    \item \code{\link{percent_mutation}} to compute the percent mutation
          in the V, D, J segments of a set of BCR or TCR sequences.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.

    \item \link[tibble]{tibble} objects implemented in the
          \pkg{tibble} package.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## ---------------------------------------------------------------------
## extract_sequence_region()
## ---------------------------------------------------------------------

query <- system.file(package="igblastr", "extdata",
                     "BCR", "heavy_sequences.fasta")
use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")

AIRR_df <- igblastn(query)

## Regions as nucleotide sequences:

V_regions <- extract_sequence_region(AIRR_df, region_type="V")
DNAStringSet(V_regions)

cdr3 <- extract_sequence_region(AIRR_df, region_type="cdr3")
DNAStringSet(cdr3)

## Regions as amino acid sequences:

extract_sequence_region(AIRR_df, region_type="cdr3", as.aa=TRUE)

extract_sequence_region(AIRR_df, region_type="fwr4", as.aa=TRUE)

## ---------------------------------------------------------------------
## extract_region_length()
## ---------------------------------------------------------------------

cdr3_lens <- extract_region_length(AIRR_df, region_type="cdr3")
summary(cdr3_lens)

D_region_lens <- extract_region_length(AIRR_df, region_type="D")
summary(D_region_lens)
}

\keyword{utilities}
