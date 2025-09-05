\name{igbrowser}

\alias{igbrowser}
\alias{bcr_browser}

\title{Display annotated BCR sequences in the browser}

\description{
  Use \code{igbrowser()} to display the annotated BCR sequences
  returned by \code{\link{igblastn}()} in the browser.
  For each sequence, the V, D, and J segments are shown as well as
  the FWR1-4 and CDR1-3 regions.
  Additionally, the C segments are shown if the C-region information
  is available.
}

\usage{
igbrowser(AIRR_df, show.full.sequence=FALSE, dna.coloring=TRUE,
          Vcolor="#FFDDD2", Dcolor="#CFC", Jcolor="#CEF", Ccolor="#EEC",
          FWRcolor="#C9D", CDRcolor="#EE4")
}

\arguments{
  \item{AIRR_df}{
    The AIRR-formatted data.frame or \link[tibble]{tibble} returned by
    \code{\link{igblastn}()}.
    Note that calling \code{igbrowser()} on a data.frame with thousands
    of rows is quite resource-intensive (it can even crash your browser!),
    so in this case we recommend subsetting the data.frame before passing
    it to \code{igbrowser()} to keep the number of rows under 2000.
  }
  \item{show.full.sequence}{
    By default, the part of the BCR sequences upstream of the V region
    is not shown. Set \code{show.full.sequence} to \code{TRUE} to show it.
  }
  \item{dna.coloring}{
    Whether the nucleotides in the BCR sequences (\code{sequence} column
    in \code{AIRR_df}) should be colored or not.
  }
  \item{Vcolor, Dcolor, Jcolor, Ccolor}{
    The background colors of the V, D, J, and C segments of the BCR
    sequences. Note that the C segments are shown only if \code{AIRR_df}
    contains C-region information.
  }
  \item{FWRcolor, CDRcolor}{
    The background colors of the Framework Regions (FWR1-4) and
    Complementarity-Determining Regions (CDR1-3), respectively.
  }
}

\value{
  0 or the error code returned by the internal call to
  \code{\link[utils]{browseURL}()}, invisibly.
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
                     "heavy_sequences.fasta")
use_germline_db("_AIRR.human.IGH+IGK+IGL.202501")

## ---------------------------------------------------------------------
## With C regions
## ---------------------------------------------------------------------

use_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")
AIRR_df <- igblastn(query)
igbrowser(AIRR_df)

## By default, the part of the sequences upstream of the V region is
## not shown. Use 'show.full.sequence=TRUE' to show the full sequences:
igbrowser(AIRR_df, show.full.sequence=TRUE)

## ---------------------------------------------------------------------
## No C regions
## ---------------------------------------------------------------------

use_c_region_db("")
AIRR_df2 <- igblastn(query)
igbrowser(AIRR_df2)
}

\keyword{utilities}
