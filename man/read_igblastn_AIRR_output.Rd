\name{read_igblastn_AIRR_output}

\alias{read_igblastn_AIRR_output}

\title{igblastn output format 19 (AIRR format)}

\description{
  Read \code{igblastn} output format 19 (AIRR format). This is the
  output produced by \code{igblastn} when \code{outfmt} is set
  to \code{"AIRR"} or 19.

  This format is sometimes called "Rearrangement summary report" or
  simply "AIRR rearrangement tabular" format. See for example IgBLAST
  web interface at \url{https://www.ncbi.nlm.nih.gov/igblast/}.
}

\usage{
read_igblastn_AIRR_output(out)
}

\arguments{
  \item{out}{
    The path to a file containing the output produced by \code{igblastn}
    when \code{outfmt} is set to \code{"AIRR"} or 19.
  }
}

\value{
  A data.frame with 1 row per query sequence and many columns.

  The columns are standard "Rearrangement Schema" fields, which
  are defined and documented by the AIRR Community at
  \url{https://docs.airr-community.org/en/latest/datarep/rearrangements.html#fields}
}

\seealso{
  \itemize{
    \item The "Rearrangement Schema" at
          \url{https://docs.airr-community.org/en/latest/datarep/rearrangements.html}

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{read_igblastn_fmt7_output}} to read and parse
          \code{igblastn} output format 7.

    \item IgBLAST web interface at \url{https://www.ncbi.nlm.nih.gov/igblast/}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## ---------------------------------------------------------------------
## Access query sequences and select germline and C-region dbs to use
## ---------------------------------------------------------------------

## Files 'heavy_sequences.fasta' and 'light_sequences.fasta' included
## in igblastr contain 250 paired heavy- and light- chain sequences (125
## sequences in each file) downloaded from OAS (the Observed Antibody
## Space database):
filenames <- paste0(c("heavy", "light"), "_sequences.fasta")
query <- system.file(package="igblastr", "extdata", "BCR", filenames)

## Keep only the first 10 sequences from each file:
query <- c(head(readDNAStringSet(query[[1L]]), n=10),
           head(readDNAStringSet(query[[2L]]), n=10))

## Select the germline and C-region dbs to use with igblastn():
use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")
use_c_region_db("_IMGT.human.IGH+IGK+IGL.202605")

## ---------------------------------------------------------------------
## Call igblastn()
## ---------------------------------------------------------------------

out <- tempfile()
igblastn(query, out=out)
AIRR_df <- read_igblastn_AIRR_output(out)

class(AIRR_df)
dim(AIRR_df)  # 1 row per query sequence

tibble(AIRR_df)
}

\keyword{manip}
\keyword{utilities}
