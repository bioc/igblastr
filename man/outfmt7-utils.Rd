\name{outfmt7-utils}

\alias{outfmt7-utils}

\alias{list_outfmt7_specifiers}
\alias{print.outfmt7_specifiers}

\alias{qseqid}
\alias{qseqid.query_details}
\alias{summary.query_details}
\alias{print.query_details}
\alias{print.VDJ_rearrangement_summary}
\alias{print.VDJ_junction_details}
\alias{print.subregion_sequence_details}
\alias{print.alignment_summary}
\alias{print.hit_table}
\alias{qseqid.fmt7record}
\alias{print.fmt7record}
\alias{print.fmt7footer}
\alias{parse_outfmt7}

\title{Handle igblastn output format 7}

\description{
  Some utilities to handle igblastn output format 7.
}

\usage{
list_outfmt7_specifiers()

parse_outfmt7(out_lines)
}

\arguments{
  \item{out_lines}{
    The character vector returned by
    \code{igblatsn(query, outfmt=7, parse.out=FALSE, ...)}.
  }
}

\value{
  \code{list_outfmt7_specifiers()} returns the list of format specifiers
  supported by \code{igblastn()} formatting option 7.

  \code{parse_outfmt7(out_lines)} returns the parsed form of \code{out_lines}
  in a list.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

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
use_germline_db("_AIRR.human.IGH+IGK+IGL.202501")
use_c_region_db("_IMGT.human.IGH+IGK+IGL.202412")

## ---------------------------------------------------------------------
## FIRST igblastn RUN: GET OUTPUT IN FORMAT 7
## ---------------------------------------------------------------------

## For this first run we specify 'outfmt=7' and 'parse.out=FALSE':
out_lines <- igblastn(query, outfmt=7, parse.out=FALSE)
out_lines  # raw output

out <- parse_outfmt7(out_lines)  # parse the output

## Output contains one record per query sequence:
length(out$records)  # 20

## Each record can have 5 or 6 sections:
##   1. query_details
##   2. VDJ_rearrangement_summary
##   3. VDJ_junction_details
##   4. subregion_sequence_details (can be missing)
##   5. alignment_summary
##   6. hit_table

## Taking a close look at the first record:
rec1 <- out$records[[1]]
rec1

qseqid(rec1)    # query sequence id associated with this record

rec1$hit_table  # data.frame with the standard columns

## ---------------------------------------------------------------------
## SECOND igblastn RUN: GET OUTPUT IN CUSTOMIZED FORMAT 7
## ---------------------------------------------------------------------

## For this second run we request a customized format 7 by supplying
## space delimited format specifiers:
outfmt <- "7 qseqid sseqid pident nident length score"
out <- igblastn(query, outfmt=outfmt)

## Taking a close look at the first record:
rec1 <- out$records[[1]]
rec1$hit_table  # data.frame with the requested columns (+ the
                # automatic "chaintype" column)
}

\keyword{manip}
