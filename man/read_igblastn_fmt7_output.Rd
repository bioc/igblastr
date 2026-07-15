\name{read_igblastn_fmt7_output}

\alias{outfmt7-utils}
\alias{outfmt7_utils}
\alias{fmt7-utils}
\alias{fmt7_utils}

\alias{read_igblastn_fmt7_output}

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

\alias{list_outfmt7_specifiers}
\alias{print.outfmt7_specifiers}

\title{igblastn output format 7}

\description{
  Read and parse \code{igblastn} output format 7. This is the output
  produced by \code{igblastn} when \code{outfmt} is set to 7.

  This format is sometimes called "Tabular with comment lines" or
  simply "Tabular" format. See for example IgBLAST web interface at
  \url{https://www.ncbi.nlm.nih.gov/igblast/}.
}

\usage{
read_igblastn_fmt7_output(out)

## Related utilities:
parse_outfmt7(out_lines)
list_outfmt7_specifiers()
}

\arguments{
  \item{out}{
    The path to a file containing the output produced by \code{igblastn}
    when \code{outfmt} is set to 7.
  }
  \item{out_lines}{
    The character vector returned by
    \code{igblatsn(query, outfmt=7, parse.out=FALSE, ...)}.
  }
}

\value{
  \code{read_igblastn_fmt7_output(out)} returns
  \code{parse_outfmt7(readLines(out))}.

  \code{parse_outfmt7(out_lines)} returns the parsed form of \code{out_lines}
  in a list.

  \code{list_outfmt7_specifiers()} returns the list of format specifiers
  supported by \code{igblastn} formatting option 7.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{read_igblastn_AIRR_output}} to read \code{igblastn}
          output format 19 (AIRR format).

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
fasta_files <- paste0(c("heavy", "light"), "_sequences.fasta")
query <- system.file(package="igblastr", "extdata", "BCR", fasta_files)

## Keep only the first 10 sequences from each file:
query <- c(head(readDNAStringSet(query[[1L]]), n=10),
           head(readDNAStringSet(query[[2L]]), n=10))

## Select the germline and C-region dbs to use with igblastn():
use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")
use_c_region_db("_IMGT.human.IGH+IGK+IGL.202605")

## ---------------------------------------------------------------------
## FIRST igblastn RUN: GET OUTPUT IN FORMAT 7
## ---------------------------------------------------------------------

parsed_out7 <- igblastn(query, outfmt=7)

## Note that the above is equivalent to:
out <- tempfile()
igblastn(query, outfmt=7, out=out)
parsed_out7b <- read_igblastn_fmt7_output(out)
stopifnot(identical(parsed_out7b, parsed_out7))

## and to:
out_lines <- igblastn(query, outfmt=7, parse.out=FALSE)
out_lines  # raw output
parsed_out7c <- parse_outfmt7(out_lines)
stopifnot(identical(parsed_out7c, parsed_out7))

## Now taking a closer look at the output...

## Output contains one record per query sequence:
length(parsed_out7$records)  # 20

## Each record can have 5 or 6 sections:
##   1. query_details
##   2. VDJ_rearrangement_summary
##   3. VDJ_junction_details
##   4. subregion_sequence_details (can be missing)
##   5. alignment_summary
##   6. hit_table

## Taking a close look at the first record:
rec1 <- parsed_out7$records[[1]]
rec1

qseqid(rec1)    # query sequence id associated with this record

rec1$hit_table  # data.frame with the standard columns

## ---------------------------------------------------------------------
## SECOND igblastn RUN: GET OUTPUT IN CUSTOMIZED FORMAT 7
## ---------------------------------------------------------------------

## For this second run we request a customized format 7 by supplying
## space delimited format specifiers. Use list_outfmt7_specifiers() to
## get the list of format specifiers supported by igblastn formatting
## option 7:
list_outfmt7_specifiers()
outfmt <- "7 qseqid sseqid pident nident length score"
parsed_out7 <- igblastn(query, outfmt=outfmt)

## Taking a close look at the first record:
rec1 <- parsed_out7$records[[1]]
rec1$hit_table  # data.frame with the requested columns (+ the
                # automatic "chaintype" column)
}

\keyword{manip}
\keyword{utilities}
