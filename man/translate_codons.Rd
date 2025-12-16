\name{translate_codons}

\alias{translate_codons}
\alias{extract_codons}

\title{Extract and translate codons from a set of DNA sequences}

\description{
  The \code{translate_codons()} function extracts and translates codons
  from a set of DNA sequences.
}

\usage{
translate_codons(dna, offset=0, with.init.codon=FALSE)

## Used internally by translate_codons():
extract_codons(dna, offset=0)
}

\arguments{
  \item{dna}{
    A \link[Biostrings]{DNAStringSet} (or \link[Biostrings]{DNAString})
    object containing the codons to translate.
  }
  \item{offset}{
    The number of nucleotides that precede the first codon to translate.
    This must be supplied as a numeric vector with one value per
    sequence in \code{dna}, or as a single value. If the latter,
    then the same offset is used for all sequences.
  }
  \item{with.init.codon}{
    Is the first codon to translate in each DNA sequence the initiation
    codon? By default, \code{with.init.codon} is set to \code{FALSE}, in
    which case \code{translate_codons()} assumes that the first codon
    to translate in each DNA sequence is \emph{not} the initiation codon.

    See the \code{no.init.codon} argument in
    \code{?\link[Biostrings]{translate}} in the \pkg{Biostrings} package
    for more information.
  }
}

\value{
  \code{translate_codons()} returns an \link[Biostrings]{AAStringSet}
  object with one amino acid sequence per input sequence.

  \code{extract_codons()} returns a \link[Biostrings]{DNAStringSet}
  object with one sequence per input sequence. The output sequences
  are obtained by trimming the original sequences as follow:
  \itemize{
    \item On their 5' end, sequences are trimmed by the amount of
          nucleotides specified in \code{offset}.
    \item On their 3' end, sequences are trimmed by the smallest
          amount of nucleotides that makes the length of the trimmed
          sequence a multiple of 3. Note that this will always be 0, 1,
          or 2 nucleotides.
  }
}

\seealso{
  \itemize{
    \item \link[Biostrings]{DNAStringSet} and \link[Biostrings]{AAStringSet}
          objects in the \pkg{Biostrings} package.

    \item The \code{\link[Biostrings]{translate}} function in the
          \pkg{Biostrings} package on which \code{translate_codons()}
          is based.

    \item \code{\link{list_germline_dbs}} to list all the \emph{cached
          germline dbs}, that is, all the germline databases currently
          installed in \pkg{igblastr}'s persistent cache.
  }
}

\examples{
## ---------------------------------------------------------------------
## translate_codons()
## ---------------------------------------------------------------------

## Load germline V gene allele sequences for human:
list_germline_dbs()
db_name <- "_AIRR.human.IGH+IGK+IGL.202410"
V_alleles <- load_germline_db(db_name, region_types="V")
V_alleles  # DNAStringSet object

## Translate them:
V_aa <- translate_codons(V_alleles)
V_aa  # AAStringSet object

## Some human germline V gene allele sequences have a stop codon:
has_stop_codon <- grepl("*", as.character(V_aa), fixed=TRUE)
V_aa[has_stop_codon]

## Handling of initiation codons (see '?translate' in the Biostrings
## package for how initiation codons are handled):
dna2 <- DNAStringSet(c("TTGTCCTTTATA", "GAATCATTTATC", "CTGTCGTTTATT"))
translate_codons(dna2)
translate_codons(dna2, with.init.codon=TRUE)

## ---------------------------------------------------------------------
## extract_codons()
## ---------------------------------------------------------------------
dna <- DNAStringSet(c("CCCAAAGGGTTT",
                      "CCAAAGGGTTT",
                      "CAAAGGGTTT",
                      "AAAGGGTTT"))
extract_codons(dna)
extract_codons(dna, offset=1)
extract_codons(dna, offset=2)
extract_codons(dna, offset=3)
extract_codons(dna, offset=4)

extract_codons(dna, offset=3:0)
extract_codons(dna, offset=4:1)
extract_codons(dna, offset=5:2)
extract_codons(dna, offset=6:3)
}

\keyword{utilities}
