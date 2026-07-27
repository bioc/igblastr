\name{net_charge}

\alias{PKA_VALUES}
\alias{net_charge}
\alias{extract_region_net_charge}

\title{Compute the net charges of amino acid sequences}

\description{
  Compute the theoretical net charges of amino acid sequences using
  the Henderson-Hasselbalch method described here:
  \url{https://www.bachem.com/knowledge-center/peptide-calculator/}
}

\usage{
net_charge(aa_sequences, pH=7.4,
           with.nTer=TRUE, with.cTer=TRUE,
           pKa_values=PKA_VALUES, as.matrix=FALSE)

extract_region_net_charge(AIRR_df, region_type,
                          pH=7.4, pKa_values=PKA_VALUES,
                          as.matrix=FALSE)
}

\arguments{
  \item{aa_sequences}{
    A character vector or \link[Biostrings]{AAStringSet} object containing
    amino acid sequences.
  }
  \item{pH}{
    A pH value. Note that this is set to 7.4 (physiological pH) by default,
    rather than 7.0 (neutral pH).
  }
  \item{with.nTer, with.cTer}{
    Set to \code{FALSE} to exclude the terminal contributions to the net
    charge. Do this only if the supplied amino acid sequences are embedded
    inside larger proteins rather than isolated peptides.
  }
  \item{pKa_values}{
    The vectors of pKa values to use for the computation.
    This must be a numeric vector with the following names on it:
    nTer, R, H, K, D, E, C, Y, cTer.

    By default, \code{pKa_values} is set to \code{PKA_VALUES}, which
    contains the values listed in the "pKa Values of Ionizable Groups"
    section of the following document:
    \url{https://www.chemverify.com/learn/amino-acid-reference-table}

    Note that the \pkg{Peptides} package contains several predefined
    vectors of pKa values that can be passed to the \code{pKa_values}
    argument of \code{net_charge()} or \code{extract_region_net_charge()}.
    See "USE PREDEFINED pKa VALUES FROM THE Peptides PACKAGE" section in
    the Examples below for some concrete examples.
  }
  \item{as.matrix}{
    Set to \code{TRUE} to return a matrix with 1 row per amino acid
    sequence and 1 column per contributing group. Note that the names
    of the contributing groups are the same as the names on the
    \code{PKA_VALUES} vector.
  }
  \item{AIRR_df}{
    An AIRR-formatted data.frame or \link[tibble]{tibble} representing
    a set of annotated BCR or TCR sequences. This is typically obtained
    with \code{\link{igblastn}()}.
  }
  \item{region_type}{
    The type of region (V, D, J, FWR1-4, or CDR1-3) to extract from
    \code{AIRR_df}.
  }
}

\value{
  \code{net_charge()} returns a numeric vector with one element per
  sequence in \code{aa_sequences}.

  \code{extract_region_net_charge()} returns a numeric vector with
  one element per row in \code{AIRR_df}.
  Note that the function simply does:
  \preformatted{  aa_sequences <- extract_sequence_region(AIRR_df, region_type,
                                          as.aa=TRUE)
  net_charge(aa_sequences, with.nTer=FALSE, with.cTer=FALSE, ...)}
}

\seealso{
  \itemize{
    \item \code{\link[Peptides]{charge}} in the \pkg{Peptides} package for
          a function similar to \code{net_charge()}.

    \item \code{\link{extract_sequence_region}} to extract a region of
          interest from a set of analyzed sequences.

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

## ---------------------------------------------------------------------
## BASIC EXAMPLE
## ---------------------------------------------------------------------

query <- system.file(package="igblastr", "extdata",
                     "BCR", "heavy_sequences.fasta")
use_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")

AIRR_df <- igblastn(query)

cdr3_aa <- extract_sequence_region(AIRR_df, region_type="cdr3", as.aa=TRUE)
cdr3_charges <- net_charge(cdr3_aa, with.nTer=FALSE, with.cTer=FALSE)
head(cdr3_charges)
summary(cdr3_charges)

## extract_region_net_charge() conveniently wraps the two steps above into
## a single one:
cdr3_charges2 <- extract_region_net_charge(AIRR_df, region_type="cdr3")

stopifnot(identical(cdr3_charges, cdr3_charges2))  # sanity check

## ---------------------------------------------------------------------
## USE PREDEFINED pKa VALUES FROM THE Peptides PACKAGE
## ---------------------------------------------------------------------

## By default net_charge() uses the following pKa values (taken from
## https://www.chemverify.com/learn/amino-acid-reference-table):
PKA_VALUES

## Note that the Peptides package contains several predefined
## vectors of pKa values that can be used instead:
library(Peptides)
data(AAdata, package="Peptides")
do.call(rbind, AAdata$pK)

## To use one of them, simply pass it to the 'pKa_values' argument
## of the net_charge() or extract_region_net_charge() function:
pKa_values <- AAdata$pK$Lehninger
pKa_values

cdr3_charges3 <- extract_region_net_charge(AIRR_df, region_type="cdr3",
                                           pKa_values=pKa_values)
head(cdr3_charges3)
summary(cdr3_charges3)

## ---------------------------------------------------------------------
## igblastr::net_charge() VS Peptides::charge()
## ---------------------------------------------------------------------

## Here are some important differences to keep in mind when comparing
## results from the two functions:
##   (a) The pH value used by default differs: it's 7.4 in net_charge(),
##       and 7 in Peptides::charge().
##   (b) The vectors of pKa values used by default differs: it's
##       PKA_VALUES in net_charge(), and the "Lehninger pKa values" in
##       Peptides::charge().
##   (c) Peptides::charge() **always** counts the terminal contributions
##       to the net charge with no option to exclude them.

fwr1_aa <- extract_sequence_region(AIRR_df, region_type="fwr1", as.aa=TRUE)

pKa_values <- AAdata$pK$Lehninger  # Lehninger pKa values

res1 <- net_charge(fwr1_aa, pH=7, pKa_values=pKa_values)
res2 <- charge(fwr1_aa, pH=7, pKscale="Lehninger")
stopifnot(all.equal(res1, res2, tolerance=1e-3))

res1 <- net_charge(fwr1_aa, pH=7.4, pKa_values=pKa_values)
res2 <- charge(fwr1_aa, pH=7.4, pKscale="Lehninger")
stopifnot(all.equal(res1, res2, tolerance=1e-3))

pKa_values <- AAdata$pK$EMBOSS  # EMBOSS pKa values

res1 <- net_charge(fwr1_aa, pH=7, pKa_values=pKa_values)
res2 <- charge(fwr1_aa, pH=7, pKscale="EMBOSS")
stopifnot(all.equal(res1, res2, tolerance=1e-3))

res1 <- net_charge(fwr1_aa, pH=7.4, pKa_values=pKa_values)
res2 <- charge(fwr1_aa, pH=7.4, pKscale="EMBOSS")
stopifnot(all.equal(res1, res2, tolerance=1e-3))
}

\keyword{utilities}
