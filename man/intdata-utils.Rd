\name{intdata-utils}

\alias{intdata-utils}
\alias{intdata_utils}
\alias{intdata}

\alias{get_intdata_path}
\alias{load_intdata}
\alias{translate_V_alleles}
\alias{V_allele_has_stop_codon}
\alias{annotate_heavy_V_alleles}
\alias{annotate_light_V_alleles}

\title{Access and manipulate IgBLAST internal data}

\description{
  IgBLAST \emph{internal data} is expected to annotate all the known
  germline V gene alleles for a given organism. It is provided by NCBI
  and is typically included in a standard IgBLAST installation.

  The \pkg{igblastr} package provides a small set of utilities to access
  and manipulate IgBLAST internal data.
}

\usage{
get_intdata_path(organism, for.aa=FALSE, domain_system=c("imgt", "kabat"),
                 which=c("live", "original"))

load_intdata(organism, for.aa=FALSE, domain_system=c("imgt", "kabat"),
             which=c("live", "original"))

translate_V_alleles(V_alleles, intdata, region=NULL)

V_allele_has_stop_codon(V_alleles, intdata)
}

\arguments{
  \item{organism}{
    A single string containing the name of an organism as
    returned by \code{\link{list_igblast_organisms}()}.
  }
  \item{for.aa}{
    By default, the data.frame returned by \code{load_intdata()} contains
    FWR/CDR start/end positions with respect to the nucleotide sequences
    of the germline V alleles.
    Setting \code{for.aa} to \code{TRUE} will return a data.frame with
    positions that are with respect to the amino acid sequences of the
    germline V alleles.
  }
  \item{domain_system}{
    Domain system to be used for segment annotation. Must be \code{"imgt"}
    (the default) or \code{"kabat"}.
  }
  \item{which}{
    By default, \code{get_intdata_path()} and \code{load_intdata()}
    access the "live IgBLAST data", that is, the IgBLAST data that
    the user has possibly updated with \code{update_live_igdata()}.
    Depending on whether updates were applied or not, the "live IgBLAST data"
    might differ from the original IgBLAST data.

    Set \code{which} to \code{"original"} if you want to access
    the original IgBLAST data instead.

    See \code{?\link{update_live_igdata}} for more information about
    "live" and "original" IgBLAST data.
  }
  \item{V_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline V
    gene allele sequences.
  }
  \item{intdata}{
    A data.frame as returned by
    \code{load_intdata(., for.aa=FALSE, domain_system="imgt")}.
  }
  \item{region}{
    The region to translate. This can be set to \code{"fwr1"}, \code{"cdr1"},
    \code{"fwr2"}, \code{"cdr2"}, or \code{"fwr3"}.
    By default (i.e. when \code{region} is omitted or set to \code{NULL}),
    the entire coding frame in each allele sequence is translated.
  }
}

\details{
  IgBLAST \emph{internal data} is typically included in a standard IgBLAST
  installation. It's located in the \code{internal_data/} directory which
  is itself a subdirectory of IgBLAST \emph{root directory}.
}

\value{
  \code{get_intdata_path()} returns a single string containing
  the path to the internal data included in the IgBLAST installation used
  by \pkg{igblastr}, for the specified organism.

  \code{load_intdata()} returns the internal data in a data.frame
  with 1 row per germline V allele sequence and the following columns:
  \enumerate{
    \item \code{allele_name}: allele name;
    \item \code{fwr1_start}, \code{fwr1_end}: FWR1 start/end
          positions (1-based);
    \item \code{cdr1_start}, \code{cdr1_end}: CDR1 start/end
          positions (1-based);
    \item \code{fwr2_start}, \code{fwr2_end}: FWR2 start/end
          positions (1-based);
    \item \code{cdr2_start}, \code{cdr2_end}: CDR2 start/end
          positions (1-based);
    \item \code{fwr3_start}, \code{fwr3_end}: FWR3 start/end
          positions (1-based);
    \item \code{chain_type}: chain type;
    \item \code{coding_frame_start}: first coding frame start
          position (0-based).
  }

  \code{translate_V_alleles()} returns a named character vector with 1
  amino acid sequence per supplied allele. The vector contains an \code{NA}
  for any allele that is not annotated in \code{intdata} or for which
  the required information is \code{NA}. The names on it are
  the names of the supplied alleles.

  \code{V_allele_has_stop_codon()} returns a named logical vector with 1
  value per supplied allele. The vector contains an \code{NA} for any
  allele that is not annotated in \code{intdata} or for which
  \code{intdata$coding_frame_start} has an \code{NA}. The names on it are
  the names of the supplied alleles.
}

\seealso{
  \itemize{
    \item \link{auxdata_utils} to access, manipulate, and generate IgBLAST
          auxiliary data.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link{translate_codons}} function on which
          \code{translate_V_alleles()} is based.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

## ---------------------------------------------------------------------
## list_igblast_organisms() and get_intdata_path()
## ---------------------------------------------------------------------

list_igblast_organisms()

get_intdata_path("rabbit")
rabbit_intdata <- load_intdata("rabbit")
head(rabbit_intdata)

rabbit_intdata2 <- load_intdata("rabbit", for.aa=TRUE)
head(rabbit_intdata2)

## The values in the "end" cols in 'rabbit_intdata' are exactly 3 times
## those in the "end" cols in 'rabbit_intdata2':
end_colnames <- grep("_end$", colnames(rabbit_intdata), value=TRUE)
stopifnot(identical(rabbit_intdata [ , end_colnames],
                    rabbit_intdata2[ , end_colnames] * 3L))

## ---------------------------------------------------------------------
## translate_V_alleles() and V_allele_has_stop_codon()
## ---------------------------------------------------------------------

human_intdata <- load_intdata("human")

db_name <- "_AIRR.human.IGH+IGK+IGL.202410"
V_alleles <- load_germline_db(db_name, region_types="V")
V_alleles  # DNAStringSet object

V_aa <- translate_V_alleles(V_alleles, human_intdata)
head(V_aa)

fwr1 <- translate_V_alleles(V_alleles, human_intdata, region="fwr1")
head(fwr1)

cdr1 <- translate_V_alleles(V_alleles, human_intdata, region="cdr1")
head(cdr1)

fwr2 <- translate_V_alleles(V_alleles, human_intdata, region="fwr2")
head(fwr2)

## No sequence in 'V_aa' should contain the letter "*" which is used
## by translate_V_alleles() to represent a stop codon. However, seven
## V alleles in _AIRR.human.IGH+IGK+IGL.202410 seem to disobey:
has_stop_codon <- grepl("*", V_aa, fixed=TRUE)
V_aa[has_stop_codon]
V_alleles[has_stop_codon]
}

\keyword{utilities}
