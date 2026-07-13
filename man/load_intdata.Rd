\name{load_intdata}

\alias{intdata}

\alias{get_intdata_path}
\alias{load_intdata}

\title{Access IgBLAST internal data}

\description{
  IgBLAST provides annotations for the germline V and J gene alleles
  of human, mouse, rabbit, rat, and rhesus monkey. Note that these are
  the organisms returned by \code{\link{list_igblast_organisms}()}.
  These annotations are provided by NCBI and are typically included in
  a standard IgBLAST installation.

  Two types of annotations are provided for each organism:
  \enumerate{
    \item The \emph{internal data}: this data is expected to annotate all
          the known germline V gene alleles for a given organism.
    \item The \emph{auxiliary data}: this data is expected to annotate all
          the known germline J gene alleles for a given organism.
  }

  The \emph{internal data} informs about the FWR/CDR boundaries on
  the V allele sequences, also known as the V gene delineations.
  This information is needed and used internally by IgBLAST.

  See \code{?\link{load_auxdata}} for more information about IgBLAST
  \emph{auxiliary data}.

  \code{get_intdata_path()} and \code{load_intdata()} can be used
  to access the \emph{internal data} included in IgBLAST or in
  one of the \emph{cached germline dbs} managed by \pkg{igblastr}.
}

\usage{
get_intdata_path(igblast_organism, for.aa=FALSE,
                 domain_system=c("imgt", "kabat"),
                 which=c("live", "original"))

load_intdata(igblast_organism, for.aa=FALSE,
             domain_system=c("imgt", "kabat"),
             which=c("live", "original"))
}

\arguments{
  \item{igblast_organism}{
    A single string containing the name of an "IgBLAST organism" as
    returned by \code{\link{list_igblast_organisms}()}.

    Alternatively, this can be the name of a cached germline db.
    Note that this works only for germline dbs that include their
    own \emph{internal data}. You can use:
    \preformatted{    list_germline_dbs(with.intdata.only=TRUE)}
    to list them.
    See \code{?\link{list_germline_dbs}} for more information.
  }
  \item{for.aa}{
    By default, the data.frame returned by \code{load_intdata()} contains
    FWR/CDR boundaries reported with respect to the nucleotide sequences
    of the germline V alleles.
    Setting \code{for.aa} to \code{TRUE} will return a data.frame where
    they are reported with respect to the amino acid sequences of the
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
}

\details{
  IgBLAST \emph{internal data} is typically included in a standard IgBLAST
  installation. It's located in the \code{internal_data/} directory which
  is itself a subdirectory of IgBLAST \emph{root directory}.
}

\value{
  \code{get_intdata_path()} returns a single string containing
  the path to the \emph{internal data} included in the IgBLAST
  installation used by \pkg{igblastr}, for the specified organism.

  \code{load_intdata()} returns the \emph{internal data} in a data.frame
  with 1 row per germline V gene allele sequence and the same columns
  as the data.frame returned by \code{read_ndm_data()}.
  See \code{?\link{read_ndm_data}} for more information.
}

\seealso{
  \itemize{
    \item \code{\link{compute_V_gene_delineations}} to annotate a set of
          germline V gene allele sequences.

    \item \code{\link{load_auxdata}} to access IgBLAST \emph{auxiliary data}.

    \item \link{intdata_misc} for miscellaneous \emph{internal data} utilities.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \code{\link{list_igblast_organisms}} to list the organisms
          supported by IgBLAST (a.k.a. "IgBLAST organisms").

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

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
## Access IgBLAST internal data for a given "IgBLAST organism"
## ---------------------------------------------------------------------

## IgBLAST only includes internal data for the following organisms:
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
## Access the internal data included in a germline db
## ---------------------------------------------------------------------

## List germline dbs with internal data:
list_germline_dbs(with.intdata.only=TRUE)

## Access the internal data included in germline db
## _OGRDB.human.IGH+IGK+IGL.202605 (this data was inferred from
## the gaps in the germline V gene allele sequences provided by
## AIRR-community/OGRDB):
db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"

get_intdata_path(db_name)

human_intdata <- load_intdata(db_name)
head(human_intdata)
}

\keyword{utilities}
