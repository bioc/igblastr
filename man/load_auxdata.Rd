\name{load_auxdata}

\alias{auxdata}

\alias{get_auxdata_path}
\alias{load_auxdata}

\alias{get_igblast_auxiliary_data}
\alias{load_igblast_auxiliary_data}

\title{Access IgBLAST auxiliary data}

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

  The \emph{auxiliary data} informs about the coding frame and CDR3/FWR4
  boundary for each J allele sequence. Although this information is not
  strictly needed by IgBLAST, the latter needs it to compute some of the
  fields in the AIRR-formatted output like \code{vj_in_frame},
  \code{productive}, \code{cdr3}, \code{fwr4}, etc...

  See \code{?\link{load_intdata}} for more information about IgBLAST
  \emph{internal data}.

  \code{get_auxdata_path()} and \code{load_auxdata()} can be used
  to access the \emph{auxiliary data} included in IgBLAST or in
  one of the \emph{cached germline dbs} managed by \pkg{igblastr}.
}

\usage{
get_auxdata_path(igblast_organism, which=c("live", "original"))

load_auxdata(igblast_organism, which=c("live", "original"))
}

\arguments{
  \item{igblast_organism}{
    A single string containing the name of an "IgBLAST organism" as
    returned by \code{\link{list_igblast_organisms}()}.

    Alternatively, this can be the name of a cached germline db.
    Note that this works only for germline dbs that include their
    own \emph{auxiliary data}. You can use:
    \preformatted{    list_germline_dbs(with.auxdata.only=TRUE)}
    to list them.
    See \code{?\link{list_germline_dbs}} for more information.
  }
  \item{which}{
    By default, \code{get_auxdata_path()} and \code{load_auxdata()}
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
  IgBLAST \emph{auxiliary data} is typically included in a standard IgBLAST
  installation. It's located in the \code{optional_file/} directory which
  is itself a subdirectory of IgBLAST \emph{root directory}.

  The data consists of one tabulated file per organism. Each file indicates
  the germline J gene coding frame start position, the J gene type, and the
  CDR3 end position for all known germline J allele sequences. See
  \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html} for
  additional details.

  IMPORTANT NOTE: Using IgBLAST with \emph{auxiliary data} that is not
  compatible with the germline J gene sequences in the germline db can
  cause it to return improper frame status or CDR3 information (other
  returned information will still be correct).
}

\value{
  \code{get_auxdata_path()} returns a single string containing
  the path to the \emph{auxiliary data} included in the IgBLAST installation
  used by \pkg{igblastr}, for the specified organism. Not necessarily
  suitable to use with \code{\link{igblastn}()} (see WARNING below).

  \code{load_auxdata()} returns the \emph{auxiliary data} in a data.frame
  with 1 row per supplied J allele sequence and the same columns as the
  data.frame returned by \code{read_auxdata()}.
  See \code{?\link{read_auxdata}} for more information.
}

\seealso{
  \itemize{
    \item \code{\link{compute_auxdata}} to annotate a set of germline J
          gene allele sequences.

    \item \code{\link{load_intdata}} to access IgBLAST \emph{internal data}.

    \item \link{auxdata_misc} for miscellaneous \emph{auxiliary data} utilities.

    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \code{\link{list_igblast_organisms}} to list the organisms
          supported by IgBLAST (a.k.a. "IgBLAST organisms").

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
          for important information about IgBLAST \emph{auxiliary data}.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## ---------------------------------------------------------------------
## Access IgBLAST auxiliary data for a given "IgBLAST organism"
## ---------------------------------------------------------------------

## IgBLAST only includes auxiliary data for the following organisms:
list_igblast_organisms()

get_auxdata_path("human")

human_auxdata <- load_auxdata("human")
head(human_auxdata)

## ---------------------------------------------------------------------
## Access the auxiliary data included in a germline db
## ---------------------------------------------------------------------

## List germline dbs with auxiliary data:
list_germline_dbs(with.auxdata.only=TRUE)

## Access the auxiliary data included in germline db
## _OGRDB.rhesus_monkey.IGH+IGK+IGL.202602 (note that this data was
## extracted from the AIRR-C JSON files provided by OGRDB for rhesus
## monkey, see '?extract_auxdata_from_ogrdb_json' for more information):
db_name <- "_OGRDB.rhesus_monkey.IGH+IGK+IGL.202602"

get_auxdata_path(db_name)

rhesus_monkey_auxdata <- load_auxdata(db_name)
head(rhesus_monkey_auxdata)
}

\keyword{utilities}
