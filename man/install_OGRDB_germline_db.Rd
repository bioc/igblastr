\name{install_OGRDB_germline_db}

\alias{install_OGRDB_germline_db}
\alias{make_intdata_files_from_ogrdb_jsons}
\alias{make_auxdata_files_from_ogrdb_jsons}

\title{Install a germline db from OGRDB germline sets}

\description{
  COMING SOON...

  If you use OGRDB in your research, please cite:
  \preformatted{  The current landscape of adaptive immune receptor genomic
  and repertoire data: OGRDB and VDJbase
  Lees, Peres, Klein, Amos et al., Nucleic Acids Research, November 2025.
  https://doi.org/10.1093/nar/gkaf1094}
}

\usage{
install_OGRDB_germline_db(organism, germline_sets,
                          suffix="", source_set=FALSE,
                          overwrite=FALSE, verbose=FALSE, ...)
}

\arguments{
  \item{organism}{
    COMING SOON...
  }
  \item{germline_sets}{
    COMING SOON...
  }
  \item{suffix}{
    COMING SOON...
  }
  \item{source_set}{
    COMING SOON...
  }
  \item{overwrite}{
    Set to \code{TRUE} to reinstall if the requested database is already
    installed.
  }
  \item{verbose}{
    Set to \code{TRUE} to have the function display some details about
    its internal operations.
  }
  \item{...}{
    Extra arguments to be passed to the internal call to \code{download.file()}.
    See \code{?\link[utils]{download.file}} in the \pkg{utils} package for
    more information.
  }
}

\details{
  COMING SOON...
}

\value{
  \code{install_OGRDB_germline_db()} returns the name to the newly
  installed germline db as an invisible string.
}

\seealso{
  \itemize{
    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

    \item \code{\link{install_custom_germline_db}} to install a germline db
          from user-supplied gene allele sequences.

    \item \code{\link{download_OGRDB_germline_sequences}} to download
          germline sequences from OGRDB.

    \item \code{\link{download_OGRDB_germline_json}} to download AIRR-C
          JSON files from OGRDB.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \code{\link{use_germline_db}} to select the cached germline db
          to use with \code{igblastn()}.

    \item The OGRDB website: \url{https://ogrdb.airr-community.org/}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## See https://ogrdb.airr-community.org/germline_sets/Mus%20musculus
## for what OGRDB germline sets are available for mouse.

## Install a germline db for mouse substrain BALB/c/ByJ:
germline_sets <- c(`BALB/c IGH`=3,
                   `BALB/c/ByJ IGKV`=1, `BALB/c/ByJ IGLV`=2,
                   `IGKJ (all strains)`=1, `IGLJ (all strains)`=1)

db_name <- install_OGRDB_germline_db("Mus musculus", germline_sets)

list_germline_dbs()

rm_germline_db(db_name)

## MORE EXAMPLES TO COME SOON...
}

\keyword{utilities}
