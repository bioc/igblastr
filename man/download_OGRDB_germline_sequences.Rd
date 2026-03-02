\name{download_OGRDB_germline_sequences}

\alias{download_OGRDB_germline_sequences}

\title{Download germline sequences from OGRDB}

\description{
  The \code{download_OGRDB_germline_sequences()} function downloads
  germline V/D/J gene allele sequences (nucleotides) from OGRDB for
  the specified germline sets.
  The sequences are downloaded in several FASTA files, one per "IMGT group".
}

\usage{
download_OGRDB_germline_sequences(organism="Homo sapiens", germline_sets,
                                  gapped=TRUE,
                                  destdir=".", overwrite=FALSE,
                                  recache=FALSE, ...)
}

\arguments{
  \item{organism}{
    A single string specifying the latin name of the organism for which
    to download the germline sequences.
    See \url{https://ogrdb.airr-community.org/germline_sets} for the
    organisms that are currently available at OGRDB.
  }
  \item{germline_sets}{
    A named integer vector. The names on the vector must be the names of
    the germline sets to download and the values in the vector must be the
    versions of the germline sets to download.
    Go to \url{https://ogrdb.airr-community.org/germline_sets}, select a
    Species, and click on "Show Germline Sets" to see the list of germline
    sets available for that species and their versions.

    Note that each OGRDB germline set corresponds to one or more "IMGT groups",
    where an "IMGT group" is a combination of locus and region type. For
    example:
    \itemize{
      \item Germline set \code{IGH_VDJ} corresponds to IMGT groups IGHV,
            IGHD, and IGHJ.
      \item Germline set \code{IGLambda_VJ} corresponds to IMGT groups IGLV
            and IGLJ.
      \item Germline set \code{C57BL/6 IGH} corresponds to IMGT groups IGHV,
            IGHD, and IGHJ.
      \item Germline set \code{C57BL/6J IGKV} corresponds to IMGT group IGKV.
    }
    For reference, the 10 "IMGT groups" for B-cell Receptors are:
    \preformatted{        IGHV, IGHD, IGHJ, IGHC,
        IGKV, IGKJ, IGKC,
        IGLV, IGLJ, IGLC}

    IMPORTANT: The germline sets in your \code{germline_sets} vector must
    correspond to distinct IMGT groups.

    Here is an example of valid \code{germline_sets} vector for Mus musculus:
    \preformatted{        germline_sets <- c(`C57BL/6 IGH`=5,
        `C57BL/6J IGKV`=1, `IGKJ (all strains)`=1
        `C57BL/6J IGLV`=1, `IGLJ (all strains)`=1)}
  }
  \item{gapped}{
    By default \code{download_OGRDB_germline_sequences()} will download
    the \emph{gapped} V allele sequences.
    Set \code{gapped} to \code{FALSE} to download the \emph{ungapped}
    V allele sequences instead.
  }
  \item{destdir}{
    A single string that is the path to the directory where the FASTA files
    containing the germline sequences should be downloaded.
  }
  \item{overwrite}{
    \code{download_OGRDB_germline_sequences()} won't replace existing
    files, unless \code{overwrite} is set to \code{TRUE}.
  }
  \item{recache}{
    \code{download_OGRDB_germline_sequences()} uses a caching mechanism for
    the data that gets downloaded from OGRDB. The data is cached on disk so
    the caching is persistent across R sessions.

    Set \code{recache} to \code{TRUE} to force a new download (and recaching)
    of the data.
  }
  \item{...}{
    Extra arguments to be passed to the internal call to \code{download.file()}.
    See \code{?\link[utils]{download.file}} in the \pkg{utils} package for
    more information.
  }
}

\value{
  The function returns the list of FASTA files that it produced in an
  invisible character vector that carries the names of the corresponding
  germline sets.
}

\seealso{
  \itemize{
    \item \code{\link{install_custom_germline_db}} to install a germline db
          from user-supplied gene allele sequences.

    \item \code{\link{download_IMGT_germline_sequences}} to download
          germline sequences from IMGT.

    \item The OGRDB website: \url{https://ogrdb.airr-community.org/}.
  }
}

\examples{
## ---------------------------------------------------------------------
## DOWNLOAD RHESUS MONKEY GERMLINE SEQUENCES FROM OGRDB
## ---------------------------------------------------------------------
## Visit https://ogrdb.airr-community.org/germline_sets/Macaca%20mulatta
## for the list of germline sets available for rhesus monkey.

## Download germline gene allele sequences for all BCR loci/regions
## for rhesus monkey to a temporary directory:
germline_sets <- c(IGH_VDJ=2, IGK_VJ=2, IGL_VJ=2)

dir.create(destdir <- tempfile("fasta_dir_"))
download_OGRDB_germline_sequences("Macaca mulatta", germline_sets,
                                  destdir=destdir)

## List the downloaded files:
list.files(destdir)  # 7 FASTA files

## Remove temporary directory:
unlink(destdir, recursive=TRUE)

## ---------------------------------------------------------------------
## DOWNLOAD MOUSE GERMLINE SEQUENCES FROM OGRDB
## ---------------------------------------------------------------------
## Visit https://ogrdb.airr-community.org/germline_sets/Mus%20musculus
## for the list of germline sets available for mouse.

## Download germline gene allele sequences for all BCR loci/regions
## for mouse strain PWD/PhJ to a temporary directory:
germline_sets <- c(`PWD/PhJ IGH`=2,
                   `PWD/PhJ IGKV`=1, `IGKJ (all strains)`=1,
                   `PWD/PhJ IGLV`=1, `IGLJ (all strains)`=1)

dir.create(destdir <- tempfile("fasta_dir_"))
download_OGRDB_germline_sequences("Mus musculus", germline_sets,
                                  destdir=destdir)

## List the downloaded files:
list.files(destdir)  # 7 FASTA files

## Remove temporary directory:
unlink(destdir, recursive=TRUE)
}

\keyword{utilities}
