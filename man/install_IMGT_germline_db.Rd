\name{install_IMGT_germline_db}

\alias{install_IMGT_germline_db}
\alias{list_IMGT_releases}
\alias{list_IMGT_organisms}
\alias{IMGT_is_up}

\title{Install a germline db from IMGT}

\description{
  The \code{install_IMGT_germline_db()} function downloads V/D/J
  germline sequences from the IMGT website for a given organism, and
  stores them in a local germline database. This local database gets
  installed in \pkg{igblastr}'s persistent cache. It can then be used
  later with \code{\link{igblastn}()}.
}

\usage{
install_IMGT_germline_db(release, organism="Homo sapiens", tcr.db=FALSE,
                         force=FALSE, ...)

## Related utilities:
list_IMGT_releases(recache=FALSE)
list_IMGT_organisms(release)
IMGT_is_up()
}

\arguments{
  \item{release}{
    A single string specifying the IMGT/V-QUEST release to get the
    germline sequences from (or to list the organisms from for
    \code{list_IMGT_organisms()}).
    Use \code{list_IMGT_releases()} to list all releases.
  }
  \item{organism}{
    A single string specifying the latin name of the organism for which
    to get the germline sequences.
  }
  \item{tcr.db}{
    Should the database be populated with allele sequences from the
    BCR (B-cell Receptor) or TCR (T-cell Receptor) germline loci?

    The BCR germline loci are: IGH, IGK, IGL.

    The TCR germline loci are: TRA, TRB, TRG, TRD.

    By default, the V/D/J allele sequences from the BCR germline loci
    are downloaded. Set \code{tcr.db} to \code{TRUE} to download the
    V/D/J allele sequences from the TCR germline loci instead.
  }
  \item{force}{
    Set to \code{TRUE} to reinstall if the requested database is already
    installed.
  }
  \item{...}{
    Extra arguments to be passed to the internal call to \code{download.file()}.
    See \code{?\link[utils]{download.file}} in the \pkg{utils} package for
    more information.
  }
  \item{recache}{
    \code{list_IMGT_releases()} uses a caching mechanism so that the
    list of IMGT/V-QUEST releases gets downloaded only once from the IMGT
    website during an R session (note that this caching is done in memory
    so it does not persist across sessions).
    Set \code{recache} to \code{TRUE} to force a new download (and
    recaching) of the list of IMGT/V-QUEST releases.
  }
}

\details{
  The following naming scheme is used to form the name of the installed
  database:
  \preformatted{    IMGT-<release>.<organism>.<loci>}
  where:
  \enumerate{
    \item <release> is the IMGT/V-QUEST release e.g. 202518-3 or 202449-1.
          Use \code{list_IMGT_releases()} to get the list of releases
          currently available at IMGT/V-QUEST.
    \item <organism> is the latin name (a.k.a. binomial name) of the
          organism with all spaces replaced with underscores (_).
          For example Homo_sapiens or Macaca_mulatta.
          Use \code{list_IMGT_organisms("<release>")} to get the list
          of organisms included in a given IMGT/V-QUEST release.
          Note that, starting with release 202405-2, IMGT/V-QUEST
          provides BCR and TCR germline sequences for mouse strain
          C57BL6J (Mus_musculus_C57BL6J).
    \item <loci> is a string obtained by concatenating the germline
          loci together separated with the + sign. For example IGH+IGK+IGL
          or TRA+TRB+TRG+TRD.
          The list of loci depends on whether the germline sequences for
          BCR or TCR were requested. See \code{tcr.db} argument above for
          more information.
          Note that for some IMGT/V-QUEST releases/organisms, only a subset
          of the loci are available. For example, in release 202343-3,
          the only TCR germline loci available for Mus_musculus_C57BL6J
          are TRA and TRB. This will be automatically reflected in the
          name of the installed germline db.
  }
}

\value{
  \code{install_IMGT_germline_db()} returns the name to the newly
  installed germline db as an invisible string.

  \code{list_IMGT_releases()} returns the list of IMGT/V-QUEST releases
  in a character vector. The releases are sorted from newest to oldest
  (latest release is first).

  \code{list_IMGT_organisms()} returns the list of organisms included
  in the specified IMGT/V-QUEST release in a character vector.

  \code{IMGT_is_up()} returns \code{TRUE} or \code{FALSE}, indicating
  whether the IMGT website at \url{https://www.imgt.org} is up and running
  or down.
}

\note{
  \code{install_IMGT_germline_db()} generates the local
  database by performing the instructions provided at
  \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}.
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{use_germline_db}} to select the cached germline db
          to use with \code{igblastn()}.

    \item The IMGT website: \url{https://www.imgt.org/}.

    \item The IMGT/V-QUEST download site:
          \url{https://www.imgt.org/download/V-QUEST/}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

if (IMGT_is_up()) {
  ## As of March 26, 2025, the latest IMGT/V-QUEST release is 202518-3:
  list_IMGT_releases()

  list_IMGT_organisms("202518-3")

  ## Download Mouse BCR germline sequences from IMGT/V-QUEST 202518-3, and
  ## store them in a cached germline database:
  install_IMGT_germline_db("202518-3", organism="Mus musculus", force=TRUE)

  ## List the cached germline databases:
  list_germline_dbs()

  ## Select newly installed germline db to use with igblastn():
  use_germline_db("IMGT-202518-3.Mus_musculus.IGH+IGK+IGL")

  ## Download Mouse TCR germline sequences from IMGT/V-QUEST 202518-3, and
  ## store them in a cached germline database:
  install_IMGT_germline_db("202518-3", organism="Mus musculus",
                           tcr.db=TRUE, force=TRUE)
  list_germline_dbs()
}
}

\keyword{utilities}
