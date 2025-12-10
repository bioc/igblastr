\name{igdata-utils}

\alias{igdata-utils}
\alias{igdata_utils}

\alias{get_igblast_auxiliary_data}
\alias{load_igblast_auxiliary_data}

\alias{update_live_igdata}
\alias{igdata_info}
\alias{time_since_live_igdata_last_checked}
\alias{reset_live_igdata}

\alias{print.aux_data_md5sum_df}
\alias{print.igdata_info}

\title{Low-level manipulation of IgBLAST auxiliary and internal data}

\description{
  A small set of low-level utility functions to access, manipulate, and
  update IgBLAST auxiliary and internal data.
}

\usage{
## Access and manipulate IgBLAST auxiliary data:
get_igblast_auxiliary_data(organism, which=c("live", "original"))
load_igblast_auxiliary_data(organism, which=c("live", "original"))

## Update IgBLAST auxiliary and internal data:
update_live_igdata(check.only=FALSE)
igdata_info()
time_since_live_igdata_last_checked(units="days")
reset_live_igdata(subdirs=c("all", "internal_data", "optional_file"))
}

\arguments{
  \item{organism}{
    A single string containing the name of an organism as
    returned by \code{\link{list_igblast_organisms}()}.
  }
  \item{which}{
    By default, \code{get_igblast_auxiliary_data()} and
    \code{load_igblast_auxiliary_data()} access the "live IgBLAST data",
    that is, the IgBLAST data that the user has possibly updated with
    \code{update_live_igdata()}. Depending on whether updates
    were applied or not, the "live IgBLAST data" might differ from
    the original IgBLAST data.

    Set \code{which} to \code{"original"} if you want to access
    the original IgBLAST data instead.
  }
  \item{check.only}{
    By default, \code{update_live_igdata()} checks for new IgBLAST
    auxiliary or internal data files available at NCBI, and it installs
    them if any are found. Set \code{check.only} to \code{TRUE} to only
    do the check without installing anything.
  }
  \item{units}{
    See \code{?base::\link{difftime}} for valid units.
  }
  \item{subdirs}{
    By default, \code{reset_live_igdata()} resets both \code{internal_data/}
    and \code{optional_file/} directories to their original states.
    Set \code{subdirs} to \code{"internal_data"} or \code{"optional_file"}
    to reset only a particular directory.
  }
}

\details{
  \subsection{Auxiliary data and internal data}{
    A standard IgBLAST installation -- like the one used by the
    \pkg{igblastr} package -- typically includes \emph{auxiliary data}
    and \emph{internal data} that are normally found in directories
    \code{internal_data/} and \code{optional_file/}, respectively.
    Both directories should be subdirectories of the \emph{root directory}
    of the IgBLAST installation, that is, of the directory returned by
    \code{\link{get_igblast_root}()}.

    We sometimes refer to this data simply as the \emph{IgBLAST data}.
  }

  \subsection{Access auxiliary data for a given organism}{
    The \emph{auxiliary data} consists of one tabulated file per organism.
    Each file indicates the germline J gene coding frame start position,
    the J gene type, and the CDR3 end position for each sequence in the
    germline J sequence database. See
    \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html} for
    the details.

    You can use \code{get_igblast_auxiliary_data()} to obtain the path
    to the file containing the auxiliary data for a given organism.
  }

  \subsection{NCBI updates}{
    NCBI occasionally updates some of the files in the
    \code{internal_data/} and \code{optional_file/} directories between
    IgBLAST releases, and it is recommended to use the new files.
    They make the new files available at
    \url{https://ftp.ncbi.nih.gov/blast/executables/igblast/release/patch/}.

    To download and install these new files, simply call
    \code{update_live_igdata()}. This will check for new IgBLAST
    auxiliary or internal data files available at NCBI, and install
    them if any are found.

    You can restore the original files at any moment with
    \code{reset_live_igdata()}.
  }
}

\value{
  \code{get_igblast_auxiliary_data()} returns a single string containing
  the path to the auxiliary data included in the IgBLAST installation used
  by \pkg{igblastr}, for the specified organism. Not necessarily suitable
  to use with \code{\link{igblastn}()} (see WARNING below).

  \code{load_igblast_auxiliary_data()} returns the auxiliary data in a
  data.frame with 1 row per germline J sequence and the following columns:
  \enumerate{
    \item \code{sseqid}: gene/allele name a.k.a. subject sequence id;
    \item \code{coding_frame_start}: first coding frame start
          position (position is 0-based);
    \item \code{chaintype}: chain type;
    \item \code{CDR3_stop}: CDR3 stop;
    \item \code{extra_bps}: extra base pairs beyond J coding end.
  }

  \code{update_live_igdata()} and \code{reset_live_igdata()} don't return
  anything (invisible \code{NULL}).

  \code{igdata_info()} returns a named list containing information about
  the state of the "live" and "original" IgBLAST data.

  \code{time_since_live_igdata_last_checked()} returns the time passed
  since the last run of \code{update_live_igdata()} in the specified
  units (days by default).
}

\section{WARNING}{
  According to \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
  the auxiliary data included in IgBLAST is specific to a particular
  NCBI or IMGT germline db. Unfortunately this means that this data is
  NOT guaranteed to be compatible with the germline db that you will
  use with \code{\link{igblastn}()}. See documentation of the
  \code{auxiliary_data} argument in \code{?\link{igblastn}} for
  more information about this.
}

\seealso{
  \itemize{
    \item \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
          for important information about the IgBLAST auxiliary data.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{install_igblast}} to perform an \emph{internal}
          IgBLAST installation.

    \item \code{\link{get_igblast_root}} to get (or set) the IgBLAST
          installation used (or to be used) by the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

igblast_info()

## ---------------------------------------------------------------------
## Access and manipulate IgBLAST auxiliary data
## ---------------------------------------------------------------------

list_igblast_organisms()

get_igblast_auxiliary_data("human")
load_igblast_auxiliary_data("human")

get_igblast_auxiliary_data("rhesus_monkey")
load_igblast_auxiliary_data("rhesus_monkey")

## ---------------------------------------------------------------------
## Check for NCBI updates
## ---------------------------------------------------------------------

igdata_info()
update_live_igdata(check.only=TRUE)
igdata_info()

## ---------------------------------------------------------------------
## "live" vs "original" IgBLAST data
## ---------------------------------------------------------------------

## By default, the "live IgBLAST data" gets accessed or returned:
get_igblast_auxiliary_data("human")
live_human_aux <- load_igblast_auxiliary_data("human")

## Access the original IgBLAST data:
get_igblast_auxiliary_data("human", which="original")
orig_human_aux <- load_igblast_auxiliary_data("human", which="original")

## "live" and "original" IgBLAST data can differ if the former was
## updated with update_live_igdata(). Otherwise, they'll be the same:
identical(live_human_aux, orig_human_aux)
}

\keyword{utilities}
