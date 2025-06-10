\name{igblastr_usage_report}

\alias{igblastr_usage_report}
\alias{BLAST_USAGE_REPORT}
\alias{Usage_reporting}
\alias{usage_reporting}
\alias{Usage_report}
\alias{usage_report}

\title{Turn "Usage Reporting" on or off}

\description{
  By default, the NCBI BLAST+ and IgBLAST programs will "call home" to
  report usage when they run on a computer with internet access. See
  \url{https://www.ncbi.nlm.nih.gov/books/NBK569851/} for the details.
  This can induce a significant slowdown in some situations e.g. when
  the \code{igblastn} \emph{standalone executable} is called in a loop
  on a small set of query sequences at each iteration.

  For this reason, the "call home" feature is disabled in \pkg{igblastr}
  by default, unless environment variable BLAST_USAGE_REPORT is set to true.

  More precisely, the "call home" feature is controlled by global option
  \code{igblastr_usage_report} in \pkg{igblastr}.
  On package startup, this option is set to \code{TRUE} if environment
  variable BLAST_USAGE_REPORT is set to true. Otherwise (i.e. if
  BLAST_USAGE_REPORT is not set, or is set to false or gibberish) it
  is set to \code{FALSE}.
}

\details{
  The user can change the value of global option \code{igblastr_usage_report}
  any time with:
  \preformatted{    options(igblastr_usage_report=TRUE)}
  or with:
  \preformatted{    options(igblastr_usage_report=FALSE)}
  To get the value of this option, use:
  \preformatted{    getOption("igblastr_usage_report")}

  Note that changing the value of a global option interactively
  with \code{options(...)} won't be remembered across R sessions.
  For a persistent change, you can either:
  \itemize{
    \item Put the \code{options(...)} command in your .Rprofile file. See
          \code{?\link{Rprofile}} for more information. Note that this is
          the standard way of setting global options persistently.

    \item In the particular case of global option \code{igblastr_usage_report}
          an alternative is to define environment variable BLAST_USAGE_REPORT
          outside R. The exact way to do this is OS-dependent e.g. on Linux
          and Mac you can define it in your user's \code{.profile} by adding
          the following line to it:
          \preformatted{    export BLAST_USAGE_REPORT=true}
  }
}

\seealso{
  \itemize{
    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
## Check current status of usage reporting:
getOption("igblastr_usage_report")

## Turn on usage reporting:
options(igblastr_usage_report=TRUE)

## Turn off usage reporting:
options(igblastr_usage_report=FALSE)
}

\keyword{utilities}
