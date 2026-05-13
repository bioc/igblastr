\name{combine_germline_dbs}

\alias{combine_germline_dbs}

\title{Combine two existing germline dbs}

\description{
  Use the \code{combine_germline_dbs()} function to combine two existing
  germline databases, typically (but not necessarily) from two different
  organisms, into a new one.
  The new germline database gets added to \pkg{igblastr}'s persistent
  cache and can then be used later with \code{\link{igblastn}()}.
}

\usage{
combine_germline_dbs(db_name, db_name1, db_name2, suffix1, suffix2,
                     overwrite=FALSE, verbose=FALSE)
}

\arguments{
  \item{db_name}{
    A single string specifying the name of the germline db to install.
    The string \emph{must} start with "comb" and cannot contain whitespace
    characters.
  }
  \item{db_name1, db_name2}{
    COMING SOON...
  }
  \item{suffix1, suffix2}{
    COMING SOON...
  }
  \item{overwrite}{
    Set to \code{TRUE} if a germline db with the specified name is
    already installed in \pkg{igblastr}'s persistent cache, in which
    case the existing db will be replaced with the new one.
  }
  \item{verbose}{
    Set to \code{TRUE} to have the function display some details about
    its internal operations.
  }
}

\value{
  \code{combine_germline_dbs()} returns the name to the newly
  installed germline db as an invisible string.
}

\seealso{
  \itemize{
    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \code{\link{use_germline_db}} to select the cached germline db
          to use with \code{igblastn()}.

    \item \code{\link{load_germline_db}} to load the nucleotide sequences
          of the gene regions stored in a cached germline db.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## TO BE EXPANDED SOON...

db_name12 <- "comb.IMGT-202616-1.human+mouse.IGH+IGK+IGL"
db_name1 <- install_IMGT_germline_db("202616-1", "Homo_sapiens",
                                     overwrite=TRUE)
db_name2 <- install_IMGT_germline_db("202616-1", "Mus_musculus",
                                     overwrite=TRUE)
combine_germline_dbs(db_name12, db_name1, db_name2, "_Hs", "_Mm")

rm_germline_db(db_name12)
}

\keyword{utilities}
