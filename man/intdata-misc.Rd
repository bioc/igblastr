\name{intdata-misc}

\alias{intdata-misc}
\alias{intdata_misc}

\alias{V_genes_with_varying_fwrcdr_boundaries}
\alias{show_intdata_disagreements}

\title{Miscellaneous internal data utilities}

\description{
  \code{V_genes_with_varying_fwrcdr_boundaries()} can be used to identify
  germline V genes for which the FWR/CDR boundaries are not the same across
  all alleles.

  \code{show_intdata_disagreements()} shows the disagreements (if any)
  between the igblastr-generated \emph{internal data} included in a
  germline db and the \emph{internal data} provided by IgBLAST for the
  corresponding organism.
}

\usage{
V_genes_with_varying_fwrcdr_boundaries(intdata, V_segment=NULL)

show_intdata_disagreements(db_name)
}

\arguments{
  \item{intdata}{
    A data.frame as returned by \code{load_intdata()}.
  }
  \item{V_segment}{
    The name of a V gene segment. This can be set to \code{"fwr1"},
    \code{"cdr1"}, \code{"fwr2"}, \code{"cdr2"}, or \code{"fwr3"}.

    By default (i.e. when \code{V_segment} is omitted or set to \code{NULL}),
    \code{V_genes_with_varying_fwrcdr_boundaries()} will identify V genes for
    which any segment has varying boundaries across alleles. Otherwise, it
    will identify V genes for which the specified segment has varying
    boundaries.
  }
  \item{db_name}{
    A single string specifying the name of a cached germline db.
  }
}

\value{
  \code{V_genes_with_varying_fwrcdr_boundaries()} returns a character
  vector containing the names of the germline V genes for which the
  FWR/CDR boundaries are not the same across all alleles.

  \code{show_intdata_disagreements()} doesn't return anything (i.e.
  invisible \code{NULL}).
}

\seealso{
  \itemize{
    \item \code{\link{load_intdata}} to access IgBLAST \emph{internal data}.

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
## V_genes_with_varying_fwrcdr_boundaries()
## ---------------------------------------------------------------------

## Note that the alleles of a given germline V gene don't necessarily
## share the same FWR/CDR boundaries. You can use utility function
## V_genes_with_varying_fwrcdr_boundaries() to identify them:
human_intdata0 <- load_intdata("human")
var_genes <- V_genes_with_varying_fwrcdr_boundaries(human_intdata0)
var_genes  # 4 genes

## Display rows associated with gene IGHV4-31:
subset(human_intdata0, allele2gene(allele_name) == "IGHV4-31")

## Human germline V genes for which the CDR1 boundaries are not the
## same across all alleles:
var_genes <- V_genes_with_varying_fwrcdr_boundaries(human_intdata0,
                                                    V_segment="cdr1")
var_genes  # 2 genes

## Display rows associated with the genes in 'var_genes':
subset(human_intdata0, allele2gene(allele_name) \%in\% var_genes)

## ---------------------------------------------------------------------
## show_intdata_disagreements()
## ---------------------------------------------------------------------

db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"
show_intdata_disagreements(db_name)

## See last example in '?install_IMGT_germline_db' for an example of
## a germline db with disagreements.
}

\keyword{utilities}
