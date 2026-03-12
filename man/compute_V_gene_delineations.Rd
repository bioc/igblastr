\name{compute_V_gene_delineations}

\alias{compute_V_gene_delineations}
\alias{compute_imgt_intdata}

\title{Compute V gene delineations}

\description{
  A utility function to annotate a set of germline V gene allele sequences
  in a way similar to how they are annotated in IgBLAST \emph{internal data}.

  Note that the annotation produced by the function can be used by IgBLAST
  as a substitute to the \emph{internal data} shipped with the IgBLAST
  software (and typically included in a standard IgBLAST installation).

  See \code{?\link{load_intdata}} for more information about IgBLAST
  \emph{internal data}.
}

\usage{
compute_V_gene_delineations(gapped_V_alleles, as.IRangesList=FALSE)
}

\arguments{
  \item{gapped_V_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline V
    gene allele \emph{gapped} sequences.
  }
  \item{as.IRangesList}{
    COMING SOON...
  }
}

\details{
  COMING SOON...
}

\value{
  Returns the computed \emph{internal data} in a data.frame with 1 row
  per supplied V allele sequence and a set of columns similar to the
  data.frame returned by \code{\link{read_ndm_data}()}. With the following
  differences:
  \itemize{
    \item The data.frame returned by \code{compute_V_gene_delineations()} has
          no \code{chain_type} column.
    \item The data.frame returned by \code{compute_V_gene_delineations()} has
          the following additional columns:
          \enumerate{
            \item \code{fwr3_trimmed}: how many nucleotides a sequence is
                  missing on its 3' end in order to have a complete FWR3;
            \item \code{starting_gap}: size (in number of nucleotides) of
                  gap block located at the very beginning of the sequence
                  if any (0 if no such block);
            \item \code{all_gaps_in_frame}: indicates whether the gap blocks
                  align with the underlying coding frame or not;
            \item \code{all_gaps_contained}: \code{TRUE} if the gap blocks
                  don't cross the FWR/CDR boundaries, and \code{FALSE}
                  otherwise.
          }
  }
}

\seealso{
  \itemize{
    \item \code{\link{load_intdata}} to access and manipulate IgBLAST
          \emph{internal data}.

    \item \code{\link{compute_auxdata}} to annotate a set of germline J
          gene allele sequences.

    \item \link[Biostrings]{DNAStringSet} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC USAGE
## ---------------------------------------------------------------------
## Visit https://ogrdb.airr-community.org/germline_sets/Homo%20sapiens
## for the list of germline sets available for human at OGRDB.

## Download OGRDB germline gene allele sequences for human locus IGK
## to a temporary directory:
germline_set <- c(`IGKappa_VJ`=4)

dir.create(destdir <- tempfile())
download_OGRDB_germline_sequences("Homo sapiens", germline_set,
                                  destdir=destdir)

## List the downloaded files:
list.files(destdir)  # 2 FASTA files (one per IMGT group)

## The V allele sequences are **gapped**:
IGKV_alleles <- readDNAStringSet(file.path(destdir, "IGKV.fasta"))
as.character(IGKV_alleles[1:3])

## Compute the V gene delineations for 'IGKV_alleles':
IGKV_gene_delinea <- compute_V_gene_delineations(IGKV_alleles)
head(IGKV_gene_delinea)

## ---------------------------------------------------------------------
## COMPARE WITH V GENE DELINEATIONS PROVIDED BY AIRR-C JSON FILE
## ---------------------------------------------------------------------
## As a sanity check, we're going to compare 'IGKV_gene_delinea' with
## the V gene delineations provided by the corresponding AIRR-C JSON
## file from OGRDB.

download_OGRDB_germline_json("Homo sapiens", germline_set,
                             destdir=destdir)

list.files(destdir)  # 1 JSON file

ogrdb_json_file <- file.path(destdir, "IGKappa_VJ.json")
IGKV_gene_delinea2 <- makeogrannote(ogrdb_json_file)
head(IGKV_gene_delinea2)

## Keep only the shared columns:
shared_colnames <- intersect(colnames(IGKV_gene_delinea),
                             colnames(IGKV_gene_delinea2))
IGKV_gene_delinea  <- IGKV_gene_delinea[ , shared_colnames]
IGKV_gene_delinea2 <- IGKV_gene_delinea2[ , shared_colnames]

## We use same_ndm_data(), a small undocumented helper function, to
## verify that 'IGKV_gene_delinea' and 'IGKV_gene_delinea2' contain
## the same delineation data (possibly with their rows in different
## order):
stopifnot(same_ndm_data(IGKV_gene_delinea, IGKV_gene_delinea2))

## Remove temporary directory:
unlink(destdir, recursive=TRUE)
}

\keyword{utilities}
