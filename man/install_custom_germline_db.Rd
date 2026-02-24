\name{install_custom_germline_db}

\alias{install_custom_germline_db}

\title{Install a germline db from user-supplied gene allele sequences}

\description{
  The \code{install_custom_germline_db()} function creates a germline
  database from a set of user-supplied FASTA files containing germline
  V/D/J gene allele sequences, and installs it in \pkg{igblastr}'s
  persistent cache. The new germline database can then be used later
  with \code{\link{igblastn}()}.
}

\usage{
install_custom_germline_db(db_name, fasta_dir, tcr.db=FALSE,
                           loci="auto", gapped=FALSE, with.intdata=FALSE,
                           overwrite=FALSE, verbose=FALSE)
}

\arguments{
  \item{db_name}{
    A single string specifying the name of the germline db to install.
    The string \emph{must} start with cu and cannot contain whitespace
    characters.
  }
  \item{fasta_dir}{
    The path to the directory containing the FASTA files to import in
    the database. Note that the files to import \emph{must} be named:
    \itemize{
      \item IGH[VDJ].fasta, IGK[VJ].fasta, and IGL[VJ].fasta for
            a BCR germline loci database;
      \item TRA[VJ].fasta, TRB[VDJ].fasta, TRG[VJ].fasta, and
            TRD[VDJ].fasta for a TCR germline loci database.
    }
    See the \code{tcr.db} argument below for more details.

    Note that the FASTA header lines are expected to contain the allele
    names:
    \itemize{
      \item either on their own if the header lines have no \code{|},
      \item or in the 2nd field if the header lines have fields
            separated with \code{|}.
    }
  }
  \item{tcr.db}{
    Should the database to install be populated with allele sequences
    from the BCR (B-cell Receptor) or TCR (T-cell Receptor) germline loci?

    The BCR germline loci are: IGH, IGK, IGL.

    The TCR germline loci are: TRA, TRB, TRG, TRD.

    By default, the V/D/J allele sequences from the BCR germline loci
    will be used. More precisely, the 7 following FASTA files will be
    expected to be present in \code{fasta_dir}:
    \preformatted{        IGHV.fasta, IGHD.fasta, IGHJ.fasta,
        IGKV.fasta, IGKJ.fasta,
        IGLV.fasta, IGLJ.fasta}
    The files that are effectively found will be imported in the database,
    with a warning if some of them are missing.

    If \code{tcr.db} is set to \code{TRUE}, then the V/D/J allele sequences
    from the TCR germline loci will be used instead. More precisely, the 10
    following FASTA files will be expected to be present in \code{fasta_dir}:
    \preformatted{        TRAV.fasta, TRAJ.fasta,
        TRBV.fasta, TRBD.fasta, TRBJ.fasta,
        TRGV.fasta, TRGJ.fasta,
        TRDV.fasta, TRDD.fasta, TRDJ.fasta}
    The files that are effectively found will be imported in the database,
    with a warning if some of them are missing.
  }
  \item{loci}{
    By default, the database to install will be populated with the allele
    sequences from all the BCR or TCR loci.
    However, if you want to restrict the database to specific loci,
    you can use the \code{loci} argument to specify these loci.
    The subset of loci can be specified as a character vector with one
    element per locus (e.g. \code{"IGH"} or \code{c("TRA", "TRB")}), or
    as a \code{+}-separated list in a single string (e.g. \code{"TRA+TRB"}).
  }
  \item{gapped}{
    Set to \code{TRUE} if the supplied V allele sequences are \emph{gapped},
    in which case the gaps will be removed when the sequences get imported
    in the database.
  }
  \item{with.intdata}{
    If the supplied V allele sequences are \emph{gapped}, you can
    set \code{with.intdata} to \code{TRUE} to have the \emph{internal data}
    associated with the V alleles computed and added to the new db.

    Note that the presence of \emph{internal data} in the new db will be
    indicated in the \code{intdata} column of the data.frame returned by
    \code{\link{list_germline_dbs}()}.

    If a cached germline db includes its own \emph{internal data}, then
    \code{\link{igblastn}()} will use it instead of the \emph{internal data}
    shipped with IgBLAST. See documentation of the \code{custom_internal_data}
    argument in \code{?\link{igblastn}} for more information.

    Also see \code{?\link{intdata_utils}} for more information about IgBLAST
    \emph{internal data}.
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

\details{
}

\value{
  \code{install_custom_germline_db()} returns the name to the newly
  installed germline db as an invisible string.
}

\note{
  \code{install_custom_germline_db()} creates the local
  database by performing the instructions provided at
  \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}.
}

\seealso{
  \itemize{
    \item \code{\link{install_IMGT_germline_db}} to install a germline db
          from IMGT.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item \code{\link{list_germline_dbs}} to list the cached germline dbs.

    \item \code{\link{load_germline_db}} to load the nucleotide sequences
          of the gene regions stored in a cached germline db.

    \item \code{\link{use_germline_db}} to select the cached germline db
          to use with \code{igblastn()}.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
if (!has_igblast()) install_igblast()

## In this admittedly artificial example, we're going to create and
## install a germline db from the gene allele sequences in built-in
## germline db _AIRR.human.IGH+IGK+IGL.202410. More precisely, we're
## going to:
##   1. Dump built-in germline db _AIRR.human.IGH+IGK+IGL.202410
##      in 7 FASTA files.
##   2. Use these FASTA files to create and install a new germline db.
##   3. Check that the new db contains the same sequences as the
##      original db.
##   4. Remove the new db.

## --- 1. Dump _AIRR.human.IGH+IGK+IGL.202410 in 7 FASTA files ---

dump_germline_db <- function(db_name, destdir) {
  for (region_type in c("V", "D", "J")) {
    alleles <- load_germline_db(db_name, region_types=region_type)
    for (locus in c("IGH", "IGK", "IGL")) {
      loc_alleles <- alleles[grep(locus, names(alleles))]
      if (length(loc_alleles) == 0L)
        next
      fasta_file <- paste0(locus, region_type, ".fasta")
      writeXStringSet(loc_alleles, file.path(test_fasta_dir, fasta_file))
    }
  }
}

db_name <- "_AIRR.human.IGH+IGK+IGL.202410"
dir.create(test_fasta_dir <- tempfile())

dump_germline_db(db_name, test_fasta_dir)

list.files(test_fasta_dir)  # 7 FASTA files

## --- 2. Use FASTA files to install new germline db ---

custom_db_name <- "cu.test.human.IGH+IGK+IGL"
install_custom_germline_db(custom_db_name, test_fasta_dir)
list_germline_dbs()

## --- 3. Check that new and original dbs have the same sequences ---

stopifnot(identical(
    as.character(load_germline_db(db_name)),
    as.character(load_germline_db(custom_db_name))
))

## --- 4. Remove the new db ---

rm_germline_db(custom_db_name)
}

\keyword{utilities}
