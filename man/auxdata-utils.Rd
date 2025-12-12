\name{auxdata-utils}

\alias{auxdata-utils}
\alias{auxdata_utils}

\alias{get_auxdata_path}
\alias{load_auxdata}
\alias{compute_auxdata}

\alias{get_igblast_auxiliary_data}
\alias{load_igblast_auxiliary_data}

\title{Access or generate IgBLAST auxiliary data}

\description{
  A small set of utilities to access or generate IgBLAST auxiliary data.
}

\usage{
get_auxdata_path(organism, which=c("live", "original"))

load_auxdata(organism, which=c("live", "original"))

compute_auxdata(J_alleles)
}

\arguments{
  \item{organism}{
    A single string containing the name of an organism as
    returned by \code{\link{list_igblast_organisms}()}.
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
  \item{J_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing the germline J
    gene allele sequences for which to compute the auxiliary data.
  }
}

\details{
  The \emph{auxiliary data} is typically included in a standard IgBLAST
  installation. It's located in the \code{optional_file/} directory which
  is itself a subdirectory of IgBLAST \emph{root directory}.
  It consists of one tabulated file per organism. Each file indicates the
  germline J gene coding frame start position, the J gene type, and the
  CDR3 end position for all known germline J allele sequences. See
  \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html} for
  additional details.
}

\value{
  \code{get_auxdata_path()} returns a single string containing
  the path to the auxiliary data included in the IgBLAST installation used
  by \pkg{igblastr}, for the specified organism. Not necessarily suitable
  to use with \code{\link{igblastn}()} (see WARNING below).

  \code{load_auxdata()} returns the auxiliary data in a data.frame
  with 1 row per germline J allele sequence and the following columns:
  \enumerate{
    \item \code{allele_name}: allele name;
    \item \code{coding_frame_start}: first coding frame start
          position (position is 0-based);
    \item \code{chain_type}: chain type;
    \item \code{CDR3_stop}: CDR3 stop (a.k.a. end) position (position
          is 0-based);
    \item \code{extra_bps}: extra base pairs beyond J coding end.
  }

  \code{compute_auxdata()} returns the auxiliary data in a data.frame
  with 1 row per supplied germline J allele sequence and the same columns
  as above, plus the \code{has_stop_codon} sequence.
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
    \item \code{\link{update_live_igdata}} for more information about "live"
          and "original" IgBLAST data.

    \item \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
          for important information about the IgBLAST auxiliary data.

    \item \link[Biostrings]{DNAStringSet} and \link[Biostrings]{AAStringSet}
          objects in the \pkg{Biostrings} package.

    \item The \code{\link[Biostrings]{translate}} function in the
          \pkg{Biostrings} package.

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
## Access and load IgBLAST auxiliary data for a given organism
## ---------------------------------------------------------------------

list_igblast_organisms()

get_auxdata_path("human")
load_auxdata("human")

get_auxdata_path("rabbit")
rabbit_aux <- load_auxdata("rabbit")

## It turns out that IgBLAST auxiliary data for rabbit matches exactly
## the set of germline J alleles for rabbit available at IMGT:
db_name <- install_IMGT_germline_db("202531-1", "Oryctolagus cuniculus",
                                    force=TRUE)
J_alleles <- load_germline_db(db_name, region_types="J")
stopifnot(setequal(names(J_alleles), rabbit_aux$allele_name))

## Note that this might change with future IMGT releases.

## Let's put the allele sequences in 'J_alleles' in the same order than
## in 'rabbit_aux':
J_alleles <- J_alleles[rabbit_aux$allele_name]
stopifnot(identical(names(J_alleles), rabbit_aux$allele_name))

## The 'coding_frame_start' column in 'rabbit_aux' contains integer
## values that are >= 0 and < 3. They indicate how many nucleotides
## precede the first codon on each allele sequence. In other words,
## this is the number of nucleotides that we need to trim on the 5'
## end of the germline J allele sequence before we start translating
## it. Let's use this information to translate the DNA sequences
## in 'J_alleles':
trimmed <- narrow(J_alleles, start=rabbit_aux$coding_frame_start + 1L)
aa <- suppressWarnings(translate(trimmed, no.init.codon=TRUE,
                                 if.fuzzy.codon="error.if.X"))

## 'aa' is an AAStringSet object:
aa

## No sequence in 'aa' should contain the letter "*" that is used
## by translate() to represent a stop codon. However, one J allele
## in IMGT-202531-1.Oryctolagus_cuniculus.IGH+IGK+IGL seems to disobey:
has_stop_codon <- vcountPattern("*", aa) != 0L
rabbit_aux[has_stop_codon, ]  # coding_frame_start = 0 for IGKJ1-2*04
J_alleles[has_stop_codon]     # first codon (TGA) is a stop codon
aa[has_stop_codon]            # indeed!

## ---------------------------------------------------------------------
## Compute auxiliary data for a set of germline J gene allele sequences
## ---------------------------------------------------------------------

## On the heavy chain, the FWR4 region is expected to start with amino
## acid pattern "WGXG" (X represents any amino acid). On the light chain,
## it's expected to start with pattern "FGXG". compute_auxdata() uses
## this to determine CDR3 end position ('CDR3_stop' column) and infer
## the 'coding_frame_start' and 'extra_bps' columns.
## It will also emit a warning if the CDR3 stop could not found be found
## for some alleles or a stop codon was found in some alleles.

computed_rabbit_aux <- compute_auxdata(J_alleles)
computed_rabbit_aux

## Alleles for which the CDR3 stop could not found:
J_alleles[is.na(computed_rabbit_aux$CDR3_stop)]

## Alleles with a stop codon:
J_alleles[which(computed_rabbit_aux$has_stop_codon)]  # IGKJ1-2*04 again

## 'computed_rabbit_aux' is in agreement with 'rabbit_aux', except for
## the 8 alleles for which the CDR3 stop could not be found:
keep_idx <- which(!is.na(computed_rabbit_aux$CDR3_stop))
stopifnot(
    identical(computed_rabbit_aux[keep_idx, -6], rabbit_aux[keep_idx, ])
)
}

\keyword{utilities}
