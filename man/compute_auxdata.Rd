\name{compute_auxdata}

\alias{compute_auxdata}
\alias{solve_cdr3_ends_using_fwr4_aa_comparisons}
\alias{solve_cdr3_ends_using_fwr4_dna_comparisons}
\alias{solve_cdr3_ends_using_fwr4_dna_PWM}

\title{Compute IgBLAST auxiliary data}

\description{
  Utility functions to generate \emph{auxiliary data} for a set of
  germline J gene allele sequences.

  See \code{?\link{load_auxdata}} for more information about IgBLAST
  \emph{auxiliary data}.
}

\usage{
## The primary utility for computing "auxiliary data":

compute_auxdata(J_alleles, codon_starts=NULL, no.warnings=FALSE)

## Additional utilities that can be used to complete the "auxiliary data"
## returned by compute_auxdata():

solve_cdr3_ends_using_fwr4_aa_comparisons(auxdata, J_alleles,
                                          max.dist=2, standout.by=5)

solve_cdr3_ends_using_fwr4_dna_comparisons(auxdata, J_alleles,
                                           max.dist=4, standout.by=15)

solve_cdr3_ends_using_fwr4_dna_PWM(auxdata, J_alleles,
                                   min.score=0.80, standout.by=0.40)
}

\arguments{
  \item{J_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline J
    gene allele sequences.
  }
  \item{codon_starts}{
    When supplied, must be a named integer vector \emph{parallel} to
    \code{J_alleles}, that is, it must have the same length and names
    as \code{J_alleles}. For each allele in \code{J_alleles},
    \code{codon_starts} should indicate the start position (1-based) of
    the first codon in the allele. This should be a number >= 1 and <= 3.
    \code{NA}s are allowed when the coding frame is not known.

    Note that the "coding frame start" is 0-based so can simply be inferred
    from the "codon start" by subtracting 1.
  }
  \item{no.warnings}{
    Set to \code{TRUE} to suppress the warnings that \code{compute_auxdata()}
    otherwise issues when the returned data.frame has \code{NA}s in its
    \code{coding_frame_start} column or negative values in its \code{cdr3_end}
    column.
  }
  \item{auxdata}{
    COMING SOON...
  }
  \item{max.dist}{
    COMING SOON...
  }
  \item{standout.by}{
    COMING SOON...
  }
  \item{min.score}{
    COMING SOON...
  }
}

\details{
  \code{compute_auxdata()} is meant to be used first, to generate
  \emph{auxiliary data} from scratch for a set of germline J gene
  allele sequences. Note that the data.frame returned by
  \code{compute_auxdata()} can contain \code{NA}s in its
  \code{cdr3_end} column. We call these \emph{unsolved CDR3 ends}.

  The \code{solve_cdr3_ends_using_fwr4_*()} functions are meant
  to be used on the \emph{auxiliary data} returned by
  \code{compute_auxdata()} to solve the unsolved CDR3 ends, if any.

  \subsection{The \code{compute_auxdata()} function}{
    \code{compute_auxdata()} uses a motif-based approach to identify
    the CDR3/FWR4 junction in the supplied J allele sequences.

    For human BCR germline J alleles, the FWR4 region is expected to start
    with the following amino acid motifs (X represents any amino acid):
    \itemize{
      \item "WGXG" on the heavy chain (IGH locus);
      \item "FGXG" on the light chain (IGK and IGL loci).
    }
    Things might deviate slightly for other organisms, or for TCR germline
    J alleles. For example, for TCR germline J alleles, the FWR4 region is
    expected to start with the "FGXG" motif, regardless of the locus.

    One complication is that the supplied J allele sequences are nucleotide
    sequences for which the coding frame might not be known.
    \code{compute_auxdata()} uses two distinct strategies to search for
    canonical motif "WGXG" or "FGXG" in a given J allele sequence, depending
    on whether the coding frame is known or not:
    \enumerate{
      \item If the coding frame is known (i.e. a vector of codon starts
            is passed thru the \code{codon_starts} argument and it reports
            a non-\code{NA} value for the allele), then this information is
            used to translate the coding region of the allele to an amino
            acid sequence. The function then tries to locate canonical
            motif "WGXG" or "FGXG" on the translated sequence.

            Note that if the allele is from the IGH locus, then canonical
            motif "WGXG" is used. For all other loci (IGK, IGL, TRA, TRB,
            TRG, TRD loci), canonical motif "FGXG" is used.

      \item If the coding frame is not known (i.e. the \code{codon_starts}
            argument is not set or it's set but the codon start reported
            for the allele is \code{NA}), then the search above is performed
            for all 3 possible coding frames, that is, for all 3 possible
            codon start values (1, 2, and 3).
    }
    Next, the CDR3/FWR4 junction is identified based on the result of the
    motif search (and the \code{cdr3_end} set accordingly):
    \itemize{
      \item If the motif search didn't return any match, then the CDR3/FWR4
            junction is considered to remain unknown. In this case, the
            \code{cdr3_end} is set to \code{NA}.

      \item If the motif search returned exactly one match, then the match
            is considered to mark the start of the FWR4. This means that
            the CDR3/FWR4 junction is now identified, and so the
            \code{cdr3_end} is set accordingly.

      \item If the motif search returned more than one match, then only the
            first match is considered to mark the start of the FWR4.
            All downstream matches are ignored.
    }

    Finally, the \code{coding_frame_start} is obtained as follow. For a given
    J allele:
    \enumerate{
      \item If a vector of codon starts was passed thru the \code{codon_starts}
            argument and it reports a non-\code{NA} value for the allele,
            then the \code{coding_frame_start} (0-based) is set to the
            reported codon start (1-based) minus 1.

      \item If the above cannot be performed, but if the CDR3/FWR4 junction
            of the allele was successfully identified by the motif search,
            then the \code{coding_frame_start} is inferred from the
            \code{cdr3_end}.

      \item If none of the two operations above can be performed, then the
            \code{coding_frame_start} is set to \code{NA} with a warning.
    }
  }

  \subsection{The \code{solve_cdr3_ends_using_fwr4_*()} functions}{
    COMING SOON...
  }
}

\value{
  All the functions documented in this man page return a data.frame
  with 1 row per supplied J allele and the same columns as the data.frame
  returned by \code{read_auxdata()}.
  See \code{?\link{read_auxdata}} for more information.
}

\seealso{
  \itemize{
    \item \code{\link{print_J_alleles}} to display a set of J alleles
          sequences where the sequences are justified w.r.t. their
          CDR3/FWR4 junction.

    \item \code{\link{load_auxdata}} to access IgBLAST \emph{auxiliary data}.

    \item \code{\link{compute_V_gene_delineations}} to annotate a set of
          germline V gene allele sequences.

    \item \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
          for important information about IgBLAST \emph{auxiliary data}.

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
## BASIC EXAMPLE
## ---------------------------------------------------------------------

## Let's load a set of human J allele sequences:
db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"
J_alleles <- load_germline_sequences(db_name, region_types="J")
J_alleles  # DNAStringSet object

computed_auxdata <- compute_auxdata(J_alleles)
head(computed_auxdata)

## The CDR3/FWR4 junction was found for all the alleles:
anyNA(computed_auxdata$cdr3_end)  # FALSE
print_J_alleles(J_alleles, computed_auxdata, translate=TRUE)

## ---------------------------------------------------------------------
## SANITY CHECK
## ---------------------------------------------------------------------

## Note that 'computed_auxdata' is in agreement with the auxiliary
## data included in IgBLAST for human, except for alleles IGHJ6*02
## and IGHJ6*03 (the 'extra_bps' column contains incorrect values
## for these two alleles, 0's instead of 1's):
human_auxdata0 <- load_auxdata("human", which="original")
bad_alleles <- c("IGHJ6*02", "IGHJ6*03")
subset(human_auxdata0, allele_name \%in\% bad_alleles)

## We manually correct this:
fixme <- human_auxdata0[ , "allele_name"] \%in\% bad_alleles
human_auxdata0[fixme, "extra_bps"] <- 1L  # replace 0L with 1L

## Now the data in 'computed_auxdata' matches exactly the corresponding
## data in 'human_auxdata0':
m <- match(computed_auxdata[ , "allele_name"],
           human_auxdata0[ , "allele_name"])
human_auxdata <- human_auxdata0[m, ]
rownames(human_auxdata) <- NULL
stopifnot(identical(computed_auxdata, human_auxdata))
}

\keyword{utilities}
