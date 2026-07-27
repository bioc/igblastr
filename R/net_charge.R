### =========================================================================
### Net charges of amino acid sequences
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### net_charge()
###

.BASIC_GROUPS  <- c("R", "H", "K")
.ACIDIC_GROUPS <- c("D", "E", "C", "Y")

### pKa values taken from
###   https://www.chemverify.com/learn/amino-acid-reference-table
### Note that Wikipedia has slightly different values:
###   https://en.wikipedia.org/wiki/Protein_pKa_calculations
### A summary of measured pK values from 2008 also has slightly different
### values (and also value for Arg is missing):
###   https://pmc.ncbi.nlm.nih.gov/articles/PMC2708032/
PKA_VALUES <- c(
    ## Basic groups:
    nTer= 8.0,   # N-terminus (alpha-NH3+)
    R   =12.5,   # Arg
    H   = 6.0,   # His
    K   =10.5,   # Lys

    ## Acidic groups:
    D   = 3.65,  # Asp
    E   = 4.25,  # Glu
    C   = 8.3,   # Cys
    Y   =10.1,   # Tyr
    cTer= 3.1    # C-terminus (alpha-COO-)
)

.normarg_pKa_values <- function(pKa_values)
{
    if (!is.numeric(pKa_values))
        stop(wmsg("'pKa_values' must be a numeric vector"))
    nms <- names(pKa_values)
    if (is.null(nms))
        stop(wmsg("'pKa_values' must have names"))
    all_groups <- c(.BASIC_GROUPS, .ACIDIC_GROUPS)
    if (!all(all_groups %in% nms)) {
        in1string <- paste(all_groups, collapse=", ")
        stop(wmsg("'pKa_values' must have at least the following ",
                  "names on it: ", in1string))
    }
    if (anyNA(pKa_values))
        stop(wmsg("'pKa_values' cannot contain NAs"))
    pKa_values
}

.compute_basic_charges <- function(pH, pKa_values)
{
    stopifnot(isSingleNumber(pH), is.numeric(pKa_values), !anyNA(pKa_values))
    +1 / (1 + 10^(pH - pKa_values))
}

.compute_acidic_charges <- function(pH, pKa_values)
{
    stopifnot(isSingleNumber(pH), is.numeric(pKa_values), !anyNA(pKa_values))
    -1 / (1 + 10^(pKa_values - pH))
}

### We choose to set 'pH' to 7.4 (physiological pH) by default, rather
### than 7.0 (neutral pH).
net_charge <- function(aa_sequences, pH=7.4,
                       with.nTer=TRUE, with.cTer=TRUE,
                       pKa_values=PKA_VALUES, as.matrix=FALSE)
{
    if (is.character(aa_sequences)) {
        aa_sequences <- AAStringSet(aa_sequences)
    } else if (!is(aa_sequences, "AAStringSet")) {
        stop(wmsg("'aa_sequences' must be an AAStringSet object"))
    }
    if (!isSingleNumber(pH))
        stop(wmsg("'pH' must be a single number"))
    if (!isTRUEorFALSE(with.nTer))
        stop(wmsg("'with.nTer' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(with.cTer))
        stop(wmsg("'with.cTer' must be TRUE or FALSE"))
    pKa_values <- .normarg_pKa_values(pKa_values)
    if (!isTRUEorFALSE(as.matrix))
        stop(wmsg("'as.matrix' must be TRUE or FALSE"))
    basic_charges  <- .compute_basic_charges(pH, pKa_values[.BASIC_GROUPS])
    acidic_charges <- .compute_acidic_charges(pH, pKa_values[.ACIDIC_GROUPS])
    taf <- t(alphabetFrequency(aa_sequences))
    pos_charges <- t(taf[.BASIC_GROUPS , , drop=FALSE] * basic_charges)
    if (with.nTer) {
        nTer_pKa <- pKa_values[["nTer"]]
        nTer_charge <- .compute_basic_charges(pH,  nTer_pKa)
        pos_charges <- cbind(nTer=nTer_charge, pos_charges)
    }
    neg_charges <- t(taf[.ACIDIC_GROUPS, , drop=FALSE] * acidic_charges)
    if (with.cTer) {
        cTer_pKa <- pKa_values[["cTer"]]
        cTer_charge <- .compute_acidic_charges(pH,  cTer_pKa)
        neg_charges <- cbind(neg_charges, cTer=cTer_charge)
    }
    all_charges <- round(cbind(pos_charges, neg_charges), digits=3L)
    if (as.matrix)
        return(all_charges)
    rowSums(all_charges)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_region_net_charge()
###

extract_region_net_charge <- function(AIRR_df, region_type,
                                      pH=7.4, pKa_values=PKA_VALUES,
                                      as.matrix=FALSE)
{
    aa_sequences <- extract_sequence_region(AIRR_df, region_type, as.aa=TRUE)
    net_charge(aa_sequences, pH=pH, with.nTer=FALSE, with.cTer=FALSE,
               pKa_values=pKa_values, as.matrix=as.matrix)
}

