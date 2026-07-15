### =========================================================================
### Miscellaneous internal data utilities
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### write_intdata_to_db()
###

### Not exported!
### Note that we don't handle custom "pdm" data at the moment. So the
### function does not need a 'for.aa' argument. Also the name of the
### function reflects the fact that it only handles "ndm" data.
write_intdata_to_db <- function(ndm_data, db_path,
                                domain_system=c("imgt", "kabat"),
                                check.and.reorder=FALSE)
{
    stopifnot(is.data.frame(ndm_data), isTRUEorFALSE(check.and.reorder))
    domain_system <- match.arg(domain_system)
    intdata_path <- make_germline_db_intdata_path(db_path, FALSE, domain_system)
    intdata_dir <- dirname(intdata_path)
    stopifnot(!dir.exists(intdata_dir))

    ## Even though write_ndm_data() will call check_ndm_data_col2class()
    ## internally, we prefer to fail **before** creating the 'intdata_dir'
    ## folder.
    check_ndm_data_col2class(ndm_data)
    if (check.and.reorder) {
        db_V_fasta_file <- get_db_fasta_file(db_path, "V")
        db_V_allele_names <- names(fasta.seqlengths(db_V_fasta_file))
        ndm_data <- check_and_reorder_igdata_rows(ndm_data, db_V_allele_names)
    }
    stopifnot(dir.create(intdata_dir))
    write_ndm_data(ndm_data, intdata_path)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### V_genes_with_varying_fwrcdr_boundaries()
###

.extract_gene_names_as_factor <- function(intdata)
{
    allele_names <- get_intdata_col(intdata, "allele_name")
    gene_names <- allele2gene(allele_names)
    unique_gene_names <- unique(gene_names)
    factor(gene_names, levels=unique_gene_names)
}

.check_V_segment <- function(V_segment)
{
    if (!isSingleNonWhiteString(V_segment))
        stop(wmsg("'V_segment' must be a single (non-empty) string"))
    if (!(V_segment %in% V_GENE_SEGMENTS)) {
        in1string <- paste0("\"", V_GENE_SEGMENTS, "\"", collapse=", ")
        stop(wmsg("'V_segment' must be one of ", in1string))
    }
}

.V_genes_with_varying_segment_boundaries <- function(intdata, V_segment)
{
    f <- .extract_gene_names_as_factor(intdata)
    .check_V_segment(V_segment)
    starts <- get_intdata_col(intdata, paste0(V_segment, "_start"))
    ends <- get_intdata_col(intdata, paste0(V_segment, "_end"))
    starts_per_gene <- unique(splitAsList(starts, f))
    ends_per_gene <- unique(splitAsList(ends, f))
    levels(f)[lengths(starts_per_gene) != 1L | lengths(ends_per_gene) != 1L]
}

V_genes_with_varying_fwrcdr_boundaries <- function(intdata, V_segment=NULL)
{
    if (!is.null(V_segment))
        return(.V_genes_with_varying_segment_boundaries(intdata, V_segment))
    found_genes <- lapply(V_GENE_SEGMENTS,
        function(V_segment)
            .V_genes_with_varying_segment_boundaries(intdata, V_segment))
    found_genes <- unique(unlist(found_genes, use.names=FALSE))
    ## Return the gene names in the same order as they show up in 'intdata'.
    unique_gene_names <- levels(.extract_gene_names_as_factor(intdata))
    unique_gene_names[unique_gene_names %in% found_genes]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show_intdata_disagreements()
###

show_intdata_disagreements <- function(db_name)
{
    check_germline_db_name(db_name)
    db_intdata <- load_intdata(db_name)
    igblast_organism <- infer_igblast_organism_from_db_name(db_name)
    if (is.na(igblast_organism))
        stop(wmsg("The specified germline db does not seem to be ",
                  "for an IgBLAST organism. Note that you can use ",
                  "list_igblast_organisms() to get the list of ",
                  "IgBLAST organisms. See '?list_igblast_organisms' ",
                  "for more information."))
    igblast_intdata <- load_intdata(igblast_organism)
    diff <- df_diff(db_intdata, igblast_intdata, "allele_name",
                    "igblastr-generated", "IgBLAST-provided")
    what <- c("the igblastr-generated \"internal data\" included in this ",
              "germline db and the \"internal data\" provided by IgBLAST ",
              "for ", igblast_organism)
    if (length(diff) == 0L) {
        msg <- c("No disagreements between ", what, ".")
        cat(wmsg2(msg, margin=0L), "\n", sep="")
    } else {
        msg <- c("Disagreements between ", what, ":")
        cat(wmsg2(msg, margin=0L), "\n", sep="")
        cat(diff, sep="")
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### count_cysteines()
###

.get_amino_acids_at <- function(V_alleles, codon_ends)
{
    stopifnot(is(V_alleles, "DNAStringSet"), is.integer(codon_ends),
              length(V_alleles) == length(codon_ends))
    codon_starts <- codon_ends - 2L
    idx <- which(codon_starts >= 1L & codon_ends <= width(V_alleles))
    codons <- subseq(V_alleles[idx], codon_starts[idx], codon_ends[idx])
    translate(codons)
}

### Not exported!
### On the ungapped germline V gene allele sequences, the codons for cysteines
### 23 and 104 are expected to be found at positions fwr1_end-11 to fwr1_end-9
### for the former and fwr3_end-2 to fwr3_end for the latter. This is a quick
### and easy way to validate the intdata of a germline db.
### Returns a data.frame with 1 row per codon and columns "nb_cys", "other",
### and "percent_cys".
count_cysteines <- function(db_name)
{
    V_alleles <- load_germline_sequences(db_name, region_types="V")
    intdata <- load_intdata(db_name)
    stopifnot(identical(names(V_alleles), intdata[ , "allele_name"]))
    codon22  <- .get_amino_acids_at(V_alleles, intdata[ , "fwr1_end"] - 12L)
    codon23  <- .get_amino_acids_at(V_alleles, intdata[ , "fwr1_end"] - 9L)
    codon24  <- .get_amino_acids_at(V_alleles, intdata[ , "fwr1_end"] - 6L)
    codon25  <- .get_amino_acids_at(V_alleles, intdata[ , "fwr1_end"] - 3L)
    codon26  <- .get_amino_acids_at(V_alleles, intdata[ , "fwr1_end"])
    codon103 <- .get_amino_acids_at(V_alleles, intdata[ , "fwr3_end"] - 3L)
    codon104 <- .get_amino_acids_at(V_alleles, intdata[ , "fwr3_end"])
    codon105 <- .get_amino_acids_at(V_alleles, intdata[ , "fwr3_end"] + 3L)
    nb_cys <- c(sum(codon22 == "C"),
                sum(codon23 == "C"),
                sum(codon24 == "C"),
                sum(codon25 == "C"),
                sum(codon26 == "C"),
                sum(codon103 == "C"),
                sum(codon104 == "C"),
                sum(codon105 == "C"))
    other  <- c(sum(codon22 != "C"),
                sum(codon23 != "C"),
                sum(codon24 != "C"),
                sum(codon25 != "C"),
                sum(codon26 != "C"),
                sum(codon103 != "C"),
                sum(codon104 != "C"),
                sum(codon105 != "C"))
    lens   <- c(length(codon22),
                length(codon23),
                length(codon24),
                length(codon25),
                length(codon26),
                length(codon103),
                length(codon104),
                length(codon105))
    percent_cys <- round(100 * nb_cys / lens, digits=2L)
    ans <- data.frame(nb_cys, other, percent_cys)
    rownames(ans) <- paste0("codon", c(22:26, 103:105))
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### count_cysteines_for_IMGT_organisms()
###

### Not exported!
### Percent cysteines obtained with igblastr 1.3.13 for IMGT 202614-2:
### ---------------------------- tcr.db=FALSE ---------------------------
###                                              db_name codon23 codon104
### 1               IMGT-202614-2.Bos_taurus.IGH+IGK+IGL   98.88    95.00
### 2   IMGT-202614-2.Canis_lupus_familiaris.IGH+IGK+IGL   97.89    47.68
### 3               IMGT-202614-2.Equus_caballus.IGH+IGK   91.01    34.83
### 4                IMGT-202614-2.Gallus_gallus.IGH+IGL   73.98    71.09
### 5  IMGT-202614-2.Gorilla_gorilla_gorilla.IGH+IGK+IGL   98.26    93.40
### 6             IMGT-202614-2.Homo_sapiens.IGH+IGK+IGL   98.47    97.41
### 7              IMGT-202614-2.Lemur_catta.IGH+IGK+IGL   96.92   100.00
### 8              IMGT-202614-2.Macaca_fascicularis.IGH   93.18     0.00
### 9           IMGT-202614-2.Macaca_mulatta.IGH+IGK+IGL   95.62    34.65
### 10            IMGT-202614-2.Mus_musculus.IGH+IGK+IGL   98.49    94.51
### 11   IMGT-202614-2.Mustela_putorius_furo.IGH+IGK+IGL   65.54    32.20
### 12                   IMGT-202614-2.Neogale_vison.IGH    0.00     0.00
### 13             IMGT-202614-2.Oncorhynchus_mykiss.IGH   98.73    94.90
### 14        IMGT-202614-2.Ornithorhynchus_anatinus.IGH   97.78     0.00
### 15   IMGT-202614-2.Oryctolagus_cuniculus.IGH+IGK+IGL   98.64   100.00
### 16          IMGT-202614-2.Pongo_pygmaeus.IGH+IGK+IGL   49.66    48.11
### 17       IMGT-202614-2.Rattus_norvegicus.IGH+IGK+IGL   95.76    42.21
### 18                     IMGT-202614-2.Salmo_salar.IGH    0.00     0.67
### 19              IMGT-202614-2.Sus_scrofa.IGH+IGK+IGL   96.92    95.38
### 20                   IMGT-202614-2.Vicugna_pacos.IGH   98.81   100.00
### ------------------------------ tcr.db=TRUE ------------------------------
###                                                  db_name codon23 codon104
### 1               IMGT-202614-2.Bos_taurus.TRA+TRB+TRG+TRD   66.58    48.56
### 2      IMGT-202614-2.Camelus_dromedarius.TRA+TRB+TRG+TRD  100.00   100.00
### 3   IMGT-202614-2.Canis_lupus_familiaris.TRA+TRB+TRG+TRD   98.77    14.81
### 4                      IMGT-202614-2.Danio_rerio.TRA+TRD    0.71     0.00
### 5              IMGT-202614-2.Felis_catus.TRA+TRB+TRG+TRD   98.85    44.83
### 6  IMGT-202614-2.Gorilla_gorilla_gorilla.TRA+TRB+TRG+TRD   53.72    45.21
### 7    IMGT-202614-2.Heterocephalus_glaber.TRA+TRB+TRG+TRD   95.65    90.43
### 8             IMGT-202614-2.Homo_sapiens.TRA+TRB+TRG+TRD   96.88    97.69
### 9                  IMGT-202614-2.Macaca_fascicularis.TRB   95.45    98.46
### 10          IMGT-202614-2.Macaca_mulatta.TRA+TRB+TRG+TRD   34.26    62.79
### 11            IMGT-202614-2.Mus_musculus.TRA+TRB+TRG+TRD   20.80    20.00
### 12   IMGT-202614-2.Mustela_putorius_furo.TRA+TRB+TRG+TRD   96.64    36.97
### 13   IMGT-202614-2.Oryctolagus_cuniculus.TRA+TRB+TRG+TRD   99.32   100.00
### 14                  IMGT-202614-2.Ovis_aries.TRA+TRB+TRD   99.47    18.72
### 15             IMGT-202614-2.Pan_troglodytes.TRA+TRG+TRD   97.50    81.25
### 16            IMGT-202614-2.Pongo_abelii.TRA+TRB+TRG+TRD   95.83    85.00
### 17                  IMGT-202614-2.Pongo_pygmaeus.TRB+TRG   95.73    84.62
### 18                      IMGT-202614-2.Sus_scrofa.TRB+TRG   92.68     0.00
count_cysteines_for_IMGT_organisms <- function(release, tcr.db=FALSE)
{
    organisms <- list_IMGT_organisms(release)
    codon_names <- c("codon23", "codon104")
    percent_cys <- lapply(organisms,
        function(organism) {
            message(organism)
            db_name <- try(suppressWarnings(suppressMessages(
                install_IMGT_germline_db(release, organism,
                                         tcr.db=tcr.db, overwrite=TRUE)
            )), silent=TRUE)
            if (inherits(db_name, "try-error"))
                return(NULL)
            percents <- count_cysteines(db_name)[codon_names, "percent_cys"]
	    ## Returns a list with 3 components: db_name, codon23, codon104.
            c(list(db_name=db_name),
	      setNames(as.list(percents), codon_names))
        })
    percent_cys <- S4Vectors:::delete_NULLs(percent_cys)
    db_name  <- vapply(percent_cys, function(x) x$db_name, character(1))
    codon23  <- vapply(percent_cys, function(x) x$codon23, numeric(1))
    codon104 <- vapply(percent_cys, function(x) x$codon104, numeric(1))
    data.frame(db_name=db_name, codon23=codon23, codon104=codon104)
}

