### =========================================================================
### Miscellaneous auxiliary data utilities
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compute_germline_db_auxdata()
###

.do_compute_germline_db_auxdata <- function(db_path, ...)
{
    J_alleles <- readDNAStringSet(get_db_fasta_file(db_path, "J"))
    compute_auxdata(J_alleles, ...)
}

compute_germline_db_auxdata <- function(db_name, ...)
{
    check_germline_db_name(db_name)
    db_path <- get_germline_db_path(db_name)
    .do_compute_germline_db_auxdata(db_path, ...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### load_and_fix_igblast_auxdata()
###

### Not exported but used in various unit tests!
### The original IgBLAST auxiliary data has known problems that this function
### fixes on the fly. Note that some of these problems -- but not all! -- have
### been corrected in the updated _gl.aux files that the NCBI folks made
### available in the patch subdirectory here:
### https://ftp.ncbi.nih.gov/blast/executables/igblast/release/patch/optional_file/
load_and_fix_igblast_auxdata <- function(igblast_organism)
{
    igblast_organism <- normalize_igblast_organism(igblast_organism)
    auxdata <- load_auxdata(igblast_organism, which="original")
    if (igblast_organism == "human") {
        allele_names <- auxdata[ , "allele_name"]
        ## We know that the human auxdata provided by IgBLAST is incorrect
        ## for alleles IGHJ6*02 and IGHJ6*03, so is discordant with what
        ## install_germline_db() will compute. So we correct these rows.
        ## Note that if the user updated their "live" IgBLAST data with
        ## update_live_igdata(), then these rows should already be correct
        ## in 'auxdata'.
        fixme <- allele_names %in% c("IGHJ6*02", "IGHJ6*03")
        auxdata[fixme, "extra_bps"] <- 1L  # replace 0 with 1
        ## It's also pretty clear that the cdr3_end for TRDJ3*01 should be
        ## 24, not 21: the FWR4 starts at position 25, which is where the
        ## FGXG motif is found. So we fix that too.
        ## TODO: Report this to the IgBLAST folks at NCBI.
        fixme <- allele_names == "TRDJ3*01"
        auxdata[fixme, "cdr3_end"] <- 24L  # replace 21 with 24
    }
    if (igblast_organism == "mouse") {
        allele_names <- auxdata[ , "allele_name"]
        ## We know that the mouse auxdata provided by IgBLAST is incorrect
        ## for alleles TRAJ31*02, TRAJ32*02, TRAJ45*02, and TRAJ59*01,
        ## so is discordant with what install_germline_db() will compute.
        ## So we correct these rows.
        ## I reported the incongruent extra_bps for the 4 alleles below
        ## to the IgBLAST folks at NCBI on June 2, 2026.
        fixme <- allele_names == "TRAJ31*02"
        auxdata[fixme, "extra_bps"] <- 1L  # replace 0 with 1
        fixme <- allele_names == "TRAJ32*02"
        auxdata[fixme, "extra_bps"] <- 1L  # replace 0 with 1
        fixme <- allele_names == "TRAJ45*02"
        auxdata[fixme, "extra_bps"] <- 0L  # replace 1 with 0
        fixme <- allele_names == "TRAJ59*01"
        auxdata[fixme, "extra_bps"] <- 2L  # replace 1 with 2
        ## IGLJ2P*01 was only added to mouse_gl.aux recently (on June 16, 2026)
        ## by the NCBI folks. They annotated it consistently with our own
        ## igblastr-generated annotations here
        ## https://github.com/HyrienLab/igblastr/wiki/Auxiliary-data-in-igblastr
        ## except that they set "coding_frame_start" to 3 instead of 0.
        ## We fix that.
        fixme <- allele_names == "IGLJ2P*01"
        auxdata[fixme, "coding_frame_start"] <-
            auxdata[fixme, "coding_frame_start"] %% 3L
    }
    ## We also remove rows for which the "coding_frame_start" / "cdr3_end"
    ## combination doesn't make sense. This is the case for example for
    ## human allele TRAJ31*01 and for mouse alleles TRAJ21*02, TRAJ24*02,
    ## and TRDJ2*02 (as of May 31, 2026).
    ## I reported this to the IgBLAST folks at NCBI on June 2, 2026.
    coding_frame_start <- auxdata[ , "coding_frame_start"]
    cdr3_end <- auxdata[ , "cdr3_end"]
    drop_idx <- which(((cdr3_end + 1L - coding_frame_start) %% 3L) != 0L)
    if (length(drop_idx) != 0L)
        auxdata <- auxdata[-drop_idx, ]
    auxdata
}

