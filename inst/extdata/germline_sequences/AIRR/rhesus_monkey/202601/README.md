This is the `extdata/germline_sequences/AIRR/rhesus_monkey/202601/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for rhesus monkey. The sequences are coming from the
following germline sets:
- `IGH_VDJ` version 1
- `IGK_VJ` version 1
- `IGL_VJ` version 1

See release announcement for these germline sets here:
https://wordpress.vdjbase.org/index.php/ogrdb/rhesus-macaque-ig-germline-sets-released/

IMPORTANT NOTE: These were the most current versions at the date indicated
by the name of the parent folder of this folder (202601 i.e. Jan 2026).

The FASTA files in this folder were obtained programmatically by running
the following code in the folder on January 8, 2026:
```r
library(Biostrings)
BASE_URL <- "https://ogrdb.airr-community.org/download_germline_set/"
ORGANISM <- "Macaca%20mulatta"
germline_sets <- c(IGH="IGH_VDJ", IGK="IGK_VJ", IGL="IGL_VJ")

for (i in seq_along(germline_sets)) {
    set <- germline_sets[[i]]
    locus <- names(germline_sets)[[i]]
    url <- paste0(BASE_URL, ORGANISM, "/", set, "/1/ungapped")
    tmp_fasta_file <- tempfile()
    download.file(url, tmp_fasta_file)
    alleles <- readDNAStringSet(tmp_fasta_file)
    groups <- substr(names(alleles), 1, 4)
    unique_groups <- unique(groups)
    stopifnot(all(substr(unique_groups, 1, 3) == locus))
    for (group in unique_groups) {
        idx <- which(groups == group)
        fasta_file <- paste0(group, ".fasta")
        message("Writing ", fasta_file, " (", length(idx), " sequences) ... ",
                appendLF=FALSE)
        writeXStringSet(alleles[idx], fasta_file)
        message("ok")
    }
}
```

