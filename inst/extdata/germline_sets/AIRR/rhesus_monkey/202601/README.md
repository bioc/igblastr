This is the `extdata/germline_sets/AIRR/rhesus_monkey/202601/` folder
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

The folder also contains the corresponding V gene FWR/CDR annotation files
for the IMGT domain system (`*.ndm.imgt` files).

The FASTA files in this folder were obtained programmatically by running
the following code in the folder on January 8, 2026:
```r
library(igblastr)
germline_sets <- c(IGH_VDJ=1, IGK_VJ=1, IGL_VJ=1)
download_OGRDB_germline_sequences("Macaca mulatta", germline_sets, gapped=FALSE)
```

Files `*.ndm.imgt` were obtained programmatically by running the
following code in the folder on Jan 12, 2026:
```r
library(igblastr)
germline_sets <- c(IGH_VDJ=1, IGK_VJ=1, IGL_VJ=1)
igblastr:::download_V_ndm_data_from_OGRDB("Macaca mulatta", germline_sets)
for (group in c("IGHV", "IGKV", "IGLV")) {
    V_names1 <- names(fasta.seqlengths(paste0(group, ".fasta")))
    V_names2 <- read_V_ndm_data(paste0(group, ".ndm.imgt"))$allele_name
    stopifnot(length(V_names1) == length(V_names2),
              setequal(V_names1, V_names2))
}
```

