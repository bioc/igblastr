This is the `extdata/germline_sets/OGRDB/rainbow_trout/202606/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for rainbow trout. The sequences are coming from the
following germline set:
- `IGH_VDJ` version 1

See release announcement for this germline set here:
https://wordpress.vdjbase.org/index.php/ogrdb/rainbow-trout-igh-set-now-available/

IMPORTANT NOTE: This was the most current version at the date indicated
by the name of the parent folder of this folder (202606 i.e. June 2026).

The FASTA files in this folder were obtained programmatically
by running the following code in the folder on July 8, 2026:
```r
library(igblastr)
germline_set <- c(IGH_VDJ=1)
download_OGRDB_germline_sequences("Oncorhynchus mykiss", germline_set)
```

