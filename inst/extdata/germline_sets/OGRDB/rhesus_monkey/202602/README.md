This is the `extdata/germline_sets/OGRDB/rhesus_monkey/202602/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for rhesus monkey. The sequences are coming from the
following germline sets:
- `IGH_VDJ` version 2
- `IGK_VJ` version 2
- `IGL_VJ` version 2

See release announcement for these germline sets here:
https://wordpress.vdjbase.org/index.php/ogrdb\_news/macaque-sets-updated/

IMPORTANT NOTE: These were the most current versions at the date indicated
by the name of the parent folder of this folder (202602 i.e. Feb 2026).

The FASTA files in this folder were obtained programmatically
by running the following code in the folder on March 2, 2026:
```r
library(igblastr)
germline_sets <- c(IGH_VDJ=2, IGK_VJ=2, IGL_VJ=2)
download_OGRDB_germline_sequences("Macaca mulatta", germline_sets)
```

