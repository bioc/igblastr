This is the `extdata/germline_sets/AIRR/mouse/PWD_PhJ/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for mouse strain PWD/PhJ. The sequences are coming
from the following germline sets:
- `PWD/PhJ IGH` version 2
- `PWD/PhJ IGKV` version 1
- `IGKJ (all strains)` version 1
- `PWD/PhJ IGLV` version 1
- `IGLJ (all strains)` version 1

Note that these FASTA files were obtained programmatically by running the
following code in the folder on Dec 29, 2024:
```r
library(igblastr)
germline_sets <- c(`PWD/PhJ IGH`=2,
                   `PWD/PhJ IGKV`=1, `IGKJ (all strains)`=1,
                   `PWD/PhJ IGLV`=1, `IGLJ (all strains)`=1)
download_OGRDB_germline_sequences("Mus musculus", germline_sets, gapped=FALSE)
```

In addition to the FASTA files, this folder has a `version` file that
contains the date of the download in YYYYMM format.

