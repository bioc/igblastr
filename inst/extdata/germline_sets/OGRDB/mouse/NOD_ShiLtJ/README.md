This is the `extdata/germline_sets/OGRDB/mouse/NOD_ShiLtJ/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for mouse strain NOD/ShiLtJ. The sequences are coming
from the following germline sets:
- `NOD/ShiLtJ IGH` version 1
- `NOD/ShiLtJ IGKV` version 1
- `IGKJ (all strains)` version 1
- `NOD/ShiLtJ IGLV` version 1
- `IGLJ (all strains)` version 1

Note that these FASTA files were obtained programmatically by running the
following code in the folder on Dec 29, 2024:
```r
library(igblastr)
germline_sets <- c(`NOD/ShiLtJ IGH`=1,
                   `NOD/ShiLtJ IGKV`=1, `IGKJ (all strains)`=1,
                   `NOD/ShiLtJ IGLV`=1, `IGLJ (all strains)`=1)
download_OGRDB_germline_sequences("Mus musculus", germline_sets, gapped=FALSE)
```

In addition to the FASTA files, this folder has a `version` file that
contains the date of the download in YYYYMM format.

