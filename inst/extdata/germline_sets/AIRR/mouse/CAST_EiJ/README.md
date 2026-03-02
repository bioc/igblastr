This is the `extdata/germline_sets/AIRR/mouse/CAST_EiJ/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for mouse strain CAST/EiJ. The sequences are coming
from the following germline sets:
- `CAST/EiJ IGH` version 1
- `CAST/EiJ IGKV` version 1
- `IGKJ (all strains)` version 1
- `CAST/EiJ IGLV` version 1
- `IGLJ (all strains)` version 1

Note that these FASTA files were obtained programmatically by running the
following code in the folder on Dec 29, 2024:
```r
library(igblastr)
germline_sets <- c(`CAST/EiJ IGH`=1,
                   `CAST/EiJ IGKV`=1, `IGKJ (all strains)`=1,
                   `CAST/EiJ IGLV`=1, `IGLJ (all strains)`=1)
download_OGRDB_germline_sequences("Mus musculus", germline_sets, gapped=FALSE)
```

In addition to the FASTA files, this folder has a `version` file that
contains the date of the download in YYYYMM format.

