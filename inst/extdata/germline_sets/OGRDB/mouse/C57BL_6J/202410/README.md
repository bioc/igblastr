This is the `extdata/germline_sets/OGRDB/mouse/C57BL_6J/202410/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for mouse strain C57BL/6J. The sequences are coming
from the following germline sets:
- `C57BL/6 IGH` version 5
- `C57BL/6J IGKV` version 1
- `C57BL/6J IGLV` version 1
- `IGKJ (all strains)` version 1
- `IGLJ (all strains)` version 1

Note that these FASTA files were obtained programmatically by running the
following code in the folder on August 5, 2026:
```r
library(igblastr)
germline_sets <- c(`C57BL/6 IGH`=5, `C57BL/6J IGKV`=1, `C57BL/6J IGLV`=1,
                   `IGKJ (all strains)`=1, `IGLJ (all strains)`=1)
download_OGRDB_germline_sequences("Mus musculus", germline_sets)
```

