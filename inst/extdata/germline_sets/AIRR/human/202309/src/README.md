This is the `extdata/germline_sets/AIRR/human/202309/src/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for Human. The sequences were extracted from the
following datasets:
- `IGH_VDJ` version 7 "Source Set"
- `IGKappa_VJ` version 2 "Source Set"
- `IGLambda_VJ` version 1 "Source Set"

IMPORTANT NOTE: These were the most current versions at the date indicated
by the name of the parent folder of this folder (202309 i.e. Sept 2023).

FWIW these datasets can be manually downloaded with:
```
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGH_VDJ/7/ungapped >Homo_sapiens_IGH_VDJ_rev_7_ungapped.fasta
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGKappa_VJ/2/ungapped >Homo_sapiens_IGKappa_VJ_rev_2_ungapped.fasta
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGLambda_VJ/1/ungapped >Homo_sapiens_IGLambda_VJ_rev_1_ungapped.fasta
```

However, the FASTA files in this folder were obtained programmatically
by running the following code in the folder on Oct 20, 2025:
```r
library(igblastr)
germline_sets <- c(IGH_VDJ=7, IGKappa_VJ=2, IGLambda_VJ=1)
download_OGRDB_germline_sequences("Homo sapiens", germline_sets, gapped=FALSE, source_set=TRUE)
```

