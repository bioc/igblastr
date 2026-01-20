This is the `extdata/germline_sets/AIRR/human/202401/src/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for Human. The sequences were extracted from the
following datasets:
- `IGH_VDJ` version 8 "Source Set"
- `IGKappa_VJ` version 3 "Source Set"
- `IGLambda_VJ` version 2 "Source Set"

IMPORTANT NOTE: These were the most current versions at the date indicated
by the name of the parent folder of this folder (202401 i.e. Jan 2024).

FWIW these datasets can be manually downloaded with:
```
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGH_VDJ/9/ungapped >Homo_sapiens_IGH_VDJ_rev_8_ungapped.fasta
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGKappa_VJ/4/ungapped >Homo_sapiens_IGKappa_VJ_rev_3_ungapped.fasta
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGLambda_VJ/3/ungapped >Homo_sapiens_IGLambda_VJ_rev_2_ungapped.fasta
```

However, the FASTA files in this folder were obtained programmatically
by running the following code in the folder on Nov 17, 2025:
```r
library(igblastr)
download_germline_sequences_from_OGRDB <- igblastr:::download_germline_sequences_from_OGRDB

Hfile_count <- download_germline_sequences_from_OGRDB("Human", set_name="IGH_VDJ", locus="IGH", release_version="8")
stopifnot(identical(Hfile_count, 3L))

Kfile_count <- download_germline_sequences_from_OGRDB("Human", set_name="IGKappa_VJ", locus="IGK", release_version="3")
stopifnot(identical(Kfile_count, 2L))

Jfile_count <- download_germline_sequences_from_OGRDB("Human", set_name="IGLambda_VJ", locus="IGL", release_version="2")
stopifnot(identical(Jfile_count, 2L))
```

