This is the `extdata/germline_sequences/AIRR/human/202309/ref/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for Human. The sequences were extracted from the
following datasets:
- `IGH_VDJ` version 7 "Reference Set" (a.k.a. "extended")
- `IGKappa_VJ` version 2 "Reference Set" (a.k.a. "extended")
- `IGLambda_VJ` version 1 "Reference Set" (a.k.a. "extended")

IMPORTANT NOTE: These were the most current versions at the date indicated
by the name of the parent folder of this folder (202309 i.e. Sept 2023).

FWIW these datasets can be manually downloaded with:
```
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGH_VDJ/7/ungapped_ex >Homo_sapiens_IGH_VDJ_rev_7_ungapped_ex.fasta
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGKappa_VJ/2/ungapped_ex >Homo_sapiens_IGKappa_VJ_rev_2_ungapped_ex.fasta
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGLambda_VJ/1/ungapped_ex >Homo_sapiens_IGLambda_VJ_rev_1_ungapped_ex.fasta
```

However, the FASTA files in this folder were obtained programmatically
by running the following code in the folder on Oct 20, 2025:
```r
library(igblastr)
download_germline_sequences_from_OGRDB <- igblastr:::download_germline_sequences_from_OGRDB

Hfile_count <- download_germline_sequences_from_OGRDB("Human", set_name="IGH_VDJ", locus="IGH", release_version="7", extended=TRUE)
stopifnot(identical(Hfile_count, 3L))

Kfile_count <- download_germline_sequences_from_OGRDB("Human", set_name="IGKappa_VJ", locus="IGK", release_version="2", extended=TRUE)
stopifnot(identical(Kfile_count, 2L))

Jfile_count <- download_germline_sequences_from_OGRDB("Human", set_name="IGLambda_VJ", locus="IGL", release_version="1", extended=TRUE)
stopifnot(identical(Jfile_count, 2L))
```

