This is the `extdata/germline_sequences/AIRR/human/202410/ref/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for Human. The sequences were extracted from the
following datasets:
- `IGH_VDJ` version 9 "Reference Set" (a.k.a. "extended")
- `IGKappa_VJ` version 4 "Reference Set" (a.k.a. "extended")
- `IGLambda_VJ` version 3 "Reference Set" (a.k.a. "extended")

IMPORTANT NOTE: These were the most current versions at the date indicated
by the name of the parent folder of this folder (202410 i.e. Oct 2024).

FWIW these datasets can be manually downloaded with:
```
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGH_VDJ/9/ungapped_ex >Homo_sapiens_IGH_VDJ_rev_9_ungapped_ex.fasta
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGKappa_VJ/4/ungapped_ex >Homo_sapiens_IGKappa_VJ_rev_4_ungapped_ex.fasta
curl https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens/IGLambda_VJ/3/ungapped_ex >Homo_sapiens_IGLambda_VJ_rev_3_ungapped_ex.fasta
```

However, the FASTA files in this folder were obtained programmatically
by running the following code in the folder on Dec 29, 2024:
```r
library(igblastr)
download_germline_sequences_from_OGRDB <- igblastr:::download_germline_sequences_from_OGRDB

Hfile_count <- download_germline_sequences_from_OGRDB("Human", set_name="IGH_VDJ", locus="IGH", release_version="9", extended=TRUE)
stopifnot(identical(Hfile_count, 3L))

Kfile_count <- download_germline_sequences_from_OGRDB("Human", set_name="IGKappa_VJ", locus="IGK", release_version="4", extended=TRUE)
stopifnot(identical(Kfile_count, 2L))

Jfile_count <- download_germline_sequences_from_OGRDB("Human", set_name="IGLambda_VJ", locus="IGL", release_version="3", extended=TRUE)
stopifnot(identical(Jfile_count, 2L))
```

