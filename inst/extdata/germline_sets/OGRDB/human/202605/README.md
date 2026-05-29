This is the `extdata/germline_sets/OGRDB/human/202605/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) extracted from
the following AIRR-community/OGRDB germline sets for Human:
- `IGH_VDJ`:     version 10 (released on 2026-05-27)
- `IGKappa_VJ`:  version  5 (released on 2026-05-27)
- `IGLambda_VJ`: version  4 (released on 2026-05-27)

Note that these were the most current versions at the date indicated by
the name of this folder (202605 i.e. May 2026).

AIRR-community/OGRDB provides two flavors of each dataset, one called
the "Reference Set" and the other one called the "Source Set". The `ref/`
subfolder contains the former and `src/` the latter.

See `download_human_germline_sequences.R` in the parent folder for how
the FASTA files in the `ref/` and `src/` subfolders were created.

