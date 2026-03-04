This is the `extdata/germline_sets/OGRDB/human/202401/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) extracted from
the following AIRR-community/OGRDB germline sets for Human:
- `IGH_VDJ`:     version 8 (released on 2024-01-11)
- `IGKappa_VJ`:  version 3 (released on 2024-01-11)
- `IGLambda_VJ`: version 2 (released on 2024-01-11)

Note that these were the most current versions at the date indicated by
the name of this folder (202401 i.e. Jan 2024).

AIRR-community/OGRDB provides two flavors of each dataset, one called
the "Reference Set" and the other one called the "Source Set". The `ref/`
subfolder contains the former and `src/` the latter.

See `download_human_germline_sequences.R` in the parent folder for how
the FASTA files in the `ref/` and `src/` subfolders were created.

