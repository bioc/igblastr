This is the `extdata/germline_sets/OGRDB/human/202410/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) extracted from
the following AIRR-community/OGRDB germline sets for Human:
- `IGH_VDJ`:     version 9 (released on 2024-10-12)
- `IGKappa_VJ`:  version 4 (released on 2024-10-12)
- `IGLambda_VJ`: version 3 (released on 2024-10-12)

Note that these were the most current versions at the date indicated by
the name of this folder (202410 i.e. Oct 2024).

AIRR-community/OGRDB provides two flavors of each dataset, one called
the "Reference Set" and the other one called the "Source Set". The `ref/`
subfolder contains the former and `src/` the latter.

See `make_all_data.R` in the parent folder for how the data files in the
`ref/` and `src/` subfolders were created.

