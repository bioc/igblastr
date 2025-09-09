This is the extdata/constant_regions/IMGT/ folder in the igblastr package.

This folder contains the C-region sequences (in FASTA format) obtained
from IMGT for the 5 official IgBLAST organisms: human, mouse, rat, rabbit,
and rhesus_monkey.

Note that these FASTA files were obtained programmatically by running the
following code in the folder on Sep 1, 2025:

    igblastr:::download_C_sequence_sets_from_IMGT()

About the 'version' files
-------------------------

The */IG/14.1/ and */TR/ subfolders have a 'version' file that we added
manually. It indicates the date of the download in YYYYMM format. This is
used by internal utility form_builtin_IMGT_c_region_db_name() to form the
names of the corresponding built-in IMGT C-region dbs.

Note that the date in IG/14.1/version is 202412 for human/rabbit and 202508
for rat. The reason for this is because the C-region sequences for these
organisms were actually downloaded from IMGT and included in the igblastr
package in December 2024 (for human/rabbit) and on Aug 21, 2025 for rat.
However, as of Sep 1, 2025, they've not changed, that is,
igblastr:::download_C_sequence_sets_from_IMGT() still retrieves
the same sequences.

