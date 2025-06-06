This is the extdata/ folder in the igblastr package.

This folder contains various data files used in the man page examples and/or
unit tests of the igblastr package.

Content
=======

- README.txt: This file.

- 1279067_1_Paired_sequences.fasta.gz: Compressed FASTA file containing 8,437
  pairs of human antibody sequences (16,874 individual sequences) downloaded
  from OAS (the Observed Antibody Space database). The file was obtained
  programmatically by running the following code in this folder on March 26,
  2025:

    library(igblastr)
    download_paired_OAS_units("Jaffe_2022", "1279067_1_Paired_All.csv.gz")
    df <- read_OAS_csv("1279067_1_Paired_All.csv.gz")
    sequences <- extract_sequences_from_paired_OAS_df(df, add.prefix=TRUE)
    writeXStringSet(sequences, "1279067_1_Paired_sequences.fasta.gz",
                    compress=TRUE)

- 1279067_1_Paired_All.json: JSON file containing the metadata associated
  with the sequences in 1279067_1_Paired_sequences.fasta.gz.
  This file was obtained by downloading 1279067_1_Paired_All.json directly
  from https://opig.stats.ox.ac.uk/webapps/ngsdb/paired/Jaffe_2022/json/

- heavy_sequences.fasta, light_sequences.fasta: Two FASTA files containing
  a small random set of heavy- and light-chain sequences extracted from
  1279067_1_Paired_sequences.fasta.gz, for use in the man page examples
  and unit tests of igblastr. Each file contains 125 sequences and
  the two files are paired, that is, the i-th heavy-chain sequence in
  heavy_sequences.fasta is paired with the i-th light-chain sequence in
  light_sequences.fasta. The two files were obtained programmatically with:

    library(igblastr)
    sequences <- readDNAStringSet("1279067_1_Paired_sequences.fasta.gz")
    num_pairs <- length(sequences) %/% 2L
    set.seed(2009)
    pair_selection <- sample(num_pairs, 125)
    heavy_sequences <- sequences[2L*pair_selection - 1L]
    light_sequences <- sequences[2L*pair_selection]
    writeXStringSet(heavy_sequences, "heavy_sequences.fasta")
    writeXStringSet(light_sequences, "light_sequences.fasta")

- constant_regions/: See README.txt in constant_regions/ subfolder.

- germline_sequences/: See README.txt in germline_sequences/ subfolder.

- ncbi_igblast_data_files/: See README.txt in ncbi_igblast_data_files/
  subfolder.

