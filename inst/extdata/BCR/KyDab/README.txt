This is the extdata/BCR/KyDab/ folder in the igblastr package.

This folder contains Kymouse BCR sequences downloaded from KyDab (Kymouse
Antibody Database): https://kydab.naturalantibody.com/

Content
=======

- README.txt: This file.

- Influenza-KV7015_heavy.fasta, Influenza-KV7015_light.fasta: Two FASTA files
  containing the Kymouse heavy- and light-chain BCR sequences from the
  KyDab Influenza #1 Project: https://kydab.naturalantibody.com/influenza-01
  The two files were obtained programmatically by running the following code
  in this folder on July 14, 2026:

    library(Biostrings)

    zipfile <- "Influenza-KV7015.zip"  # source dataset .zip
    zip_url <- paste0("https://static-ag-specific-data.naturalantibody.com/assets/data/influenza-kv7015/", zipfile)
    download.file(zip_url, zipfile)
    unzip(zipfile)

    sequences <- read.table("Influenza-KV7015.csv", header=TRUE, sep=",")
    heavy_sequences <- sequences$HEAVY_CHAIN_SEQ
    names(heavy_sequences) <- add_prefix(sequences$CLONE_NAME, "heavy_")
    light_sequences <- sequences$LIGHT_CHAIN_SEQ
    names(light_sequences) <- add_prefix(sequences$CLONE_NAME, "light_")

    heavy_sequences <- DNAStringSet(heavy_sequences[nzchar(heavy_sequences)])
    light_sequences <- DNAStringSet(light_sequences[nzchar(light_sequences)])

    writeXStringSet(heavy_sequences, "Influenza-KV7015_heavy.fasta")
    writeXStringSet(light_sequences, "Influenza-KV7015_light.fasta")

  Influenza-KV7015_heavy.fasta contains 3950 BCR nucleotide sequences.
  Influenza-KV7015_light.fasta contains 3852 BCR nucleotide sequences.

- Influenza-KV7015-metadata.json: Metadata associated with KyDab Influenza #1
  Project. The file is included in the source dataset Influenza-KV7015.zip
  downloaded and unzipped above. 

