## Differences between _AIRR.human.IGH+IGK+IGL.202309 and
## _AIRR.human.IGH+IGK+IGL.202401:
##   o No alleles were added or removed.
##   o The 5 following alleles were renamed (but their sequences remained
##     the same):
##       IGHV3-13*i01  -> IGHV3-13*06
##       IGHV3-43D*i02 -> IGHV3-43D*05
##       IGKV1-33*01   -> IGKV1D-33*01
##       IGKV2D-28*01  -> IGKV2-28*01
##       IGKV6-21*01   -> IGKV6D-21*01
##   o Zero change in allele sequences.
##
## This can be seen with the code below.

library(igblastr)

## Load the two databases as DNAStringSet objects:
old <- load_germline_db("_AIRR.human.IGH+IGK+IGL.202309")
new <- load_germline_db("_AIRR.human.IGH+IGK+IGL.202401")

## 5 alleles were replaced:
old_names <- setdiff(names(old), names(new))
old_names
# [1] "IGHV3-13*i01"  "IGHV3-43D*i02" "IGKV1-33*01"   "IGKV2D-28*01"
# [5] "IGKV6-21*01"
new_names <- setdiff(names(new), names(old))
new_names
# [1] "IGHV3-13*06"  "IGHV3-43D*05" "IGKV1D-33*01" "IGKV2-28*01"  "IGKV6D-21*01"

## However, this replacement didn't alter the sequences:
old[old_names] == new[new_names]
# [1] TRUE TRUE TRUE TRUE TRUE

## so it's fair to call this a renaming rather than a replacement.

## No other sequence has changed:
shared_names <- intersect(names(old), names(new))
all(old[shared_names] == new[shared_names])
# [1] TRUE

## So overall, no allele sequence has changed.

