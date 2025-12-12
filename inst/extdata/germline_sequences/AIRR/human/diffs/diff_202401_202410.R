## Differences between _AIRR.human.IGH+IGK+IGL.202401 and
## _AIRR.human.IGH+IGK+IGL.202410:
##   o No alleles were added or removed.
##   o The 2 following alleles were renamed (but their sequences remained
##     the same):
##       IGKV1D-39*01 -> IGKV1-39*01
##       IGKV2-40*01 -> IGKV2D-40*01
##   o Zero change in allele sequences.
##
## This can be seen with the code below.

library(igblastr)

## Load the two databases as DNAStringSet objects:
old <- load_germline_db("_AIRR.human.IGH+IGK+IGL.202401")
new <- load_germline_db("_AIRR.human.IGH+IGK+IGL.202410")

## 2 alleles were replaced:
old_names <- setdiff(names(old), names(new))
old_names
# [1] "IGKV1D-39*01" "IGKV2-40*01"
new_names <- setdiff(names(new), names(old))
new_names
# [1] "IGKV1-39*01"  "IGKV2D-40*01"

## However, this replacement didn't alter the sequences:
old[old_names] == new[new_names]
# [1] TRUE TRUE

## so it's fair to call this a renaming rather than a replacement.

## No other sequence has changed:
shared_names <- intersect(names(old), names(new))
all(old[shared_names] == new[shared_names])
# [1] TRUE

## So overall, no allele sequence has changed.

