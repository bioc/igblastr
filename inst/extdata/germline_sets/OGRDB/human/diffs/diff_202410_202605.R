## Differences between _OGRDB.human.IGH+IGK+IGL.202410 and
## _OGRDB.human.IGH+IGK+IGL.202605:
##   o 25 alleles were added, all V alleles.
##   o No allele sequence was touched.
##
## This can be seen with the code below.

library(igblastr)

## Load the two databases as DNAStringSet objects:
old <- load_germline_db("_OGRDB.human.IGH+IGK+IGL.202410")
new <- load_germline_db("_OGRDB.human.IGH+IGK+IGL.202605")

length(old)
# [1] 396

length(new)
# [1] 421

## All the alleles in the previous db version are in the new one:
all(names(old) %in% names(new))
# [1] TRUE

## and their sequences are still the same:
all(old == new[names(old)])
# [1] TRUE

## See 25 new alleles (they're all V alleles):
setdiff(names(new), names(old))

