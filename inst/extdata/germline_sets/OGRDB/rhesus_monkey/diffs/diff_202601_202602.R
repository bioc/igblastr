## Differences between _OGRDB.rhesus_monkey.IGH+IGK+IGL.202601 and
## _OGRDB.rhesus_monkey.IGH+IGK+IGL.202602:
##   o The two germline dbs contain exactly the same set of alleles with
##     the exact same ungapped sequences.
##   o Only difference is in the intdata.

library(igblastr)
intdata1 <- load_intdata("_OGRDB.rhesus_monkey.IGH+IGK+IGL.202601")
intdata2 <- load_intdata("_OGRDB.rhesus_monkey.IGH+IGK+IGL.202602")

### 187/2294 V alleles are annotated differently in the 2 dbs:
idx <- c(1872:1882, 1999:2087, 2121:2136, 2139:2149, 2158:2159, 2164:2192,
         2198:2201, 2212:2215, 2222:2226, 2256:2261, 2265:2272, 2284:2285)

identical(intdata1[-idx, ], intdata2[-idx, ])  # TRUE

## These 187 V alleles can be divided in 3 groups:
## - 11 alleles have different fwr2_end/cdr2_start;
## - 174 alleles have different cdr2_end/fwr3_start;
## - 2 alleles have different fwr3_end: IGLV25-QYEU*01 and IGLV25-QYEU*02

intdata1[2158:2159, c("allele_name", "fwr3_start", "fwr3_end")]
#         allele_name fwr3_start fwr3_end
# 2158 IGLV25-QYEU*01        157      267
# 2159 IGLV25-QYEU*02        157      267
intdata2[2158:2159, c("allele_name", "fwr3_start", "fwr3_end")]
#         allele_name fwr3_start fwr3_end
# 2158 IGLV25-QYEU*01        157      264
# 2159 IGLV25-QYEU*02        157      264

## So not exactly true that the Feb 2026 update of the macaque sets
## did not affect the delineation of the start of CDR3 as claimed here:
## https://wordpress.vdjbase.org/index.php/ogrdb_news/macaque-sets-updated/

