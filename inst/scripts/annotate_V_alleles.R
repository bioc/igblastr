## On 2026-01-07, I added experimental functions annotate_heavy_V_alleles()
## and annotate_light_V_alleles() to generate the CDR/FWR segmentation
## information of a set of V allele sequences. Note that this information
## is included in the internal data shipped with IgBLAST but only for about
## half of all known human V alleles!


## --------------------------------------------------------------------------
## TYPICAL USE
##

## Let's say we want to annotate the V alleles included in germline db
## IMGT-202531-1.Homo_sapiens.IGH+IGK+IGL:
library(igblastr)
db_name <- "IMGT-202531-1.Homo_sapiens.IGH+IGK+IGL"
if (!(db_name %in% list_germline_dbs(names.only=TRUE)))
    install_IMGT_germline_db("202531-1", "Homo_sapiens")
all_V_alleles <- load_germline_db(db_name, region_types="V")
all_V_alleles  # 705 alleles

## Let's split them in 3 groups, one per IG locus:
has_prefix <- igblastr:::has_prefix
IGHV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGH")]
IGKV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGK")]
IGLV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGL")]

##
## --- Annotate all human heavy V alleles ---
##

## Functions annotate_heavy_V_alleles()/annotate_light_V_alleles() actually
## use igblastn() to analyze/annotate the V alleles. igblastn() will use
## _OGRDB.human.IGH+IGK+IGL.202410 for the annalysis:
db_name0 <- "_OGRDB.human.IGH+IGK+IGL.202410"
system.time(IGHV_ann <- annotate_heavy_V_alleles(IGHV_alleles, db_name0))
#     user   system  elapsed
# 1953.872    4.628  491.682

## Does 'IGHV_ann' agree with the V-allele annotations shipped with IgBLAST?
check_ann <- function(ann, intdata0, strict=FALSE) {
    m <- match(ann[ , "allele_name"], intdata0[ , "allele_name"])
    ann1 <- ann[!is.na(m), ]
    perc <- format(100 * nrow(ann1) / nrow(ann), digits=2)
    cat(nrow(ann1), "/", nrow(ann), " (", perc, "%) of the ",
        "V alleles in 'ann' are already annotated in 'intdata0' ",
        "(IgBLAST internal data).\n", sep="")
    ann2 <- intdata0[m[!is.na(m)], colnames(ann)]
    rownames(ann1) <- rownames(ann2) <- NULL
    if (strict) {
        ok <- identical(ann1, ann2)
    } else {
        ok <- identical(ann1[ , -2L], ann2[ , -2L])
    }
    if (ok) {
        msg <- "For all of them, 'ann' and 'intdata0' are in agreement"
        if (!strict) {
            msg <- c(msg, " (modulo the \"fwr1_start\" column in 'ann' ",
                     "which can contain values != 1)")
        }
    } else {
        msg <- "For some of them, 'ann' and 'intdata0' disagree"
    }
    cat(msg, ".\n", sep="")
    ok
}

## V-allele annotations shipped with IgBLAST:
human_intdata0 <- load_intdata("human")
check_ann(IGHV_ann, human_intdata0, strict=TRUE)  # FALSE!
stopifnot(check_ann(IGHV_ann, human_intdata0))

##
## --- Annotate all human light V alleles ---
##

system.time(IGKV_ann <- annotate_light_V_alleles(IGKV_alleles, db_name0, "IGK"))
#    user  system elapsed
#  14.311   0.063   3.934

## Does 'IGKV_ann' agree with the V-allele annotations shipped with IgBLAST?
stopifnot(check_ann(IGKV_ann, human_intdata0, strict=TRUE))

system.time(IGLV_ann <- annotate_light_V_alleles(IGLV_alleles, db_name0, "IGL"))
#    user  system elapsed
#  19.399   0.055   5.179

## Does 'IGLV_ann' agree with the V-allele annotations shipped with IgBLAST?
stopifnot(check_ann(IGLV_ann, human_intdata0, strict=TRUE))


## --------------------------------------------------------------------------
## REPEATING THE ABOVE ANALYSIS FOR mouse, rat, rhesus monkey, and rabbit
##

## --- mouse ---

db_name <- "IMGT-202531-1.Mus_musculus.IGH+IGK+IGL"
if (!(db_name %in% list_germline_dbs(names.only=TRUE)))
    install_IMGT_germline_db("202531-1", "Mus_musculus")
all_V_alleles <- load_germline_db(db_name, region_types="V")
all_V_alleles  # 865 alleles

IGHV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGH")]
IGKV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGK")]
IGLV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGL")]

db_name0 <- "_OGRDB.mouse.PWD_PhJ.IGH+IGK+IGL.202501"
system.time(IGHV_ann <- annotate_heavy_V_alleles(IGHV_alleles, db_name0))
#    user  system elapsed
# 430.251   0.903 108.513

system.time(IGKV_ann <- annotate_light_V_alleles(IGKV_alleles, db_name0, "IGK"))
#    user  system elapsed
#  27.366   0.076   7.158

system.time(IGLV_ann <- annotate_light_V_alleles(IGLV_alleles, db_name0, "IGL"))
#    user  system elapsed
#   2.157   0.019   0.739

mouse_ann <- rbind(IGHV_ann, IGKV_ann, IGLV_ann)
## Compare with V-allele annotations shipped with IgBLAST:
mouse_intdata0 <- load_intdata("mouse")
## Note that only one allele in 'all_V_alleles' is already annotated in
## IgBLAST's internal data!
intersect(names(all_V_alleles), mouse_intdata0$allele_name)  # "IGKV20-101-2*01"
stopifnot(check_ann(mouse_ann, mouse_intdata0, strict=TRUE))

## --- rat ---

db_name <- "IMGT-202531-1.Rattus_norvegicus.IGH+IGK+IGL"
if (!(db_name %in% list_germline_dbs(names.only=TRUE)))
    install_IMGT_germline_db("202531-1", "Rattus_norvegicus")
all_V_alleles <- load_germline_db(db_name, region_types="V")
all_V_alleles  # 403 alleles
IGHV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGH")]
IGKV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGK")]
IGLV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGL")]

db_name0 <- db_name
system.time(IGHV_ann <- annotate_heavy_V_alleles(IGHV_alleles, db_name0))
#    user  system elapsed
# 585.088   1.728 147.804

system.time(IGKV_ann <- annotate_light_V_alleles(IGKV_alleles, db_name0, "IGK"))
#    user  system elapsed
#  19.918   0.053   5.241

system.time(IGLV_ann <- annotate_light_V_alleles(IGLV_alleles, db_name0, "IGL"))
#    user  system elapsed
#   0.933   0.012   0.580

rat_ann <- rbind(IGHV_ann, IGKV_ann, IGLV_ann)
## Compare with V-allele annotations shipped with IgBLAST:
rat_intdata0 <- load_intdata("rat")
stopifnot(check_ann(rat_ann, rat_intdata0, strict=TRUE))

## --- rhesus monkey ---

db_name <- "IMGT-202531-1.Macaca_mulatta.IGH+IGK+IGL"
if (!(db_name %in% list_germline_dbs(names.only=TRUE)))
    install_IMGT_germline_db("202531-1", "Macaca_mulatta")
all_V_alleles <- load_germline_db(db_name, region_types="V")
all_V_alleles  # 457 alleles
IGHV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGH")]
IGKV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGK")]
IGLV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGL")]

db_name0 <- db_name
system.time(IGHV_ann <- annotate_heavy_V_alleles(IGHV_alleles, db_name0))
#     user   system  elapsed
# 1446.023    4.330  365.055

system.time(IGKV_ann <- annotate_light_V_alleles(IGKV_alleles, db_name0, "IGK"))
#   user  system elapsed
# 15.519   0.059   4.132

system.time(IGLV_ann <- annotate_light_V_alleles(IGLV_alleles, db_name0, "IGL"))
#   user  system elapsed
# 24.178   0.069   6.338

rhesus_monkey_ann <- rbind(IGHV_ann, IGKV_ann, IGLV_ann)
## Compare with V-allele annotations shipped with IgBLAST:
rhesus_monkey_intdata0 <- load_intdata("rhesus_monkey")
stopifnot(check_ann(rhesus_monkey_ann, rhesus_monkey_intdata0, strict=TRUE))

## --- rabbit ---

db_name <- "IMGT-202531-1.Oryctolagus_cuniculus.IGH+IGK+IGL"
if (!(db_name %in% list_germline_dbs(names.only=TRUE)))
    install_IMGT_germline_db("202531-1", "Oryctolagus_cuniculus")
all_V_alleles <- load_germline_db(db_name, region_types="V")
all_V_alleles  # 148 alleles
IGHV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGH")]
IGKV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGK")]
IGLV_alleles <- all_V_alleles[has_prefix(names(all_V_alleles), "IGL")]

db_name0 <- db_name
system.time(IGHV_ann <- annotate_heavy_V_alleles(IGHV_alleles, db_name0))
#   user  system elapsed
# 75.613   0.331  19.393

system.time(IGKV_ann <- annotate_light_V_alleles(IGKV_alleles, db_name0, "IGK"))
#   user  system elapsed
# 27.088   0.070   7.044

system.time(IGLV_ann <- annotate_light_V_alleles(IGLV_alleles, db_name0, "IGL"))
#   user  system elapsed
#  2.262   0.022   0.739

rabbit_ann <- rbind(IGHV_ann, IGKV_ann, IGLV_ann)
## Compare with V-allele annotations shipped with IgBLAST:
rabbit_intdata0 <- load_intdata("rabbit")
stopifnot(check_ann(rabbit_ann, rabbit_intdata0, strict=TRUE))


## --------------------------------------------------------------------------
## CONCLUSION
##

## annotate_heavy_V_alleles() and annotate_light_V_alleles() produce
## annotations that are **always** in agreement with IgBLAST's internal data!

