test_that("use_c_region_db()", {
    use_c_region_db("")
    expect_identical(use_c_region_db(), "")

    db_name <- "_IMGT.rabbit.IGH.202412"
    use_c_region_db(db_name)
    expect_identical(use_c_region_db(), db_name)
})

test_that("load_c_region_db()", {
    db_name <- "_IMGT.rabbit.IGH.202412"
    object <- load_c_region_db(db_name)
    expect_true(is(object, "DNAStringSet"))
    expect_equal(length(object), 28)
})

