test_that("use_germline_db()", {
    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
    use_germline_db(db_name)
    expect_identical(use_germline_db(), db_name)
})

test_that("load_germline_db()", {
    db_name <- "_OGRDB.human.IGH+IGK+IGL.202410"
    object <- load_germline_db(db_name)
    expect_true(is(object, "DNAStringSet"))
    expect_equal(length(object), 396)
    object <- load_germline_db(db_name, region_types="V")
    expect_true(is(object, "DNAStringSet"))
    expect_equal(length(object), 342)
    object <- load_germline_db(db_name, region_types="D")
    expect_true(is(object, "DNAStringSet"))
    expect_equal(length(object), 31)
    object <- load_germline_db(db_name, region_types="J")
    expect_true(is(object, "DNAStringSet"))
    expect_equal(length(object), 23)
    object <- load_germline_db(db_name, region_types="JV")
    expect_true(is(object, "DNAStringSet"))
    expect_equal(length(object), 365)
})

