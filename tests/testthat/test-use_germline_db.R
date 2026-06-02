test_that("use_germline_db()", {
    db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"
    use_germline_db(db_name)
    expect_identical(use_germline_db(), db_name)
})

test_that("load_germline_sequences()", {
    db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"
    object <- load_germline_sequences(db_name)
    expect_true(is(object, "DNAStringSet"))
    expect_equal(length(object), 421)
    object <- load_germline_sequences(db_name, region_types="V")
    expect_true(is(object, "DNAStringSet"))
    expect_equal(length(object), 367)
    object <- load_germline_sequences(db_name, region_types="D")
    expect_true(is(object, "DNAStringSet"))
    expect_equal(length(object), 31)
    object <- load_germline_sequences(db_name, region_types="J")
    expect_true(is(object, "DNAStringSet"))
    expect_equal(length(object), 23)
    object <- load_germline_sequences(db_name, region_types="JV")
    expect_true(is(object, "DNAStringSet"))
    expect_equal(length(object), 390)
})

