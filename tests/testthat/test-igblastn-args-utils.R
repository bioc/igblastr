test_that("make_igblastn_command_line_args()", {
    ## make_igblastn_command_line_args() is not exported.
    make_igblastn_command_line_args <-
        igblastr:::make_igblastn_command_line_args
    CORE_ARGNAMES <- c("query", "outfmt",
                       paste0("germline_db_", c("V", "D", "J")),
                       "organism", "domain_system", "ig_seqtype", "out")
    organism_idx <- match("organism", CORE_ARGNAMES)

    ## --- no germline db selection ---

    ## For this test, we use invalid V-, D-, J-region dbs but
    ## make_igblastn_command_line_args() should still accept them.
    ## Should always work, even if no germline db is currently in use.
    db_path <- system.file(package="igblastr", "extdata")
    germline_db_V <- file.path(db_path, "V")
    germline_db_D <- file.path(db_path, "D")
    germline_db_J <- file.path(db_path, "J")
    cmd_args <- make_igblastn_command_line_args("path/to/query",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J,
                     organism="rhesus_monkey",
                     c_region_db=NULL,
                     custom_internal_data=NULL,
                     auxiliary_data=NULL,
                     ig_seqtype="Ig")
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expect_identical(names(cmd_args), CORE_ARGNAMES)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$organism, "rhesus_monkey")
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "Ig")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    ## Error: 'organism' must be specified
    expect_error(make_igblastn_command_line_args("path/to/query",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J),
                 regexp="organism.*must be specified")
    ## Error: 'ig_seqtype' must be specified
    expect_error(make_igblastn_command_line_args("path/to/query",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J,
                     organism="rabbit"),
                 regexp="ig_seqtype.*must be specified")

    cmd_args <- make_igblastn_command_line_args("path/to/query",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J,
                     organism="rabbit",
                     ig_seqtype="tcr")
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    opt_argnames <- "auxiliary_data"
    expected_argnames <- append(CORE_ARGNAMES, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)

    ## --- selecting _AIRR.human.IGH+IGK+IGL.202410 ---

    ## For these tests we're not specifying any of the 'germline_db_[VDJ]'
    ## arguments so we need to select a cached germline db. Once we do
    ## this, we don't need to specify 'organism' either.

    db_name <- "_AIRR.human.IGH+IGK+IGL.202410"
    use_germline_db(db_name)
    cmd_args <- make_igblastn_command_line_args("path/to/query",
                                                c_region_db=NULL,
                                                custom_internal_data=NULL,
                                                auxiliary_data=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expect_identical(names(cmd_args), CORE_ARGNAMES)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$organism, "human")
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "Ig")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    cmd_args <- make_igblastn_command_line_args("path/to/query",
                                                c_region_db=NULL,
                                                auxiliary_data=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    opt_argnames <- "custom_internal_data"
    expected_argnames <- append(CORE_ARGNAMES, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$organism, "human")
    expect_identical(cmd_args$custom_internal_data, get_intdata_path(db_name))
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "Ig")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    cmd_args <- make_igblastn_command_line_args("path/to/query",
                                                c_region_db=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    opt_argnames <- c("custom_internal_data", "auxiliary_data")
    expected_argnames <- append(CORE_ARGNAMES, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$organism, "human")
    expect_identical(cmd_args$custom_internal_data, get_intdata_path(db_name))
    expect_identical(cmd_args$auxiliary_data, get_auxdata_path("human"))
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "Ig")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    use_c_region_db("")
    cmd_args <- make_igblastn_command_line_args("path/to/query",
                                                auxiliary_data=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    opt_argnames <- "custom_internal_data"
    expected_argnames <- append(CORE_ARGNAMES, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)

    db_name <- "_IMGT.rabbit.IGH.202412"
    use_c_region_db(db_name)
    cmd_args <- make_igblastn_command_line_args("path/to/query",
                                                auxiliary_data=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    opt_argnames <- c("c_region_db", "custom_internal_data")
    expected_argnames <- append(CORE_ARGNAMES, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)

    cmd_args <- make_igblastn_command_line_args("path/to/query")
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    opt_argnames <- c("c_region_db", "custom_internal_data", "auxiliary_data")
    expected_argnames <- append(CORE_ARGNAMES, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)

    ## --- selecting IMGT-202531-1.Mus_musculus.IGH+IGK+IGL ---

    db_name <- install_IMGT_germline_db("202531-1", "Mus musculus",
                                        tcr.db=TRUE, force=TRUE)
    db_name <- use_germline_db(db_name)
    cmd_args <- make_igblastn_command_line_args("path/to/query",
                                                c_region_db=NULL,
                                                auxiliary_data=NULL)
    expect_identical(names(cmd_args), CORE_ARGNAMES)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$organism, "mouse")
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "TCR")
    expect_true(attr(cmd_args$out, "safe_to_remove"))
})

