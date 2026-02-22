test_that("prepare_igblastn_cmdline_args()", {
    ## prepare_igblastn_cmdline_args() is not exported.
    prepare_igblastn_cmdline_args <- igblastr:::prepare_igblastn_cmdline_args
    CORE_ARGNAMES <- c("query", "outfmt",
                       paste0("germline_db_", c("V", "D", "J")),
                       "organism", "domain_system", "ig_seqtype", "out")
    organism_idx <- match("organism", CORE_ARGNAMES)
    use_c_region_db("")

    ## --- no germline db selection ---

    ## For this test, we use invalid V-, D-, J-region dbs but
    ## prepare_igblastn_cmdline_args() should still accept them.
    ## Should always work, even if no germline db is currently in use.
    db_path <- system.file(package="igblastr", "extdata")
    germline_db_V <- file.path(db_path, "V")
    germline_db_D <- file.path(db_path, "D")
    germline_db_J <- file.path(db_path, "J")
    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
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
    expect_error(prepare_igblastn_cmdline_args("path/to/query",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J),
                 regexp="organism.*must be specified")
    ## Error: 'ig_seqtype' must be specified
    expect_error(prepare_igblastn_cmdline_args("path/to/query",
                     germline_db_V=germline_db_V,
                     germline_db_D=germline_db_D,
                     germline_db_J=germline_db_J,
                     organism="rabbit"),
                 regexp="ig_seqtype.*must be specified")

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
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
    suppressMessages(use_germline_db(db_name))
    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
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

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              c_region_db=NULL,
                                              auxiliary_data=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    expect_identical(names(cmd_args), expected_argnames)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$custom_internal_data, get_intdata_path(db_name))
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "Ig")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              c_region_db=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    opt_argnames <- "auxiliary_data"
    expected_argnames <- append(expected_argnames, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$custom_internal_data, get_intdata_path(db_name))
    expect_identical(cmd_args$auxiliary_data, get_auxdata_path("human"))
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "Ig")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              auxiliary_data=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    expect_identical(names(cmd_args), expected_argnames)

    db_name <- "_IMGT.rabbit.IGH.202412"
    use_c_region_db(db_name)
    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              auxiliary_data=NULL)
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "c_region_db"
    opt_argnames <- "custom_internal_data"
    expected_argnames <- append(expected_argnames, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query")
    expect_true(is.list(cmd_args))
    expect_true(all(lengths(cmd_args) == 1L))
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "c_region_db"
    opt_argnames <- c("custom_internal_data", "auxiliary_data")
    expected_argnames <- append(expected_argnames, opt_argnames, organism_idx)
    expect_identical(names(cmd_args), expected_argnames)

    ## --- selecting IMGT-202531-1.Mus_musculus.IGH+IGK+IGL ---

    db_name <- install_IMGT_germline_db("202531-1", "Mus musculus",
                                        tcr.db=TRUE, force=TRUE)
    suppressMessages(use_germline_db(db_name))
    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              c_region_db=NULL,
                                              auxiliary_data=NULL)
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    expect_identical(names(cmd_args), expected_argnames)
    expect_identical(cmd_args$query, "path/to/query")
    expect_identical(cmd_args$outfmt, "AIRR")
    expect_identical(cmd_args$domain_system, "imgt")
    expect_identical(cmd_args$ig_seqtype, "TCR")
    expect_true(attr(cmd_args$out, "safe_to_remove"))

    ## --- selecting IMGT-202518-3.Sus_scrofa.IGH+IGK+IGL ---
    ## Note that pig is not an IgBLAST organism!

    db_name <- install_IMGT_germline_db("202518-3", "Sus_scrofa", force=TRUE)
    suppressMessages(use_germline_db(db_name))
    use_c_region_db("")

    regexp <- "what auxiliary data to use"
    expect_error(prepare_igblastn_cmdline_args("path/to/query", regexp=regexp))

    cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                              auxiliary_data=NULL)
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    expect_identical(names(cmd_args), expected_argnames)

    db_name <- install_IMGT_germline_db("202518-3", "Sus_scrofa",
                                        without.intdata=TRUE, force=TRUE)
    suppressMessages(use_germline_db(db_name))

    errmsg <- "Don't know how to infer 'organism' from germline db name"
    expect_error2(prepare_igblastn_cmdline_args("path/to/query"), errmsg)

    custom_internal_data <- get_intdata_path("_AIRR.human.IGH+IGK+IGL.202410")
    regexp <- "what auxiliary data to use"
    expect_error(
        suppressWarnings(
            prepare_igblastn_cmdline_args(
                    "path/to/query",
                    custom_internal_data=custom_internal_data)
        ),
        regexp=regexp
    )

    expect_warning(
        cmd_args <- prepare_igblastn_cmdline_args("path/to/query",
                                 custom_internal_data=custom_internal_data,
                                 auxiliary_data=NULL),
        regexp="Incomplete custom internal data"
    )
    expected_argnames <- CORE_ARGNAMES
    expected_argnames[[organism_idx]] <- "custom_internal_data"
    expect_identical(names(cmd_args), expected_argnames)
})

