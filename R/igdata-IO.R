### =========================================================================
### Low-level helpers to read/write IgBLAST annotation files
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###
### IMPORTANT NOTES:
###
### 1) IgBLAST annotation files (i.e. internal and auxiliary data files) are
###    supposedly "tab-delimited" but they are broken in many ways:
###      o they use a mix of whitespaces for the field separators;
###      o each line contains a variable number of fields;
###      o some lines contain trailing whitespaces.
###    So we cannot simply read them with read.table().
###
### 2) There's no documented way to represent missing values in these
###    files. However, some digging into IgBLAST C++ code reveals that
###    the following special values are used to represent missing values:
###      o "N/A" for missing character values;
###      o -1 for missing integer values.
###    IgBLAST annotation files only contain character and integer data so
###    that's all we need to know.
###    See lines 148-184 in include/algo/blast/igblast/igblast.hpp for the
###    details:
###
###      https://github.com/ncbi/ncbi-cxx-toolkit-public/blob/main/include/algo/blast/igblast/igblast.hpp
###
###    This means that write_igdata() and read_igdata() must take care of
###    mapping back and forth between NAs to these special values.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### right_pad_and_unlist()
### matrix2df()
###
### Two very low-level helpers used by read_igdata() below in this file and
### other functions in other files.
###

### 'x' must be a list of character vectors of variable length.
### Conceptually right-pads the list elements with the specified padding
### string to make the list "constant-width" before unlisting it.
### Returns a character vector of length 'length(x) * width'.
right_pad_and_unlist <- function(x, padding_string, width=NA)
{
    stopifnot(is.list(x), isSingleStringOrNA(padding_string),
              isSingleNumberOrNA(width))
    x_len <- length(x)
    if (x_len == 0L)
        return(character(0))
    x_lens <- lengths(x)
    max_x_lens <- max(x_lens)
    if (is.na(width)) {
        width <- max_x_lens
    } else {
        width <- as.integer(width)
        stopifnot(width >= max_x_lens)
    }
    y_lens <- width - x_lens
    x_seqalong <- seq_along(x)
    f <- rep.int(x_seqalong, y_lens)
    attributes(f) <- list(levels=as.character(x_seqalong), class="factor")
    y <- split(rep.int(padding_string, length(f)), f)
    collate_subscript <- rep(x_seqalong, each=2L)
    collate_subscript[2L * x_seqalong] <- x_seqalong + x_len
    unlist(c(x, y)[collate_subscript], recursive=FALSE, use.names=FALSE)
}

matrix2df <- function(m, col2class)
{
    stopifnot(is.matrix(m), is.character(m),
              is.character(col2class), !is.null(names(col2class)),
              ncol(m) == length(col2class))
    cols <- lapply(setNames(seq_along(col2class), names(col2class)),
        function(i) {
            col <- m[ , i]  # character vector
            Class <- col2class[[i]]
            if (Class %in% c("integer", "numeric", "double"))
                col[col %in% "NA"] <- NA_character_
            as(col, Class)
        }
    )
    as.data.frame(cols, optional=TRUE, fix.empty.names=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .replace_IgBLAST_special_values_with_NAs(df)
### .replace_NAs_with_IgBLAST_special_values(df)
###

.IGBLAST_MISSING_STRING <- "N/A"
.IGBLAST_MISSING_INT <- -1L

.get_IgBLAST_missing_value <- function(Class)
{
    stopifnot(isSingleNonWhiteString(Class))
    switch(Class,
        character=.IGBLAST_MISSING_STRING,
        integer=.IGBLAST_MISSING_INT,
        stop(wmsg("data.frame can only have character or integer columns")))
}

### Empty strings and strings with whitespaces are not allowed in the
### character columns of a data.frame representing IgBLAST annotations.
.invalid_character_col <- function(col)
    !all(nzchar(col)) || any(has_whitespace(col))

### Negative values are not allowed in the integer columns of a data.frame
### representing IgBLAST annotations.
.invalid_integer_col <- function(col) any(col < 0L, na.rm=TRUE)

### NAs, the "N/A" string, the empty string, and strings with whitespaces
### are not allowed in the "allele_name" column (1st column).
.check_igdata_first_col <- function(df)
{
    stopifnot(is.data.frame(df), ncol(df) >= 2L,
              identical(colnames(df)[[1L]], "allele_name"))
    col1 <- df[[1L]]
    stopifnot(identical(class(col1)[[1L]], "character"))
    error <- anyNA(col1) || .IGBLAST_MISSING_STRING %in% col1 ||
             .invalid_character_col(col1)
    if (error)
        stop(wmsg("NAs, the \"", .IGBLAST_MISSING_STRING, "\" string, ",
                  "the empty string, and strings with whitespaces ",
                  "are not allowed in the \"allele_name\" column"))
}

.replace_IgBLAST_special_values_with_NAs <- function(df)
{
    .check_igdata_first_col(df)
    for (j in 2:ncol(df)) {
        col <- df[[j]]
        missing_val <- .get_IgBLAST_missing_value(class(col)[[1L]])
        col[col %in% missing_val] <- NA
        df[[j]] <- col
    }
    df
}

### Raises an error if 'df' contains empty strings or strings with
### whitespaces or negative values.
### Returns a data.frame that is guaranteed to NOT contain any NAs or empty
### strings or strings with whitespaces or negative values.
.replace_NAs_with_IgBLAST_special_values <- function(df)
{
    .check_igdata_first_col(df)
    for (j in 2:ncol(df)) {
        col <- df[[j]]
        ## Note that this must be checked **before** the replacement of
        ## NAs with 'missing_val'.
        if (is.integer(col) && .invalid_integer_col(col))
            stop(wmsg("integer columns are not allowed to contain ",
                      "negative values"))
        missing_val <- .get_IgBLAST_missing_value(class(col)[[1L]])
        col[is.na(col)] <- missing_val
        ## Note that it doesn't matter whether we check this before or after
        ## the replacement of NAs with 'missing_val'.
        if (is.character(col) && .invalid_character_col(col))
            stop(wmsg("character columns are not allowed to contain ",
                      "the empty string or strings with whitespaces"))
        df[[j]] <- col
    }
    df
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### read_igdata()
###

### Returns the table in a character matrix.
.read_jagged_table_as_character_matrix <- function(filepath)
{
    lines <- readLines(filepath)
    lines <- lines[nzchar(lines) & !has_prefix(lines, "#")]
    data <- strsplit(lines, split="[ \t]+")
    data <- right_pad_and_unlist(data, padding_string=NA_character_)
    matrix(data, nrow=length(lines), byrow=TRUE)
}

### Read IgBLAST broken annotation file.
### See IMPORTANT NOTES at the top of this file.
.read_broken_igdata_file <- function(filepath, col2class)
{
    stopifnot(isSingleNonWhiteString(filepath),
              is.character(col2class), !is.null(names(col2class)))
    m <- .read_jagged_table_as_character_matrix(filepath)
    if (ncol(m) != length(col2class))
        stop(wmsg("error loading ", filepath, ": unexpected ",
                  "number of fields"))
    matrix2df(m, col2class)
}

read_igdata <- function(filepath, col2class)
{
    df <- .read_broken_igdata_file(filepath, col2class)
    .replace_IgBLAST_special_values_with_NAs(df)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### write_igdata()
###

write_igdata <- function(df, file="")
{
    df <- .replace_NAs_with_IgBLAST_special_values(df)
    stopifnot(isSingleNonWhiteString(file))
    header <- paste0("#", paste(colnames(df), collapse=", "))
    cat(header, "\n", sep="", file=file, append=FALSE)
    ## Thanks to .replace_NAs_with_IgBLAST_special_values(), 'df' is
    ## guaranteed to be free of NAs, empty strings, and strings with
    ## whitespaces.
    ## This is important because:
    ##   (a) Having NAs in the resulting file would crash IgBLAST internal
    ##       parsing code!
    ##   (b) Strings are stored **unquoted** in IgBLAST annotation files, and
    ##       the field delimiters are blocks of whitespaces. This means that
    ##       there's no way to store an empty or white string in such file.
    ##       Also storing a string with inner whitespaces in such file would
    ##       be seen as two separate values by the parsing code.
    write.table(df, file, append=TRUE, quote=FALSE,
                sep="\t", row.names=FALSE, col.names=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### check_and_reorder_igdata_rows()
###

### Used in write_ndm_data_to_db() and write_auxdata_to_db() to put
### the rows of the data.frame in the same order as the alleles in the
### corresponding region db.
check_and_reorder_igdata_rows <- function(df, db_allele_names)
{
    .check_igdata_first_col(df)
    allele_names <- df[ , "allele_name"]
    stopifnot(is.character(db_allele_names),
              !anyDuplicated(allele_names), !anyDuplicated(db_allele_names),
              nrow(df) == length(db_allele_names),  # not really necessary
              setequal(allele_names, db_allele_names))
    m <- match(db_allele_names, allele_names)
    S4Vectors:::extract_data_frame_rows(df, m)
}

