# igblastr [<img src="https://bioconductor.org/images/logo/jpg/bioconductor_logo_rgb.jpg" width="200" align="right"/>](https://bioconductor.org/)

[![in_Bioc](https://bioconductor.org/shields/years-in-bioc/igblastr.svg)](https://bioconductor.org/packages/igblastr)
[![rank](https://bioconductor.org/shields/downloads/release/igblastr.svg)](https://bioconductor.org/packages/stats/bioc/igblastr)

**igblastr** is an R/Bioconductor package that provides functions to conveniently install and use a local IgBLAST installation from within R.

To cite the **igblastr** package in publications, please use:

Pagès H, Hyrien O, MacPhee K, Duff M, Taylor J (2025).
_igblastr: User-friendly R Wrapper to IgBLAST_.
doi:10.18129/B9.bioc.igblastr
<https://doi.org/10.18129/B9.bioc.igblastr>.
<https://bioconductor.org/packages/igblastr>.

A BibTeX entry for LaTeX users is
```
@Manual{,
  title = {igblastr: User-friendly R Wrapper to IgBLAST},
  author = {Hervé Pagès and Ollivier Hyrien and Kellie MacPhee and Michael Duff and Jason Taylor},
  year = {2025},
  doi = {10.18129/B9.bioc.igblastr},
  url = {https://bioconductor.org/packages/igblastr},
}
```
Please also cite the original IgBLAST tool if you use it in your analysis.

See <https://bioconductor.org/packages/igblastr> for more information,
including for how to install the release version of the package (please refrain
from installing directly from GitHub).

|BioC version|Status|
|------------|------|
| release    | [![build](https://bioconductor.org/shields/build/release/bioc/igblastr.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/igblastr/) [![dependencies](https://bioconductor.org/shields/dependencies/release/igblastr.svg)](https://bioconductor.org/packages/release/igblastr#since) |
| devel      | [![build](https://bioconductor.org/shields/build/devel/bioc/igblastr.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/igblastr/) [![dependencies](https://bioconductor.org/shields/dependencies/devel/igblastr.svg)](https://bioconductor.org/packages/devel/igblastr#since) |

IgBLAST is described at <https://pubmed.ncbi.nlm.nih.gov/23671333/>

IgBLAST web interface: <https://www.ncbi.nlm.nih.gov/igblast/>

Please use <https://github.com/HyrienLab/igblastr/issues> to report bugs,
provide feedback, request features, etc...


## Quick start


### 1. Install and load igblastr

#### Install igblastr

First install and load **BiocManager**:
```r
if (!require("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

library(BiocManager)
```

Loading **BiocManager** should display a message that indicates the
version of Bioconductor that you are using.

If you are using Bioconductor >= 3.22, simply install **igblastr** with:
```r
BiocManager::install("igblastr")
```
This will install the version of **igblastr** that matches your version
of Bioconductor.

If you are using Bioconductor 3.21, then you need to install directly
from GitHub (usually discouraged). That's because **igblastr** was only
added to Bioconductor >= 3.22 so cannot be installed from the Bioconductor
3.21 repositories:
```r
if (!require("remotes", quietly=TRUE))
    BiocManager::install("remotes")

BiocManager::install("HyrienLab/igblastr@RELEASE_3_22")
```

#### Load igblastr

```r
library(igblastr)
```


### 2. Set up IgBLAST

#### Install a pre-compiled IgBLAST

```r
install_igblast()
```

This will download the arch-specific pre-compiled IgBLAST from
<https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/>
and install it in a persistent location. See `?install_igblast`.

#### or use an existing installation of IgBLAST

```r
Sys.setenv(IGBLAST_ROOT="path/to/igblast/root")
```

Note: Environment variable `IGBLAST_ROOT` should preferrably be set
in a more persistent manner outside R e.g. by defining it in the
user's `.profile` (on Unix/Mac). See `?IGBLAST_ROOT`.

#### Check IgBLAST

```r
igblast_info()
```

This will make sure that the installed IgBLAST is functional. It will
also display basic information about it (e.g. `igblastn` version).


### 3. Install germline dbs and select db to use

#### Install germline dbs

At least one germline db must be installed.

Several specialized functions will be provided for that e.g.:

```r
install_IMGT_germline_db()   # see '?install_IMGT_germline_db'
install_OGRDB_germline_db()  # not ready yet!
# etc...
```

#### Select the germline db to use with igblastn()

```r
list_germline_dbs()
use_germline_db(<db-name>)
```

See `?use_germline_db`.


### 4. Select constant region db (optional)

The C-regions from IMGT are pre-installed for a few organisms:
```r
list_c_region_dbs()
```

Select the C-region db to use with `igblastn()`:
```r
use_c_region_db(<db-name>)
```

See `?use_c_region_db`.


### 5. Use igblastn()

See `?igblastn`.

