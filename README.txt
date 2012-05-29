How to build the R package
==========================

1. Run roxygen2 to generate/update the R documentation:

    > setwd("/path/to/Rngqs/parent/dir")
    > roxygenize("Rngq")

2. Build the package

    cd /path/to/Rngqs/parent/dir
    R CMD build Rngq

3. Check the package (you must supply the correct version number):

    R CMD check Rngq_x.y.tar.gz

4. Install the package locally

    R CMD INSTALL Rngq_x.y.tar.gz

-M
