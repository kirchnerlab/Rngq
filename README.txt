How to build/install the R package
==========================

1. Run roxygen2 to generate/update the R documentation:

    > setwd("/path/to/Rngqs/parent/dir")
    > roxygenize("Rngq")

2. Build the package

    cd /path/to/Rngqs/parent/dir
    R CMD build Rngq

3. Check the package (you must supply the correct version number):

    R CMD check Rngq_x.y.tar.gz

   This may throw a few errors related to XLConnect which is not yet available
   for R 2.15 which I use on the build system. Simply ignore these errors.
   However, you should make sure that XLConnect is installed on your system
   if you would like to store peptide/protein results to MS Excel XML files
   (.xlsx).

4. Install the package locally

    R CMD INSTALL Rngq_x.y.tar.gz

-M
