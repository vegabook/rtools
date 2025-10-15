# T Browne financial analysis tools for R using Bloomberg Terminal #

## Installation
* Clone the repo
* In R:  
  `source("path/to/repo/tools.r")`  
  `source("path/to/repo/common.r")`  
  For each library that does not exist, you will have to `install.packages("offending_library")`  
* Have a look at the code but the main usefull tools are `blpConnect`, `bbdh` for daily historical data, `bbdt` for ticks, `bds` reference data`, and `regress` for regression charting.
* There is a bunch of other linear algebra stuff in there specifically PCA-related like `makePCs`, and also convenience functions like `logret`. 



