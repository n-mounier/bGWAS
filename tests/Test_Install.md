# Test installation of the package on different OS / R versions:


## Test on my laptop
macOS, version 10.13.6    
- R version 3.4.3 (with `TwoSampleMR` already installed)   
ok

## Test using travis CI

Note, August 2019: "release"=3.6.1 - "oldrel"=3.5.3    
<!---
rversions::r_release()
rversions::r_oldrel()
--->



## Test on rstudio cloud   
Ubuntu, version 6.04.6 LTS (Xenial Xerus)   
- R version 3.6.0   
ok   
- R version 3.5.3   
ok   
- R version 3.4.4   
Failed to install 'TwoSampleMR' from GitHub:   
 (converted from warning) package ‘metafor’ is not available (for R version 3.4.4)   
- R version 3.3.3   
Failed to install 'TwoSampleMR' from GitHub:      
 (converted from warning) package ‘metafor’ is not available (for R version 3.3.3)   
- R version 3.2.5   
Failed to install 'TwoSampleMR' from GitHub:      
 (converted from warning) package ‘metafor’ is not available (for R version 3.3.3)   


## Test on HPC1
CentOS Linux, version 7
- R version 3.4.3 (with `TwoSampleMR` already installed)   
ok


<!---

install.packages("remotes")
remotes::install_github("n-mounier/bGWAS")
library(bGWAS)



cat /etc/os-release
--->


# Dependencies and versions

- DEPENDS    
`dplyr` (>= 3.3.0) -> [R 3.2.0](https://cran.r-project.org/web/packages/dplyr/dplyr.pdf)    
`magrittr` (>= 1.5) -> [no version](https://cran.r-project.org/web/packages/magrittr/magrittr.pdf)    

- SUGGESTS    
`testthat` (...) -> [R 3.1](https://cran.r-project.org/web/packages/testthat/testthat.pdf)    
`knitr` (...) -> [R 3.2.3](https://cran.r-project.org/web/packages/knitr/knitr.pdf)     
`remotes` (...) -> [R 3.0](https://cran.r-project.org/web/packages/remotes/remotes.pdf)    
`rmarkdown` (...) -> [R 3.0](https://cran.r-project.org/web/packages/rmarkdown/rmarkdown.pdf)    

- IMPORTS
`calibrate` (>= 1.7.2) -> [R 1.8.0](https://cran.r-project.org/web/packages/calibrate/calibrate.pdf)    
`data.table` (>= 1.12.0) -> [R 3.1.0](https://cran.r-project.org/web/packages/data.table/data.table.pdf)    
`ggplot2` (>= 2.2.1) -> [R 3.2.1](https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf)    
`gplots` (>= 3.0.1) -> [R 3.0](https://cran.r-project.org/web/packages/gplots/gplots.pdf)     
`qqman` (>= 0.1.4) -> [R 3.0.0](https://cran.r-project.org/web/packages/qqman/qqman.pdf)    
`readr` (>= 1.3.1) -> [R 3.1](https://cran.r-project.org/web/packages/readr/readr.pdf)    
`Rcpp` (>= 0.12.15) -> [R 3.0.0](https://cran.r-project.org/web/packages/Rcpp/Rcpp.pdf)      
`rlang` (>= 0.4.0)  -> [R 3.2](https://cran.r-project.org/web/packages/rlang/rlang.pdf)   
`R.utils` (>= 2.9.0) -> [R 2.14.0](https://cran.r-project.org/web/packages/R.utils/R.utils.pdf)     
`stringr` (>= 1.4.0) -> [R 3.1](https://cran.r-project.org/web/packages/stringr/stringr.pdf)    
`tibble` (>= 2.1.1) -> [R 3.1.0](https://cran.r-project.org/web/packages/tibble/tibble.pdf)    
`tidyr` (>= 0.8.3) -> [R 3.1](https://cran.r-project.org/web/packages/tidyr/tidyr.pdf)     
`TwoSampleMR` (>= 0.3.0) -> [R 3.1](https://github.com/MRCIEU/TwoSampleMR/blob/master/DESCRIPTION)  


According to this, the lowest version needed would be R 3.2.3 ... but actually, `TwoSampleMR` requires `meta` that requires `metafor` ([R 3.5.0](https://cran.r-project.org/web/packages/metafor/metafor.pdf))...

