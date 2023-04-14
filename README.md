<!-- README.md is generated from README.Rmd. Please edit that file -->

# sdmgamindex <a href="https://emilymarkowitz-noaa.github.io/sdmgamindex/"><img src="man/figures/logo.png" align="right" height="133" /></a>

(formally {surveyIndex})

R package for calculating survey indices by age from DATRAS exchange
data.

[![](https://img.shields.io/github/last-commit/EmilyMarkowitz-NOAA/sdmgamindex.svg)](https://github.com/EmilyMarkowitz-NOAA/sdmgamindex/commits/main)

> Code is still in development at
> <https://github.com/EmilyMarkowitz-NOAA/sdmgamindex>

*Code was originally developed by:*

**Casper W. Berg** (@casperwberg)

National Institute of Aquatic Resources,

Technical University of Denmark

[**Berg et al. (2014): “Evaluation of alternative age-based methods for
estimating relative abundance from survey data in relation to assessment
models”, Fisheries Research 151(2014)
91-99.**](https://doi.org/10.1016/j.fishres.2013.10.005)

*And then modified and adapted for the AFSC by:*

**Emily Markowitz** (@EmilyMarkowitz-noaa; Emily.Markowitz AT noaa.gov)

**Margaret Siple** (@MargaretSiple-noaa; Margaret.Siple AT noaa.gov)

Alaska Fisheries Science Center

National Marine Fisheries Service

National Oceanic and Atmospheric Administration

Seattle, WA 98195

## Installation

To install you need the DATRAS package

``` r
library(remotes)
remotes::install_github("DTUAqua/DATRAS/DATRAS")
# remotes::install_github("casperwberg/sdmgamindex/sdmgamindex")
remotes::install_github("emilymarkowitz-noaa/sdmgamindex")
```

## Usage

See examples in the [pkgdown
site](https://EmilyMarkowitz-NOAA.github.io/sdmgamindex/) and in the
[?get_surveyidx()
documentation](https://emilymarkowitz-noaa.github.io/sdmgamindex/reference/get_surveyidx.html).

## Metadata

This package was last produced using:

    ## R version 4.2.3 (2023-03-15 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] knitr_1.42        badger_0.2.3      sdmgamindex_0.1.0 tweedie_2.3.5    
    ##  [5] MASS_7.3-58.3     marmap_1.0.10     RANN_2.6.1        maptools_1.1-6   
    ##  [9] sp_1.6-0          mapdata_2.3.1     maps_3.4.1        mgcv_1.8-42      
    ## [13] nlme_3.1-162      DATRAS_1.01       RODBC_1.3-20      roxygen2_7.2.3   
    ## [17] devtools_2.4.5    usethis_2.1.6     here_1.0.1       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] colorspace_2.1-0    gitcreds_0.1.2      ellipsis_0.3.2     
    ##   [4] rprojroot_2.0.3     fs_1.6.1            rstudioapi_0.14    
    ##   [7] remotes_2.4.2       gh_1.4.0            bit64_4.0.5        
    ##  [10] fansi_1.0.4         xml2_1.3.3          codetools_0.2-19   
    ##  [13] splines_4.2.3       ncdf4_1.21          cachem_1.0.7       
    ##  [16] pkgload_1.3.2       jsonlite_1.8.4      icesDatras_1.4.0   
    ##  [19] shiny_1.7.4         BiocManager_1.30.20 compiler_4.2.3     
    ##  [22] rvcheck_0.2.1       Matrix_1.5-3        fastmap_1.1.1      
    ##  [25] cli_3.6.1           later_1.3.0         htmltools_0.5.5    
    ##  [28] prettyunits_1.1.1   tools_4.2.3         igraph_1.4.1       
    ##  [31] gtable_0.3.3        glue_1.6.2          reshape2_1.4.4     
    ##  [34] dplyr_1.1.1         rappdirs_0.3.3      Rcpp_1.0.10        
    ##  [37] pkgdown_2.0.7       raster_3.6-20       vctrs_0.6.1        
    ##  [40] xfun_0.38           stringr_1.5.0       ps_1.7.3           
    ##  [43] brio_1.1.3          testthat_3.1.7      mime_0.12          
    ##  [46] miniUI_0.1.1.1      lifecycle_1.0.3     sys_3.4.1          
    ##  [49] dlstats_0.1.6       terra_1.7-18        scales_1.2.1       
    ##  [52] credentials_1.3.2   promises_1.2.0.1    httr2_0.2.2        
    ##  [55] RColorBrewer_1.1-3  gert_1.9.2          yaml_2.3.7         
    ##  [58] curl_5.0.0          memoise_2.0.1       ggplot2_3.4.1      
    ##  [61] yulab.utils_0.0.6   stringi_1.7.12      RSQLite_2.3.0      
    ##  [64] desc_1.4.2          pkgbuild_1.4.0      shape_1.4.6        
    ##  [67] rlang_1.1.0         pkgconfig_2.0.3     evaluate_0.20      
    ##  [70] lattice_0.21-8      purrr_1.0.1         htmlwidgets_1.6.2  
    ##  [73] bit_4.0.5           processx_3.8.0      tidyselect_1.2.0   
    ##  [76] plyr_1.8.8          magrittr_2.0.3      R6_2.5.1           
    ##  [79] generics_0.1.3      profvis_0.3.7       DBI_1.1.3          
    ##  [82] pillar_1.9.0        foreign_0.8-84      withr_2.5.0        
    ##  [85] tibble_3.2.1        crayon_1.5.2        utf8_1.2.3         
    ##  [88] rmarkdown_2.20      urlchecker_1.0.1    grid_4.2.3         
    ##  [91] blob_1.2.4          callr_3.7.3         digest_0.6.31      
    ##  [94] xtable_1.8-4        gdistance_1.6       adehabitatMA_0.3.16
    ##  [97] httpuv_1.6.9        openssl_2.0.6       munsell_0.5.0      
    ## [100] sessioninfo_1.2.2   askpass_1.1

## NOAA README

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.

## NOAA License

Software code created by U.S. Government employees is not subject to
copyright in the United States (17 U.S.C. §105). The United
States/Department of Commerce reserve all rights to seek and obtain
copyright protection in countries other than the United States for
Software authored in its entirety by the Department of Commerce. To this
end, the Department of Commerce hereby grants to Recipient a
royalty-free, nonexclusive license to use, copy, and create derivative
works of the Software outside of the United States.

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" alt="NOAA Fisheries" height="75"/>

[U.S. Department of Commerce](https://www.commerce.gov/) \| [National
Oceanographic and Atmospheric Administration](https://www.noaa.gov) \|
[NOAA Fisheries](https://www.fisheries.noaa.gov/)
