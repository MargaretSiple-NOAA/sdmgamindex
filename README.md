<!-- README.md is generated from README.Rmd. Please edit that file -->

## surveyIndex

R package for calculating survey indices by age from DATRAS exchange
data.

[![](https://img.shields.io/github/last-commit/EmilyMarkowitz-NOAA/surveyIndex.svg)](https://github.com/EmilyMarkowitz-NOAA/surveyIndex/commits/main)

> Code is still in development

Code was originally developed by: **Casper W. Berg** (@casperwberg)

> Berg et al. (2014): “Evaluation of alternative age-based methods for
> estimating relative abundance from survey data in relation to
> assessment models”, Fisheries Research 151(2014) 91-99.

And then modified and adapted for the AFSC by:

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
# remotes::install_github("casperwberg/surveyIndex/surveyIndex")
remotes::install_github("emilymarkowitz-noaa/surveyIndex")
```

## Usage

See examples in the [pkgdown
site](https://EmilyMarkowitz-NOAA.github.io/surveyIndex/) and in the
[?get_surveyidx()
documentation](https://emilymarkowitz-noaa.github.io/surveyIndex/reference/get_surveyidx.html).

## Metadata

This package was last produced using:

    ## R version 4.2.1 (2022-06-23 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] knitr_1.40        badger_0.2.1      surveyIndex_0.1.0 marmap_1.0.6      RANN_2.6.1        maptools_1.1-4   
    ##  [7] sp_1.5-0          mapdata_2.3.0     maps_3.4.0        mgcv_1.8-40       nlme_3.1-157      DATRAS_1.01      
    ## [13] RODBC_1.3-19      roxygen2_7.2.1    devtools_2.4.4    usethis_2.1.6     here_1.0.1       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] colorspace_2.0-3    gitcreds_0.1.2      ellipsis_0.3.2      rprojroot_2.0.3     fs_1.5.2           
    ##   [6] rstudioapi_0.14     remotes_2.4.2       gh_1.3.1            bit64_4.0.5         fansi_1.0.3        
    ##  [11] xml2_1.3.3          codetools_0.2-18    splines_4.2.1       ncdf4_1.19          cachem_1.0.6       
    ##  [16] pkgload_1.3.0       jsonlite_1.8.0      icesDatras_1.4.0    shiny_1.7.2         httr_1.4.4         
    ##  [21] BiocManager_1.30.18 compiler_4.2.1      rvcheck_0.2.1       assertthat_0.2.1    Matrix_1.5-0       
    ##  [26] fastmap_1.1.0       cli_3.4.0           later_1.3.0         htmltools_0.5.3     prettyunits_1.1.1  
    ##  [31] tools_4.2.1         gtable_0.3.1        glue_1.6.2          reshape2_1.4.4      dplyr_1.0.10       
    ##  [36] Rcpp_1.0.9          raster_3.5-29       vctrs_0.4.1         xfun_0.32           stringr_1.4.1      
    ##  [41] ps_1.7.1            brio_1.1.3          testthat_3.1.4      mime_0.12           miniUI_0.1.1.1     
    ##  [46] lifecycle_1.0.2     sys_3.4             dlstats_0.1.5       terra_1.6-17        MASS_7.3-57        
    ##  [51] scales_1.2.1        promises_1.2.0.1    credentials_1.3.2   RColorBrewer_1.1-3  gert_1.8.0         
    ##  [56] yaml_2.3.5          curl_4.3.2          memoise_2.0.1       ggplot2_3.3.6       yulab.utils_0.0.5  
    ##  [61] stringi_1.7.8       RSQLite_2.2.17      desc_1.4.2          pkgbuild_1.3.1      shape_1.4.6        
    ##  [66] rlang_1.0.5         pkgconfig_2.0.3     evaluate_0.16       lattice_0.20-45     purrr_0.3.4        
    ##  [71] htmlwidgets_1.5.4   bit_4.0.4           processx_3.7.0      tidyselect_1.1.2    plyr_1.8.7         
    ##  [76] magrittr_2.0.3      R6_2.5.1            generics_0.1.3      profvis_0.3.7       DBI_1.1.3          
    ##  [81] pillar_1.8.1        foreign_0.8-82      withr_2.5.0         tibble_3.1.8        crayon_1.5.1       
    ##  [86] utf8_1.2.2          rmarkdown_2.16      urlchecker_1.0.1    grid_4.2.1          blob_1.2.3         
    ##  [91] callr_3.7.2         digest_0.6.29       xtable_1.8-4        adehabitatMA_0.3.14 httpuv_1.6.6       
    ##  [96] openssl_2.0.2       munsell_0.5.0       tweedie_2.3.5       sessioninfo_1.2.2   askpass_1.1

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
