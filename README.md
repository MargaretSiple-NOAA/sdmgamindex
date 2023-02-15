<!-- README.md is generated from README.Rmd. Please edit that file -->

## surveyIndex

R package for calculating survey indices by age from DATRAS exchange
data.

[![](https://img.shields.io/github/last-commit/EmilyMarkowitz-NOAA/surveyIndex.svg)](https://github.com/EmilyMarkowitz-NOAA/surveyIndex/commits/main)

> Code is still in development at
> <https://github.com/EmilyMarkowitz-NOAA/surveyIndex>

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

    ## R version 4.2.2 (2022-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] knitr_1.42        badger_0.2.3      surveyIndex_0.1.0 RODBC_1.3-20      roxygen2_7.2.3    devtools_2.4.5    here_1.0.1        gstat_2.1-0       purrr_1.0.1       ggplot2_3.4.0     dplyr_1.0.10     
    ## [12] visreg_2.7.0      gratia_0.8.1      sdmTMB_0.3.0      tweedie_2.3.5     MASS_7.3-58.2     marmap_1.0.9      RANN_2.6.1        maptools_1.1-6    sp_1.6-0          mapdata_2.3.1     maps_3.4.1       
    ## [23] mgcv_1.8-41       nlme_3.1-161      DATRAS_1.01       usethis_2.1.6    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] lwgeom_0.2-11       inlabru_2.7.0       plyr_1.8.8          igraph_1.3.5        splines_4.2.2       digest_0.6.31       yulab.utils_0.0.6   htmltools_0.5.4     fansi_1.0.4         magrittr_2.0.3     
    ##  [11] memoise_2.0.1       gert_1.9.2          remotes_2.4.2       credentials_1.3.2   xts_0.12.2          askpass_1.1         pkgdown_2.0.7       prettyunits_1.1.1   colorspace_2.1-0    blob_1.2.3         
    ##  [21] gitcreds_0.1.2      xfun_0.36           dlstats_0.1.6       callr_3.7.3         crayon_1.5.2        jsonlite_1.8.4      zoo_1.8-11          glue_1.6.2          stars_0.6-0         gtable_0.3.1       
    ##  [31] pkgbuild_1.4.0      shape_1.4.6         abind_1.4-5         scales_1.2.1        DBI_1.1.3           ggeffects_1.1.5     miniUI_0.1.1.1      Rcpp_1.0.10         xtable_1.8-4        units_0.8-1        
    ##  [41] foreign_0.8-84      bit_4.0.5           proxy_0.4-27        intervals_0.15.2    profvis_0.3.7       htmlwidgets_1.6.1   httr_1.4.4          RColorBrewer_1.1-3  FNN_1.1.3.1         ellipsis_0.3.2     
    ##  [51] urlchecker_1.0.1    pkgconfig_2.0.3     utf8_1.2.2          tidyselect_1.2.0    rlang_1.0.6         reshape2_1.4.4      later_1.3.0         munsell_0.5.0       tools_4.2.2         cachem_1.0.6       
    ##  [61] cli_3.6.0           generics_0.1.3      RSQLite_2.2.20      evaluate_0.20       stringr_1.5.0       fastmap_1.1.0       yaml_2.3.7          sys_3.4.1           processx_3.8.0      bit64_4.0.5        
    ##  [71] fs_1.6.0            ncdf4_1.21          gh_1.3.1            whisker_0.4.1       mime_0.12           mvnfast_0.2.7       xml2_1.3.3          brio_1.1.3          compiler_4.2.2      rstudioapi_0.14    
    ##  [81] curl_5.0.0          testthat_3.1.6      e1071_1.7-12        tibble_3.1.8        spacetime_1.2-8     stringi_1.7.12      gdistance_1.6       ps_1.7.2            desc_1.4.2          lattice_0.20-45    
    ##  [91] Matrix_1.5-3        classInt_0.4-8      vctrs_0.5.2         pillar_1.8.1        lifecycle_1.0.3     BiocManager_1.30.19 raster_3.6-14       httpuv_1.6.8        patchwork_1.1.2     R6_2.5.1           
    ## [101] promises_1.2.0.1    icesDatras_1.4.0    KernSmooth_2.23-20  sessioninfo_1.2.2   codetools_0.2-18    assertthat_0.2.1    pkgload_1.3.2       openssl_2.0.5       rprojroot_2.0.3     withr_2.5.0        
    ## [111] adehabitatMA_0.3.15 terra_1.7-3         grid_4.2.2          tidyr_1.3.0         class_7.3-21        rvcheck_0.2.1       rmarkdown_2.20      sf_1.0-9            shiny_1.7.4

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
