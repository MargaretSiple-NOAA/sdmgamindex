<!-- README.md is generated from README.Rmd. Please edit that file -->

# [sdmgamindex](link_repo) <img src="https://avatars.githubusercontent.com/u/91760178?s=96&amp;v=4" alt="Logo." align="right" width="139" height="139"/>

This repository was previously forked from `casperwberg/surveyIndex` and
previously named `emilymarkowitz-noaa/surveyIndex`.

R package for calculating survey indices by age from DATRAS exchange
data.

[![](https://img.shields.io/github/last-commit/EmilyMarkowitz-NOAA/sdmgamindex.svg)](https://github.com/EmilyMarkowitz-NOAA/sdmgamindex/commits/main)

> This code is always in development. Find code used for various reports
> in the code
> [releases](https://github.com/EmilyMarkowitz-NOAA/sdmgamindex//releases).

### Code has been modified and adapted by:

**Emily Markowitz** (Emily.Markowitz AT noaa.gov;
[@EmilyMarkowitz-NOAA](https://github.com/EmilyMarkowitz-NOAA))

**Margaret Siple** (Margaret.Siple AT noaa.gov;
[@MargaretSiple-noaa](https://github.com/MargaretSiple-noaa))

Alaska Fisheries Science Center

National Marine Fisheries Service

National Oceanic and Atmospheric Administration

Seattle, WA 98195

### Code was originally developed by:

**Casper W. Berg** (([**casperwberg?**](#ref-casperwberg))) Berg et al.
([2014](#ref-Berg2014))

National Institute of Aquatic Resources,

Technical University of Denmark

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

    ## R version 4.3.1 (2023-06-16 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] knitr_1.44   badger_0.2.3
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.4        httr2_0.2.3         xfun_0.40           ggplot2_3.4.3       htmlwidgets_1.6.2  
    ##  [6] devtools_2.4.5      remotes_2.4.2.1     gh_1.4.0            processx_3.8.2      callr_3.7.3        
    ## [11] yulab.utils_0.1.0   vctrs_0.6.3         tools_4.3.1         ps_1.7.5            generics_0.1.3     
    ## [16] curl_5.0.2          tibble_3.2.1        fansi_1.0.4         pkgconfig_2.0.3     RColorBrewer_1.1-3 
    ## [21] lifecycle_1.0.3     compiler_4.3.1      stringr_1.5.0       credentials_2.0.1   dlstats_0.1.7      
    ## [26] munsell_0.5.0       httpuv_1.6.11       sys_3.4.2           htmltools_0.5.6     usethis_2.2.2      
    ## [31] yaml_2.3.7          later_1.3.1         pkgdown_2.0.7.9000  pillar_1.9.0        crayon_1.5.2       
    ## [36] urlchecker_1.0.1    ellipsis_0.3.2      openssl_2.1.0       cachem_1.0.8        sessioninfo_1.2.2  
    ## [41] mime_0.12           tidyselect_1.2.0    digest_0.6.33       stringi_1.7.12      dplyr_1.1.3        
    ## [46] purrr_1.0.2         rprojroot_2.0.3     fastmap_1.1.1       grid_4.3.1          here_1.0.1         
    ## [51] colorspace_2.1-0    cli_3.6.1           magrittr_2.0.3      pkgbuild_1.4.2      utf8_1.2.3         
    ## [56] rappdirs_0.3.3      prettyunits_1.2.0   scales_1.2.1        promises_1.2.1      rmarkdown_2.25     
    ## [61] gitcreds_0.1.2      rvcheck_0.2.1       askpass_1.2.0       memoise_2.0.1       shiny_1.7.5        
    ## [66] evaluate_0.22       miniUI_0.1.1.1      profvis_0.3.8       rlang_1.1.1         gert_2.0.0         
    ## [71] Rcpp_1.0.11         xtable_1.8-4        glue_1.6.2          BiocManager_1.30.22 pkgload_1.3.3      
    ## [76] rstudioapi_0.15.0   jsonlite_1.8.7      R6_2.5.1            fs_1.6.3

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

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-Berg2014" class="csl-entry">

Berg, C. W., Nielsen, A., and Kristensen, K. (2014). Evaluation of
alternative age-based methods for estimating relative abundance from
survey data in relation to assessment models. *Fish. Res.*, *151*,
91–99. <https://doi.org/10.1016/j.fishres.2013.10.005>

</div>

</div>
