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

**Casper W. Berg**
(<https://orbit.dtu.dk/en/persons/casper-willestofte-berg>;
[@casperwberg](https://github.com/casperwberg))

National Institute of Aquatic Resources,

Technical University of Denmark

*Based on the work published in:* Berg et al. ([2014](#ref-Berg2014))

Repository:
`remotes::install_github("casperwberg/surveyIndex/surveyIndex")`

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
    ## [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] knitr_1.45     badger_0.2.3   pkgdown_2.0.7  tweedie_2.3.5  MASS_7.3-60    sf_1.0-15      sp_2.1-3       RANN_2.6.1    
    ##  [9] marmap_1.0.10  mapdata_2.3.1  maps_3.4.2     tidyr_1.3.1    mgcv_1.9-1     nlme_3.1-164   RODBC_1.3-23   roxygen2_7.3.1
    ## [17] devtools_2.4.5 usethis_2.2.3  here_1.0.1    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] DBI_1.2.2           httr2_1.0.0         remotes_2.4.2.1     rlang_1.1.3         magrittr_2.0.3      e1071_1.7-14       
    ##  [7] compiler_4.3.1      RSQLite_2.3.5       vctrs_0.6.5         reshape2_1.4.4      stringr_1.5.1       profvis_0.3.8      
    ## [13] crayon_1.5.2        pkgconfig_2.0.3     shape_1.4.6.1       fastmap_1.1.1       ellipsis_0.3.2      utf8_1.2.4         
    ## [19] promises_1.2.1      rmarkdown_2.25      sessioninfo_1.2.2   purrr_1.0.2         bit_4.0.5           gert_2.0.1         
    ## [25] xfun_0.42           rvcheck_0.2.1       cachem_1.0.8        jsonlite_1.8.8      blob_1.2.4          later_1.3.2        
    ## [31] terra_1.7-71        R6_2.5.1            stringi_1.8.3       RColorBrewer_1.1-3  pkgload_1.3.4       gdistance_1.6.4    
    ## [37] Rcpp_1.0.12         gitcreds_0.1.2      dlstats_0.1.7       httpuv_1.6.14       Matrix_1.6-5        splines_4.3.1      
    ## [43] igraph_2.0.2        tidyselect_1.2.0    rstudioapi_0.15.0   yaml_2.3.8          codetools_0.2-19    miniUI_0.1.1.1     
    ## [49] curl_5.2.0          pkgbuild_1.4.3      lattice_0.22-5      tibble_3.2.1        plyr_1.8.9          shiny_1.8.0        
    ## [55] askpass_1.2.0       evaluate_0.23       adehabitatMA_0.3.16 units_0.8-5         proxy_0.4-27        urlchecker_1.0.1   
    ## [61] xml2_1.3.6          pillar_1.9.0        BiocManager_1.30.22 KernSmooth_2.23-22  ncdf4_1.22          generics_0.1.3     
    ## [67] credentials_2.0.1   rprojroot_2.0.4     ggplot2_3.5.0       munsell_0.5.0       scales_1.3.0        xtable_1.8-4       
    ## [73] class_7.3-22        glue_1.7.0          tools_4.3.1         sys_3.4.2           fs_1.6.3            grid_4.3.1         
    ## [79] gh_1.4.0            colorspace_2.1-0    raster_3.6-26       cli_3.6.2           rappdirs_0.3.3      fansi_1.0.6        
    ## [85] dplyr_1.1.4         gtable_0.3.4        yulab.utils_0.1.4   digest_0.6.34       classInt_0.4-10     htmlwidgets_1.6.4  
    ## [91] memoise_2.0.1       htmltools_0.5.7     lifecycle_1.0.4     mime_0.12           openssl_2.1.1       bit64_4.0.5

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
