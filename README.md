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
    ## [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] knitr_1.44         badger_0.2.3       pkgdown_2.0.7.9000 usethis_2.2.2     
    ##  [5] here_1.0.1         sdmgamindex_0.1.1  tweedie_2.3.5      MASS_7.3-60       
    ##  [9] sf_1.0-14          sp_2.0-0           RANN_2.6.1         marmap_1.0.10     
    ## [13] mapdata_2.3.1      maps_3.4.1         tidyr_1.3.0        mgcv_1.8-42       
    ## [17] nlme_3.1-162      
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3      sys_3.4.2               rstudioapi_0.15.0      
    ##   [4] dlstats_0.1.7           jsonlite_1.8.7          shape_1.4.6            
    ##   [7] magrittr_2.0.3          rmarkdown_2.25          fs_1.6.3               
    ##  [10] ragg_1.2.5              vctrs_0.6.3             memoise_2.0.1          
    ##  [13] askpass_1.2.0           terra_1.7-46            gh_1.4.0               
    ##  [16] janitor_2.2.0           htmltools_0.5.6         curl_5.0.2             
    ##  [19] raster_3.6-23           KernSmooth_2.23-22      httr2_0.2.3            
    ##  [22] plyr_1.8.8              stars_0.6-4             lubridate_1.9.2        
    ##  [25] cachem_1.0.8            uuid_1.1-1              igraph_1.5.1           
    ##  [28] mime_0.12               lifecycle_1.0.3         pkgconfig_2.0.3        
    ##  [31] Matrix_1.5-4.1          R6_2.5.1                fastmap_1.1.1          
    ##  [34] shiny_1.7.5             snakecase_0.11.1        digest_0.6.33          
    ##  [37] colorspace_2.1-0        patchwork_1.1.3         rprojroot_2.0.3        
    ##  [40] textshaping_0.3.6       RSQLite_2.3.1           fansi_1.0.4            
    ##  [43] timechange_0.2.0        gdistance_1.6.4         abind_1.4-5            
    ##  [46] compiler_4.3.1          proxy_0.4-27            bit64_4.0.5            
    ##  [49] fontquiver_0.2.1        DBI_1.1.3               gratia_0.8.1           
    ##  [52] openssl_2.1.0           rappdirs_0.3.3          classInt_0.4-10        
    ##  [55] gfonts_0.2.0            tools_4.3.1             units_0.8-3            
    ##  [58] zip_2.3.0               httpuv_1.6.11           glue_1.6.2             
    ##  [61] promises_1.2.1          grid_4.3.1              reshape2_1.4.4         
    ##  [64] generics_0.1.3          gtable_0.3.4            class_7.3-22           
    ##  [67] data.table_1.14.8       xml2_1.3.5              utf8_1.2.3             
    ##  [70] pillar_1.9.0            stringr_1.5.0           yulab.utils_0.1.0      
    ##  [73] later_1.3.1             splines_4.3.1           dplyr_1.1.3            
    ##  [76] lattice_0.21-8          bit_4.0.5               tidyselect_1.2.0       
    ##  [79] rvcheck_0.2.1           fontLiberation_0.1.0    gitcreds_0.1.2         
    ##  [82] fontBitstreamVera_0.1.1 crul_1.4.0              xfun_0.40              
    ##  [85] adehabitatMA_0.3.16     credentials_2.0.1       stringi_1.7.12         
    ##  [88] yaml_2.3.7              evaluate_0.22           codetools_0.2-19       
    ##  [91] httpcode_0.3.0          officer_0.6.2           gdtools_0.3.3          
    ##  [94] tibble_3.2.1            BiocManager_1.30.22     cli_3.6.1              
    ##  [97] xtable_1.8-4            systemfonts_1.0.4       munsell_0.5.0          
    ## [100] Rcpp_1.0.11             gert_2.0.0              ellipsis_0.3.2         
    ## [103] ggplot2_3.4.3           blob_1.2.4              mvnfast_0.2.8          
    ## [106] scales_1.2.1            e1071_1.7-13            ncdf4_1.21             
    ## [109] purrr_1.0.2             crayon_1.5.2            flextable_0.9.3        
    ## [112] rlang_1.1.1

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
