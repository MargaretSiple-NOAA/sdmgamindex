
# Connect to oracle ------------------------------------------------------------

# This has a specific username and password because I DONT want people to have access to this!
# source("https://raw.githubusercontent.com/afsc-gap-products/metadata/main/code/functions_oracle.R")

library(magrittr)
library(readr)
library(dplyr)

locations <- c("Z:/Projects/ConnectToOracle.R")
for (i in 1:length(locations)){
  if (file.exists(locations[i])){
    source(locations[i])
  }
}

# Load column metadata table ---------------------------------------------------

metadata_table_comment <- dplyr::bind_rows(
  # tables
  RODBC::sqlQuery(
    channel = channel_products,
    query = "SELECT table_name, comments
FROM all_tab_comments
WHERE owner = 'GAP_PRODUCTS'
ORDER BY table_name") %>% 
    data.frame(), 
  # materialized view
  RODBC::sqlQuery(
    channel = channel_products,
    query = "SELECT *FROM user_mview_comments") %>% 
    data.frame() %>% 
    dplyr::rename(TABLE_NAME = MVIEW_NAME) ) 

metadata_colname <- RODBC::sqlQuery(
  channel = channel_products, 
  query = "SELECT * FROM GAP_PRODUCTS.METADATA_COLUMN") %>% 
  janitor::clean_names()

# Load new FOSS data -----------------------------------------------------------

noaa_afsc_public_foss <- RODBC::sqlQuery(
    channel = channel_products, 
    query = 
      "SELECT 
hh.YEAR,
hh.SRVY,                 
hh.SURVEY,
hh.SURVEY_DEFINITION_ID,
-- hh.SURVEY_NAME,
hh.CRUISE,
hh.CRUISEJOIN,           
hh.HAUL,
hh.HAULJOIN,
hh.STRATUM,
hh.STATION,
hh.VESSEL_ID,
hh.VESSEL_NAME,          
hh.DATE_TIME,
hh.LATITUDE_DD_START, 
hh.LONGITUDE_DD_START, 
hh.LATITUDE_DD_END,
hh.LONGITUDE_DD_END, 
hh.BOTTOM_TEMPERATURE_C,
hh.SURFACE_TEMPERATURE_C,
hh.DEPTH_M,
cc.SPECIES_CODE,
ss.ITIS,
ss.WORMS,
ss.COMMON_NAME,     
ss.SCIENTIFIC_NAME,
ss.ID_RANK,
CASE WHEN cc.CPUE_KGKM2 IS NULL THEN 0 ELSE cc.CPUE_KGKM2 END AS CPUE_KGKM2,
CASE WHEN cc.CPUE_NOKM2 IS NULL THEN 0 ELSE cc.CPUE_NOKM2 END AS CPUE_NOKM2,
CASE WHEN cc.COUNT IS NULL THEN 0 ELSE cc.COUNT END AS COUNT,
CASE WHEN cc.WEIGHT_KG IS NULL THEN 0 ELSE cc.WEIGHT_KG END AS WEIGHT_KG,
CASE WHEN cc.TAXON_CONFIDENCE IS NULL THEN NULL ELSE cc.TAXON_CONFIDENCE END AS TAXON_CONFIDENCE,
hh.AREA_SWEPT_KM2,       
hh.DISTANCE_FISHED_KM,
hh.DURATION_HR,          
hh.NET_WIDTH_M,
hh.NET_HEIGHT_M,
hh.PERFORMANCE 
FROM GAP_PRODUCTS.FOSS_SURVEY_SPECIES sv
FULL OUTER JOIN GAP_PRODUCTS.FOSS_SPECIES ss
ON sv.SPECIES_CODE = ss.SPECIES_CODE
FULL OUTER JOIN GAP_PRODUCTS.FOSS_HAUL hh
ON sv.SURVEY_DEFINITION_ID = hh.SURVEY_DEFINITION_ID
FULL OUTER JOIN GAP_PRODUCTS.FOSS_CATCH cc
ON sv.SPECIES_CODE = cc.SPECIES_CODE
AND hh.HAULJOIN = cc.HAULJOIN
WHERE sv.SURVEY_DEFINITION_ID = 98 
AND sv.SPECIES_CODE IN (21740, 10210, 69322) 
AND hh.YEAR >= 1982 --2015
--GROUP BY ss.COMMON_NAME, hh.HAULJOIN
") %>% 
  janitor::clean_names()

# Save table to local directory
save(noaa_afsc_public_foss, file = here::here("data", "noaa_afsc_public_foss.rda"))
  
column <- metadata_colname %>%
  dplyr::filter(metadata_colname %in% toupper(names(noaa_afsc_public_foss))) %>%
  dplyr::mutate(metadata_colname = tolower(metadata_colname)) %>%
  dplyr::distinct()

str0 <- paste0("#' @title Public data from FOSS for EBS walleye pollock, yellowfin sole, and red king crab from 1982 to present
#' @description ",metadata_table_comment$COMMENT[metadata_table_comment$TABLE_NAME == "FOSS_CATCH"]," 
#' @usage data('noaa_afsc_public_foss')
#' @author Emily Markowitz (Emily.Markowitz AT noaa.gov)
#' @format A data frame with ",nrow(noaa_afsc_public_foss)," observations on the following ",
ncol(noaa_afsc_public_foss)," variables.
#' \\describe{
",
               paste0(paste0("#'   \\item{\\code{",column$metadata_colname,"}}{", column$metadata_colname_long, ". ", column$metadata_colname_desc,"}"), collapse = "\n"),
               "#'   }
#' @source https://github.com/afsc-gap-products/gap_products
#' @keywords species code data
#' @examples
#' data(noaa_afsc_public_foss)
#' @details The Resource Assessment and Conservation Engineering (RACE) Division Groundfish Assessment Program (GAP) of the Alaska Fisheries Science Center (AFSC) conducts fisheries-independent bottom trawl surveys to assess the populations of demersal fish and crab stocks of Alaska. 

'noaa_afsc_public_foss'")

write.table(str0, 
            file = here::here("R","noaa_afsc_public_foss.R"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Load Biomass data ------------------------------------------------------------

noaa_afsc_biomass_estimates <- RODBC::sqlQuery(
  channel = channel_products, 
  query = 
    "SELECT DISTINCT 
bb.SURVEY_DEFINITION_ID,
bb.SPECIES_CODE,
bb.YEAR,
bb.BIOMASS_MT,
bb.BIOMASS_VAR,
bb.POPULATION_COUNT,
bb.POPULATION_VAR
FROM GAP_PRODUCTS.AKFIN_BIOMASS bb
WHERE bb.SURVEY_DEFINITION_ID = 98 
AND bb.SPECIES_CODE IN (21740, 10210, 69322) 
AND AREA_ID = 99901
AND bb.YEAR >= 1982
") %>% 
  janitor::clean_names()

# Save table to local directory
save(noaa_afsc_biomass_estimates, file = here::here("data", "noaa_afsc_biomass_estimates.rda"))

column <- metadata_colname %>%
  dplyr::filter(metadata_colname %in% toupper(names(noaa_afsc_biomass_estimates))) %>%
  dplyr::mutate(metadata_colname = tolower(metadata_colname)) %>%
  dplyr::distinct()

str0 <- paste0("#' @title Biomass Estimates from AKFIN for EBS walleye pollock, yellowfin sole, and red king crab from 1982 to present
#' @description ",metadata_table_comment$COMMENT[metadata_table_comment$TABLE_NAME == "AKFIN_BIOMASS"]," 
#' @usage data('noaa_afsc_biomass_estimates')
#' @author Emily Markowitz (Emily.Markowitz AT noaa.gov)
#' @format A data frame with ",nrow(noaa_afsc_biomass_estimates)," observations on the following ",
               ncol(noaa_afsc_biomass_estimates)," variables.
#' \\describe{
",
               paste0(paste0("#'   \\item{\\code{",column$metadata_colname,"}}{", column$metadata_colname_long, ". ", column$metadata_colname_desc,"}"), collapse = "\n"),
               "#'   }
#' @source https://github.com/afsc-gap-products/gap_products
#' @keywords species code data biomass
#' @examples
#' data(noaa_afsc_biomass_estimates)
#' @details The Resource Assessment and Conservation Engineering (RACE) Division Groundfish Assessment Program (GAP) of the Alaska Fisheries Science Center (AFSC) conducts fisheries-independent bottom trawl surveys to assess the populations of demersal fish and crab stocks of Alaska. 

'noaa_afsc_biomass_estimates'")

write.table(str0, 
            file = here::here("R","noaa_afsc_biomass_estimates.R"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Document and create Package --------------------------------------------------

# library(remotes)
# remotes::install_github("DTUAqua/DATRAS/DATRAS")

.rs.restartR()

# options(rmarkdown.html_vignette.check_title = FALSE)
Sys.setenv('PATH' = paste0('C:/Program Files/qpdf-10.3.1/bin;', Sys.getenv('PATH')))
library(here)
library(devtools)
library(usethis)
library(roxygen2)
library(RODBC)
devtools::document()
setwd("..")
install("sdmgamindex")
3
setwd(here::here())
# devtools::check()

## Create Documentation GitHub-Pages -------------------------------------------

.rs.restartR()
# devtools::install_github("rstudio/fontawesome", force = T)
# library(fontawesome)
library(here)
library(usethis)
library(pkgdown)
rmarkdown::render(here::here("inst", "misc", "README.Rmd"),
                  output_dir = "./",
                  output_file = "README.md")


# devtools::install_github("r-lib/pkgdown")
# pkgdown::build_favicons()
# devtools::build_vignettes()
# usethis::use_pkgdown(config_file = "./pkgdown/_pkgdown.yml")

# pkgdown::clean_site()
# pkgdown::build_site(pkg = here::here())
pkgdown::build_site()
# usethis::use_github_action("pkgdown")

