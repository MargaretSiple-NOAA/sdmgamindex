
# Load new FOSS data -----------------------------------------------------------

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
    dplyr::rename(TABLE_NAME = MVIEW_NAME)
) %>% 
  dplyr::filter(TABLE_NAME == "FOSS_CATCH")

metadata_table_comment <- paste0(metadata_table_comment$COMMENTS, 
                                 " Data was pulled ", Sys.Date(), ". ")

metadata_colname <- RODBC::sqlQuery(
  channel = channel_products, 
  query = "SELECT * FROM GAP_PRODUCTS.METADATA_COLUMN") %>% 
  janitor::clean_names()

noaa_afsc_public_foss <- RODBC::sqlQuery(
    channel = channel_products, 
    query = 
      "SELECT DISTINCT 
hh.YEAR,
hh.SRVY,                 
hh.SURVEY,
hh.SURVEY_DEFINITION_ID,
hh.SURVEY_NAME,
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
cc.TAXON_CONFIDENCE,
cc.WEIGHT_KG,
cc.COUNT,
cc.CPUE_KGKM2,
cc.CPUE_NOKM2,
hh.AREA_SWEPT_KM2,       
hh.DISTANCE_FISHED_KM,
hh.DURATION_HR,          
hh.NET_WIDTH_M,
hh.NET_HEIGHT_M,
hh.PERFORMANCE 
FROM GAP_PRODUCTS.FOSS_CATCH cc
LEFT JOIN GAP_PRODUCTS.FOSS_HAUL hh
ON cc.HAULJOIN = hh.HAULJOIN
FULL JOIN GAP_PRODUCTS.FOSS_SPECIES_TEST ss
ON cc.SPECIES_CODE = ss.SPECIES_CODE
WHERE hh.SRVY = 'EBS' 
AND cc.SPECIES_CODE IN (21740, 10210, 69323) 
AND hh.YEAR >= 2015") %>% 
  janitor::clean_names()

# Save table to local directory
save(noaa_afsc_public_foss, file = here::here("data", "noaa_afsc_public_foss.rda"))
  
column <- metadata_colname %>%
  dplyr::filter(metadata_colname %in% toupper(names(noaa_afsc_public_foss))) %>%
  dplyr::mutate(metadata_colname = tolower(metadata_colname)) %>%
  dplyr::distinct()

# column <- read_csv(file = "../notforgit/metadata_column_current.csv")
# column <- column[column$TABLE_NAME == "FOSS_CPUE_PRESONLY",]
# column$col_name <- names(public_data)
# table <- read_csv(file = "../notforgit/metadata_table_current.csv")

table0 <- "noaa_afsc_public_foss"

str0 <- paste0("#' @title Public data from FOSS
#' @description ",metadata_table_comment," 
#' @usage data('noaa_afsc_public_foss')
#' @author Emily Markowitz (Emily.Markowitz AT noaa.gov)
#' @format A data frame with ",nrow(noaa_afsc_public_foss)," observations on the following ",
ncol(noaa_afsc_public_foss)," variables.
#' \\describe{
",
               paste0(paste0("#'   \\item{\\code{",column$metadata_colname,"}}{", column$metadata_colname_long, ". ", column$metadata_colname_desc,"}"), collapse = "\n"),
               "#'   }
#' @source https://github.com/afsc-gap-products/gap_public_data
#' @keywords species code data
#' @examples
#' data(noaa_afsc_public_foss)
#' @details DETAILS
'noaa_afsc_public_foss'")

write.table(str0, 
            file = here::here("R","noaa_afsc_public_foss.R"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Document and create Package --------------------------------------------------
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

## Create Package Readme.md -----------------------------------------------------
rmarkdown::render(here::here("inst", "misc", "README.Rmd"),
                  output_dir = "./",
                  output_file = "README.md")

## Create Documentation GitHub-Pages -------------------------------------------

.rs.restartR()
# devtools::install_github("rstudio/fontawesome", force = T)
# library(fontawesome)
library(here)
library(usethis)
library(pkgdown)
rmarkdown::render(input = "./inst/r/README.Rmd",
                  output_dir = "./",
                  output_file = "README.md")


# devtools::install_github("r-lib/pkgdown")
# pkgdown::build_favicons()
# devtools::build_vignettes()
# usethis::use_pkgdown(config_file = "./pkgdown/_pkgdown.yml")

pkgdown::build_site(pkg = here::here())
# usethis::use_github_action("pkgdown")