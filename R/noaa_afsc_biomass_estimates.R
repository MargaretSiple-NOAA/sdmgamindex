#' @title Biomass Estimates from AKFIN for EBS walleye pollock, yellowfin sole, and red king crab from 1982 to present
#' @description This table is a copy of GAP_PRODUCTS.BIOMASS and does not have any other object dependencies. These data are produced by the Resource Assessment and Conservation Engineering Division (RACE) Groundfish Assessment Program (GAP) of the Alaska Fisheries Science Center (AFSC). There are legal restrictions on access to the data. These data are not intended for public dissemination and should not be shared without the explicit written consent of the data managers and owners (NOAA Fisheries). The GitHub repository for the scripts that created this code can be found at https://github.com/afsc-gap-products/gap_products. These data were last updated April 10, 2024. 
#' @usage data('noaa_afsc_biomass_estimates')
#' @author AFSC Groundfish Assessment Program (nmfs.afsc.gap.metadata AT noaa.gov)
#' @format A data frame with 146 observations on the following 8 variables.
#' \describe{
#'   \item{\code{area_id}}{Area ID. Area ID key code for each statistical area used to produce production estimates (e.g., biomass, population, age comps, length comps). Each area ID is unique within each survey.}
#'   \item{\code{biomass_mt}}{Estimated biomass. The estimated total biomass.}
#'   \item{\code{biomass_var}}{Estimated biomass variance. The estimated variance associated with the total biomass.}
#'   \item{\code{population_count}}{Estimated population. The estimated population caught in the survey for a species, group, or total for a given survey.}
#'   \item{\code{population_var}}{Estimated population variance. The estimated population variance caught in the survey for a species, group, or total for a given survey.}
#'   \item{\code{species_code}}{Taxon code. The species code of the organism associated with the 'common_name' and 'scientific_name' columns. For a complete species list, review the [code books](https://www.fisheries.noaa.gov/resource/document/groundfish-survey-species-code-manual-and-data-codes-manual).}
#'   \item{\code{survey_definition_id}}{Survey ID. The survey definition ID key code uniquely identifies a survey/survey design. Integer code that uniquely identifies survey. Full list of survey definition IDs are in RACE_DATA.SURVEY_DEFINITIONS. IDs used in GAP_PRODUCTS are: 47 (Gulf of Alaska); 52 (Aleutian Islands); 78 (Bering Sea Slope); 98 (Eastern Bering Sea Shelf); 143 (Northern Bering Sea Shelf). The column 'survey_definition_id' is associated with the 'srvy' and 'survey' columns. For a complete list of surveys, review the [code books](https://www.fisheries.noaa.gov/resource/document/groundfish-survey-species-code-manual-and-data-codes-manual).}
#'   \item{\code{year}}{Survey year. Year the observation (survey) was collected.}#'   }
#' @source https://github.com/afsc-gap-products/gap_products
#' @keywords species code data biomass
#' @examples
#' data(noaa_afsc_biomass_estimates)
#' @details The Resource Assessment and Conservation Engineering (RACE) Division Groundfish Assessment Program (GAP) of the Alaska Fisheries Science Center (AFSC) conducts fisheries-independent bottom trawl surveys to assess the populations of demersal fish and crab stocks of Alaska. 

'noaa_afsc_biomass_estimates'
