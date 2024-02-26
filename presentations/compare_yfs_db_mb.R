library(gapindex)
library(googlesheets4)
library(dplyr)

# Get db index for YFS
sql_channel <- gapindex::get_connected()
dat <- gapindex::get_data(
  year_set = 1995:2023,
  survey_set = c("EBS"),
  spp_codes = c(
    10210 # yfs
  ),
  abundance_haul = "Y", sql_channel = sql_channel
)

cpue <- gapindex::calc_cpue(racebase_tables = dat)
## Calculate stratum-level biomass, population abundance, mean CPUE and
## associated variances
biomass_stratum <- gapindex::calc_biomass_stratum(
  racebase_tables = dat,
  cpue = cpue
)

## Calculate aggregated biomass and population abundance across subareas,
## management areas, and regions
biomass_subareas <- gapindex::calc_biomass_subarea(
  racebase_tables = dat,
  biomass_strata = biomass_stratum
)

db_yfs <- biomass_subareas |>
  filter(AREA_ID == 99901) |> # EBS
  mutate(BIOMASS_CV_mcs = sqrt(BIOMASS_VAR) / BIOMASS_MT, index_type = "db")

# Compare YFS design- and model-based
# Get mb index from google drive

mb_yfs <- googlesheets4::read_sheet(ss = "https://docs.google.com/spreadsheets/d/1N3Y0kLsRthR34zdOCvzUb5w3rBF5kArn6YBqlZRxkOM/edit#gid=1695026758")
mb_yfs <- mb_yfs |>
  mutate_at(
    .vars = c("Estimate", "Std. Error for Estimate"),
    .funs = list(mt = function(x) x / 1000)
  ) |>
  tibble::add_column(index_type = "VAST") |>
  mutate(BIOMASS_CV_mcs = `Std. Error for Estimate_mt` / Estimate_mt) |>
  rename(
    "YEAR" = "Time",
    "BIOMASS_MT" = "Estimate_mt"
  ) |>
  dplyr::filter(Stratum == "EBS")

plot(mb_yfs$YEAR, mb_yfs$BIOMASS_MT, type = "l", xlab = "Year", ylab = "Index (mt)", main = "Yellowfin sole")
lines(db_yfs$YEAR, db_yfs$BIOMASS_MT, col = "red")
legend(legend = c("VAST", "design-based"), "topright", lty = c(1, 1), col = c("black", "red"))
