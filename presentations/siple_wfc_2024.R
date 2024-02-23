# Homespun GAMs with CIs

# Load libraries ----------------------------------------------------------
PKG <- c(
  "sdmTMB", # install.packages("sdmTMB", dependencies = TRUE)
  "mgcv",
  "gratia",
  "visreg",
  "gstat",
  "dplyr",
  "ggplot2",
  "INLA",
  "prediction",
  "inlabru",
  "purrr",
  "dplyr"
)

for (p in PKG) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    require(p, character.only = TRUE)
  }
}

options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")
theme_set(theme_light())


# Load map backgrounds ----------------------------------------------------

reg_dat_ebs <- akgfmaps::get_base_layers(
  select.region = "ebs",
  set.crs = "EPSG:3338"
)
reg_dat_ebs$survey.area <- reg_dat_ebs$survey.area |>
  dplyr::mutate(
    SRVY = "EBS",
    color = scales::alpha(colour = "grey80", 0.7),
    SURVEY = "Eastern Bering Sea"
  )


# Map of strata for EBS ---------------------------------------------------
png("output/EBS_with_strata.png", width = 5, height = 5, units = "in", res = 300)
ggplot() +
  geom_sf(data = reg_dat_ebs$akland) +
  geom_sf(data = reg_dat_ebs$bathymetry) +
  geom_sf(data = reg_dat_ebs$survey.area, fill = NA) +
  geom_sf(data = reg_dat_ebs$graticule, color = "grey70", alpha = 0.5) +
  geom_sf(data = reg_dat_ebs$survey.strata, aes(fill = Stratum)) +
  coord_sf(
    xlim = reg_dat_ebs$plot.boundary$x,
    ylim = reg_dat_ebs$plot.boundary$y
  ) +
  scale_fill_manual(values = MetBrewer::met.brewer(name = "Renoir", n = length(unique(reg_dat_ebs$survey.strata$Stratum)))) +
  scale_x_continuous(
    name = "Longitude",
    breaks = reg_dat_ebs$lon.breaks
  ) +
  scale_y_continuous(
    name = "Latitude",
    breaks = reg_dat_ebs$lat.breaks
  )
dev.off()

# Load CPUE data ----------------------------------------------------------
load(file = "data/cpue_allspps.rds") # object: cpue_tab

# Plot to check it out
png("presentations/figs/threespps_rawcpue.png", width = 12, height = 5, units = "in", res = 200)
cpue_tab |>
  filter(YEAR == 2023 & CPUE_KGKM2==0) |>
  ggplot() +
  geom_point(mapping = aes(
    x = LONGITUDE_DD_START, y = LATITUDE_DD_START,
    size = CPUE_KGKM2, color = BOTTOM_TEMPERATURE_C
  ),
    color = "red",
    shape = 4,
    size = 3) +
  geom_point(
    data = filter(cpue_tab, YEAR == 2023 & CPUE_KGKM2 > 0),
    aes(x = LONGITUDE_DD_START, y = LATITUDE_DD_START,
        size = CPUE_KGKM2,
        color = BOTTOM_TEMPERATURE_C),
    alpha = 0.6
  ) +
  scale_color_viridis_c(expression(Temperature~(degree*C))) +
  scale_size_area(max_size = 15,expression(CPUE~(kg/km^2))) +
  # geom_sf(
  #   data = reg_dat_ebs$akland,
  #   color = NA,
  #   fill = "grey40"
  # ) +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~species_name) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
dev.off()


# Load and format grid ----------------------------------------------------

# This is another form of the same grid:
load(here::here("data/pred_grid_ebs.rda")) # object: pred_grid_ebs

get_crs(dat = pred_grid_ebs, ll_names = c("lon", "lat"))

grid <- add_utm_columns(pred_grid_ebs,
  # ll_crs = 32603,
  ll_names = c("lon", "lat")
)

grid <- purrr::map_dfr(unique(cpue_tab$YEAR), ~ tibble(grid, year = .x))
range(grid$X)

# Fit basic GAM with not that many effects --------------------------------

start.time <- Sys.time()

fit_gam_10120 <- gam(
  formula = CPUE_KGKM2 ~ as.factor(YEAR) +
    s(X, Y, k = 50),
  family = tw(link = "log"),
  data = dplyr::filter(cpue_tab, SPECIES_CODE == 10120)
)

cat(
  "The GAM took ",
  difftime(Sys.time(), start.time, units = "mins"),
  " mins to run"
)

fit_gam_21740 <- gam(
  formula = CPUE_KGKM2 ~ as.factor(YEAR) +
    s(X, Y, k = 50),
  family = tw(link = "log"),
  data = dplyr::filter(cpue_tab, SPECIES_CODE == 21740)
)

fit_gam_69322 <- gam(
  formula = CPUE_KGKM2 ~ as.factor(YEAR) + # temporal
    s(X, Y, k = 50), # spatial
  family = tw(link = "log"),
  data = dplyr::filter(cpue_tab, SPECIES_CODE == 69322)
)


png("QQ_all_spatial_temporal.png", width = 8, height = 3.1, units = "in", res = 200)
par(mfrow = c(1, 3))
qq.gam(fit_gam_10120, pch = 20)
qq.gam(fit_gam_21740, pch = 20)
qq.gam(fit_gam_69322, pch = 20)
dev.off()


# Fit a GAM analogous to a regular spatiotemporal model -------------------
# This one has no covariates!

start.time <- Sys.time()

fit_gam_s_t_st_10120 <- gam(
  formula = CPUE_KGKM2 ~ as.factor(YEAR) + # temporal
    s(X, Y, bs = c("ts", k = 375)) + # spatial
    s(X, Y, bs = c("ts"), k = 50, by = as.factor(YEAR), id = 1), # spatiotemporal
  family = tw(link = "log"),
  data = dplyr::filter(cpue_tab, SPECIES_CODE == 10120)
)

cat(
  "The GAM took ",
  difftime(Sys.time(), start.time, units = "mins"),
  "mins to run"
)

start.time <- Sys.time()

fit_gam_s_t_st_21740 <- gam(
  formula = CPUE_KGKM2 ~ as.factor(YEAR) + # temporal
    s(X, Y, bs = c("ts", k = 375)) + # spatial
    s(X, Y, bs = c("ts"), k = 50, by = as.factor(YEAR), id = 1), # spatiotemporal
  family = tw(link = "log"),
  data = dplyr::filter(cpue_tab, SPECIES_CODE == 21740)
)

cat(
  "The GAM took ",
  difftime(Sys.time(), start.time, units = "mins"),
  "mins to run"
)

start.time <- Sys.time()

fit_gam_s_t_st_69322 <- gam(
  formula = CPUE_KGKM2 ~ as.factor(YEAR) + # temporal
    s(X, Y, bs = c("ts", k = 375)) + # spatial
    s(X, Y, bs = c("ts"), k = 50, by = as.factor(YEAR), id = 1), # spatiotemporal
  family = tw(link = "log"),
  data = dplyr::filter(cpue_tab, SPECIES_CODE == 69322)
)

cat(
  "The GAM took ",
  difftime(Sys.time(), start.time, units = "mins"),
  "mins to run"
)

png("QQ_all_spatial_temporal_st.png", width = 8, height = 3.1, units = "in", res = 200)
par(mfrow = c(1, 3))
qq.gam(fit_gam_s_t_st_10120, pch = 20)
qq.gam(fit_gam_s_t_st_21740, pch = 20)
qq.gam(fit_gam_s_t_st_69322, pch = 20)
dev.off()

# Formula from Casper that fits better
# fm = "Year + s(sx,sy,bs=c('ts'),k=376) +
# s(sx,sy,bs=c('ts'),k=50,by=Year,id=1) +
# s(BOTTOM_DEPTH,bs='ts',k=10) +
# s(log(GEAR_TEMPERATURE+3),bs='ts',k=10)  + offset(log(EFFORT))"


# Take GAM results and turn them into an index ----------------------------

# don't forget that the effort units are now in km^2
pred_gam <- predict(fit_gam, type = "response", newdata = grid) # This takes a long time.
pred_gam_df <- cbind(grid, pred_gam)
pred_gam_df$Shape_Area_km2 <- pred_gam_df$Shape_Area * 1e-6 #* 0.0001 for hectares; original Shape_area is in m^2
pred_gam_df$predicted_tot_grid <- pred_gam_df$Shape_Area_km2 * pred_gam_df$pred_gam

gam_idx_mt <- pred_gam_df |>
  dplyr::group_by(year) |>
  summarize(total_wt_mt = sum(predicted_tot_grid) / 1000) # convert kg --> mt



# FIGS FOR PRES -----------------------------------------------------------
load("output/gamresults_2023.RDS") # results.df
# These are the results from Emily's GAMs using Casper's pkg

results.df$index_type <- "GAM"
results.df <- dplyr::rename(results.df, idx_cv = "value")
results.df2 <- results.df |>
  tidyr::pivot_longer(
    cols = `arrowtooth flounder`:`red king crab`,
    names_to = "species", values_to = "value"
  )

# Load the VAST indices and the design-based
rkc <- read.csv("data/indices/EBS/bbrkc_2022_comparison.csv", header = TRUE) # this one was updated by Jon R after Thorson saw the comparison at WKUSER
rkc$species <- "red king crab"
yfs <- read.csv("data/indices/EBS/YFS_10210_estimate_summary.csv")
yfs$species <- "yellowfin sole"
wep <- read.csv("data/indices/EBS/WEP_21740_estimate_summary.csv")
wep$species <- "walleye pollock"

all_spps <- dplyr::bind_rows(yfs, wep, rkc)

png("output/VAST_vs_design.png", width = 8, height = 3, units = "in", res = 200)
all_spps %>%
  mutate_at(
    .vars = c("design_mt", "VAST_mt", "design_se", "VAST_se"),
    .funs = function(x) x / 1e6
  ) |>
  ggplot(aes(x = design_mt, y = VAST_mt, color = Year)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "A") +
  facet_wrap(~species, scales = "free") +
  theme_light(base_size = 14) +
  xlab("Design-based index (millions mt)") +
  ylab("Model-based index (millions mt)") +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = c(0, 0.4)),
      scale_x_continuous(limits = c(1, 8.5)),
      scale_x_continuous(limits = c(1, 3.75))
    ),
    y = list(
      scale_y_continuous(limits = c(0, 0.4)),
      scale_y_continuous(limits = c(1, 8.5)),
      scale_y_continuous(limits = c(1, 3.75))
    )
  ) +
  geom_abline(slope = 1, color = "red", linetype = "dashed")
dev.off()

png("output/VAST_vs_design_CV.png", width = 8, height = 3, units = "in", res = 200)
all_spps %>%
  mutate_at(
    .vars = c("design_mt", "VAST_mt", "design_se", "VAST_se"),
    .funs = function(x) x / 1e6
  ) |>
  ggplot(aes(x = design_CV, y = VAST_CV, color = Year)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "A") +
  facet_wrap(~species, scales = "free") +
  theme_light(base_size = 14) +
  xlab("CV of design-based index") +
  ylab("CV of model-based index") +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = c(0, 0.7)),
      scale_x_continuous(limits = c(0, 0.26)),
      scale_x_continuous(limits = c(0.04, 0.13))
    ),
    y = list(
      scale_y_continuous(limits = c(0, 0.7)),
      scale_y_continuous(limits = c(0, 0.26)),
      scale_y_continuous(limits = c(0.04, 0.13))
    )
  ) +
  geom_abline(slope = 1, color = "red", linetype = "dashed")
dev.off()


# Plot time series of each: VAST, GAM, design-based -----------------------
# all_spps2 <- all_spps |> mutate (GAM=NA,GAM_cv=NA)
all.results <- results.df2 |>
  tidyr::pivot_wider(names_from = idx_cv) |>
  rename(Year = "year", GAM = "index", GAM_cv = "cv") |>
  select(-index_type) |>
  left_join(all_spps) |>
  filter(!is.na(design_mt)) |>
  janitor::clean_names() |>
  rename(gam_mt = "gam") |>
  tidyr::pivot_longer(gam_mt:vast_cv) |>
  tidyr::separate(name, into = c("index_type", "value2")) |>
  mutate(value2 = tolower(value2)) |>
  tidyr::pivot_wider(id_cols = year:index_type, names_from = value2, values_from = value)

png("output/VAST_vs_design_ts.png", width = 8, height = 6, units = "in", res = 200)
p1 <- all.results |>
  filter(index_type != "gam") |> # species != "arrowtooth flounder" &
  mutate_at(.vars = c("mt", "se"), .funs = function(x) x / 1e6) |>
  ggplot(aes(x = year, y = mt, color = index_type, fill = index_type, group = index_type)) +
  geom_point(size = 2.5) +
  geom_line(lwd = 1.2) +
  geom_ribbon(aes(ymin = mt - se, ymax = mt + se), alpha = 0.3, color = NA) +
  facet_wrap(~species, scales = "free", ncol = 1) +
  MetBrewer::scale_color_met_d(name = "Klimt") +
  MetBrewer::scale_fill_met_d(name = "Klimt") +
  ylab("Estimated biomass (millions mt)") +
  theme_light(base_size = 14)
print(p1)
dev.off()


# Load results and stuff from vignette D ----------------------------------
# Em results - load into nested list

model_path_list <- list()
paths <- c(
  "inst/vigD_model_fits_fm_1_s_t_st.Rdata", # object:models
  "inst/vigD_model_fits_fm_2_cov.Rdata",
  "inst/vigD_model_fits_fm_3_s_t_st_cov.Rdata"
)
for (m in 1:length(paths)) {
  load(paths[[m]])
  model_path_list[[m]] <- models
}

# this one has scaling issues or something:  load("inst/VigD_model_fits_cov_vs_st.Rdata")

# Em code for tidying
YEARS <- 2015:2023
head(all_spps)

# GAM stuff
dat_combined <- data.frame()

for (m in 1:length(paths)) {
  dat <- data.frame()
  for (i in 1:length(models)) {
    temp <- model_path_list[[m]][[i]]
    dat0 <- data.frame(
      idx = temp$idx[, 1] / 1e2, # /1e2 clearly having an issue with units
      lo = temp$lo[, 1] / 1e2,
      up = temp$up[, 1] / 1e2,
      Year = rownames(temp$idx),
      group = names(models)[i],
      formula = paste0(
        "cpue_kgkm2 ~ ",
        as.character(temp$pModels[[1]]$formula)[[3]]
      )
    )
    dat <- dplyr::bind_rows(dat, dat0)
  }
  any(is.na(dat$index_type))

  dat_combined <- bind_rows(dat_combined, dat)
}

dat <- dat_combined

dat$common_name <- paste0(sapply(X = strsplit(x = dat$group, split = " fm"), `[`, 1))

dat2 <- dat |>
  dplyr::rename(species = "common_name") |>
  dplyr::select(-group) |>
  mutate(
    index_type = case_when(
      formula == 'cpue_kgkm2 ~ Year + s(sx, sy, bs = c("ts"), k = 376) + s(sx, sy, bs = c("ts"), k = 10, by = Year)' ~ "GAM_s_t_st",
      formula == 'cpue_kgkm2 ~ s(BOTTOM_DEPTH, bs = "ts", k = 10) + s(log(GEAR_TEMPERATURE + 3), bs = "ts", k = 10)' ~ "GAM_covariates",
      formula == 'cpue_kgkm2 ~ Year + s(sx, sy, bs = c("ts"), k = 376) + s(sx, sy, bs = c("ts"), k = 10, by = Year) + s(BOTTOM_DEPTH, bs = "ts", k = 10) + s(log(GEAR_TEMPERATURE + 3), bs = "ts", k = 10)' ~ "GAM_s_t_st_covariates"
    ),
    year = as.numeric(Year)
  ) |>
  dplyr::select(-formula) |>
  rename(mt = "idx") |>
  filter(year != 2020) |>
  tibble::remove_rownames() |>
  filter(up < 50e6)

png("output/VAST_vs_design_vs_GAM_ts.png", width = 8, height = 6, units = "in", res = 200)
p1 +
  geom_point(data = dat2, aes(x = year, y = mt / 1e6), size = 2.5) +
  geom_line(data = dat2, aes(x = year, y = mt / 1e6), lwd = 1.2) +
  geom_ribbon(data = dat2, aes(ymin = lo / 1e6, ymax = up / 1e6), alpha = 0.3, color = NA) +
  guides(
    color = guide_legend(title = "Index type"),
    fill = guide_legend(title = "Index type")
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 0.2)),
      scale_y_continuous(limits = c(2, 10.5)),
      scale_y_continuous(limits = c(0.5,3.5))
    )
  )
dev.off()



# AIC for all models ------------------------------------------------------
all_aic <- data.frame()
for (m in 1:length(model_path_list)) {
  AICs <- lapply(model_path_list[[m]], FUN = get_surveyidx_aic) |>
    unlist()
  modelname <- names(AICs)
  all_aic <- bind_rows(all_aic, data.frame(modelname = modelname, AIC = AICs))
}

rownames(all_aic) <- NULL
all_aic$stock <- sapply(stringr::str_split(all_aic$model, " fm_"), `[`, 1)

all_aic |>
  group_by(stock) |>
  arrange(AIC) |>
  select(stock, modelname,  AIC)


# Update pcod-crab overlap ------------------------------------------------
# For Erin Fedewa (I asked her if she had looked at an updated overlap graph)
sql_channel <- gapindex::get_connected()

dat <- gapindex::get_data(year_set = 2005:2023,
                          survey_set = c("EBS","NBS"),
                          spp_codes = c(68580, 21720),
                          abundance_haul = "Y",sql_channel = sql_channel)

cpue <- gapindex::calc_cpue(racebase_tables = dat)

write.csv(cpue,row.names = FALSE,file = "presentations/snowcrab_pcod_cpue_table.csv")
