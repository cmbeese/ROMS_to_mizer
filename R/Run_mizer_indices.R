# Environmental indices for the Gulf of Alaska mizer model (NMFS areas 610-650,
# 1000 m isobath). Follows "ROMS index generation example.qmd" exactly: same
# input files, same delta_correction() call and arguments (ref_yrs = 1990:2014,
# lognormal = FALSE, use_sd = TRUE), only the variables and areas differ.
#
# Temperature comes from the depth/space-averaged files (nep_avg_*), the
# plankton variables from the depth-integrated files (nep_sum_*, depthclass
# "All", mg C m^-2), as in the example's copepod index.
#
# Scenarios follow the Topic 2 run specification (GFDL earth system model):
#   ssp126, ssp245, ssp585   delta-corrected GFDL-ESM2M projections
#   stable_climate           "persistence": the 1991-2020 hindcast monthly mean
#                            repeated for every projection year
# The hindcast is spliced in up to 2020 and projections start in 2021.
#
# Run from the repo root:  Rscript R/Run_mizer_indices.R

library(dplyr)
library(lubridate)
library(tidyr)
source("R/Delta_correction.R")

isobath <- "1000"
projection_start <- 2021
climate_reference_years <- 1991:2020   # "historical climate years" in the run spec
dir_in  <- "Data/NEP_10k_revised_indices"
dir_out <- "Output"

temp_vars <- "temp"
zoo_vars  <- c("Cop", "NCa", "Eup", "MZS", "MZL")  # small/large copepods, euphausiids, micro-zooplankton
phyto_vars <- c("PhS", "PhL")                       # small/large phytoplankton
sum_vars  <- c(zoo_vars, phyto_vars)

load_sims <- function(kind, vars) {
  files <- c(hindcast   = sprintf("nep_%s_hind_%s.csv", kind, isobath),
             historical = sprintf("nep_%s_wb_hist_%s.csv", kind, isobath),
             ssp126     = sprintf("nep_%s_wb_ssp126_%s.csv", kind, isobath),
             ssp245     = sprintf("nep_%s_wb_ssp245_%s.csv", kind, isobath),
             ssp585     = sprintf("nep_%s_wb_ssp585_%s.csv", kind, isobath))
  out <- lapply(names(files), function(s) {
    d <- read.csv(file.path(dir_in, files[[s]]))
    d <- d[d$varname %in% vars, ]
    d$simulation <- s
    d
  })
  names(out) <- names(files)
  lapply(out, function(d) {
    d %>% mutate(date = lubridate::as_date(date),
                 month = lubridate::month(date),
                 year = lubridate::year(date))
  })
}

correct <- function(sims, scenario, include_hindcast) {
  delta_correction(
    hindcast   = sims$hindcast,
    historical = sims$historical,
    projection = sims[[scenario]],
    ref_yrs = 1990:2014,
    lognormal = FALSE,
    use_sd = TRUE,
    include_hindcast = include_hindcast) %>%
    mutate(simulation = scenario, hindcast_spliced = include_hindcast)
}

# Persistence scenario: every projection year repeats the hindcast's monthly
# mean over the climate reference years, separately for each area, depth class
# and variable.
stable_climate <- function(hindcast, last_year = 2099) {
  monthly_mean <- hindcast %>%
    filter(year %in% climate_reference_years) %>%
    group_by(NMFS_AREA, depthclass, varname, unit, month) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  projection_years <- expand.grid(year = projection_start:last_year, month = 1:12)
  projection <- merge(monthly_mean, projection_years, by = "month") %>%
    mutate(date = as.Date(sprintf("%d-%02d-15", year, month)),
           value_dc = value,
           simulation = "stable_climate", hindcast_spliced = TRUE) %>%
    select(NMFS_AREA, depthclass, varname, year, month, date, value, value_dc,
           unit, simulation, hindcast_spliced)

  observed <- hindcast %>%
    filter(year < projection_start) %>%
    mutate(value_dc = value, simulation = "stable_climate", hindcast_spliced = TRUE) %>%
    select(names(projection))
  bind_rows(observed, projection)
}

run_group <- function(kind, vars) {
  sims <- load_sims(kind, vars)
  scenarios <- c("ssp126", "ssp245", "ssp585")
  res <- bind_rows(lapply(scenarios, function(sc) {
    bind_rows(correct(sims, sc, FALSE), correct(sims, sc, TRUE))
  }))
  hind <- sims$hindcast %>%
    mutate(value_dc = value, simulation = "hindcast", hindcast_spliced = NA) %>%
    select(NMFS_AREA, depthclass, varname, year, month, date, value, value_dc,
           unit, simulation, hindcast_spliced)
  bind_rows(res, hind, stable_climate(sims$hindcast))
}

dir.create(dir_out, showWarnings = FALSE)
temperature <- run_group("avg", temp_vars)
plankton    <- run_group("sum", sum_vars) %>% filter(depthclass == "All")

write.csv(temperature, file.path(dir_out, "GOA_mizer_temperature_610_650_1000m.csv"), row.names = FALSE)
write.csv(plankton,    file.path(dir_out, "GOA_mizer_plankton_610_650_1000m.csv"),    row.names = FALSE)
cat("temperature rows:", nrow(temperature), " plankton rows:", nrow(plankton), "\n")
