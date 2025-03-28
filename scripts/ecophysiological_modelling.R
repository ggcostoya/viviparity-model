
## Ecophysiological modelling ## 

## Packages ----

library(tidyverse)
library(NicheMapR)

## Load data ----

# load data on adult thermal physiology and LH traits and additional info
database <- read.csv("data/adult_data.csv")

# load embryonic thermal physiology data
tphys_e <- read.csv("data/embryonic_thermal_phys.csv")

# define embryonic thermal physiology traits
topt_e <- mean(tphys_e$topt_e, na.rm = T)
ctmin_e <- mean(tphys_e$ctmin_e, na.rm = T)
ctmax_e <- mean(tphys_e$ctmax_e, na.rm = T)

## Run ecophysiological modelling ----

# holder data set
eco_data <- tibble(species = c(), lat = c(), lon = c(), elev = c(),
                   month = c(), tg = c(), tap = c(), tep = c(), 
                   shade = c(), depth = c(), tenv = c())

# loop to predict temperatures
for(i in 1:nrow(database)){
  
  # slightly modify coordinates in populations too close to sea
  add <- ifelse(paste(database$species[i], database$lat[i], sep = "_") %in% 
                  c("Hemiergis_peronii_-34.5", "Uta_stansburiana_30.436"),
                0.25, 0)
  
  # run microhabitat model
  micro <- micro_global(loc = c(database$lon[i] + add, database$lat[i]),
                        elev = database$elev[i],
                        nyears = 10,
                        timeinterval = 365,
                        run.gads = 2,
                        Refhyt = 2)
  
  # determine if diurnal, nocturnal or cathemeral
  diurnal <- ifelse(database$activity[i] != "Nocturnal", 1, 0)
  nocturnal <- ifelse(database$activity[i] != "Diurnal", 1, 0)
  
  # determine if it can climb for thermoregulation
  climb <- ifelse(database$substrate[i] %in% 
                    c("Terrestrial", "Fossorial&Terrestrial", "Cryptic"), 0, 1)
  
  # run ecotherm model during gestation
  gest <- ectotherm(Ww_g = database$za[i],
                    T_pref = topt_e,
                    CT_max = ctmax_e,
                    CT_min = ctmin_e,
                    shape = 3, # lizard shaped
                    shade_seek = 1, # can seek the shade
                    climb = climb, # can climb to thermoregulate
                    diurn = diurnal, # has diurnal activity
                    nocturn = nocturnal, # also has nocturnal activity
                    mindepth = 1, # surface
                    maxdepth = 4) # 10 cm below ground
  
  # process data from ectotherm model output during gestation
  gest <- as.data.frame(gest$environ) |> 
    mutate(month = ceiling(DOY/30.417)) |> 
    group_by(month) |>
    summarise(tg = mean(TC)) 
  
  # run ectotherm model post gestation
  post <- ectotherm(Ww_g = database$za[i], 
                    T_pref = database$topt_a[i], # topt of the adult
                    CT_max = database$ctmax_a[i], # ctmax of the adult
                    CT_min = database$ctmin_a[i], # ctmin of the adult
                    shape = 3, # lizard shaped
                    shade_seek = 1, # can seek the shade
                    climb = climb, # can climb to thermoregulate
                    diurn = diurnal, # has diurnal activity
                    nocturn = nocturnal, # also has nocturnal activity
                    mindepth = 1, # surface
                    maxdepth = 4) # 10 cm below ground
  
  # extract data from ectotherm model post gestation
  post <- as.data.frame(post$environ) |> 
    mutate(month = ceiling(DOY/30.417)) |> 
    group_by(month) |> 
    summarise(tap = mean(TC)) 
  
  # get temperature experienced by egg in the sun
  soil_temp_sun <- as.data.frame(micro$soil) |> 
    mutate(month = ceiling(DOY/30.417)) |>
    group_by(month) |> 
    summarise(d0 = mean(D0cm), d5 = mean(D5cm), 
              d10 = mean(D10cm), d15 = mean(D15cm)) |> 
    pivot_longer(cols = -month, names_to = "depth", values_to = "tep") |>
    mutate(depth = as.numeric(gsub("d", "", depth))) |>
    mutate(shade = 0)
  
  # get temperature experienced by egg in shade
  soil_temp_shade <- as.data.frame(micro$shadsoil) |> 
    mutate(month = ceiling(DOY/30.417)) |>
    group_by(month) |> 
    summarise(d0 = mean(D0cm), d5 = mean(D5cm), 
              d10 = mean(D10cm), d15 = mean(D15cm)) |> 
    pivot_longer(cols = -month, names_to = "depth", values_to = "tep") |>
    mutate(depth = as.numeric(gsub("d", "", depth))) |>
    mutate(shade = 1)
  
  # get temperature experienced by embryo average shade
  soil_temp_half <- tibble(month = soil_temp_sun$month,
                           depth = soil_temp_sun$depth,
                           tep = (soil_temp_sun$tep + soil_temp_shade$tep)/2,
                           shade = 0.5)
  
  # merge soil data
  soil <- rbind(soil_temp_sun, soil_temp_half, soil_temp_shade)
  
  # get average environmental temperature
  tenv <- as.data.frame(micro$shadmet) |> 
    mutate(month = ceiling(DOY/30.417)) |> 
    group_by(month) |>
    summarise(tenv = mean(TAREF)) # get air temperature at 2m height in shade
  
  # merge data
  temp_list <- list(gest, post, soil, tenv)
  sps_eco_data <- purrr::reduce(.x = temp_list, merge, by = c("month"), all = T)
  
  # add species and coordinates
  sps_eco_data <- sps_eco_data |> 
    mutate(species = database$species[i],
           lat = database$lat[i],
           lon = database$lon[i],
           elev = database$elev[i]) |> 
    select(species, lat, lon, elev, month, tg, tap, tep, shade, depth, tenv)
  
  # bind to holder data set
  eco_data <- rbind(eco_data, sps_eco_data)
  
  print(paste(i, "out of ", nrow(database), " populations processed."))
  
}

## Save data ----

# write a .csv
write.csv(eco_data, "data/eco_data.csv")


