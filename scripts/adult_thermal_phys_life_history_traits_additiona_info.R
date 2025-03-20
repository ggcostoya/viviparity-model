
## Assemble adult thermal physiology and life-history traits database ## 

## Packages ----

library(tidyverse)
library(readxl)
library(elevatr)
library(sf)

## Adult thermal physiology data ----

# Buckley et al. 2022 data is available at: 
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.vhhmgqnwq
# Some species names where modified in the original file to ease processing

# read Buckley et al. TPC data
tpc <- read.csv("raw_data/buckley_et_al_2022/tpcs.csv")

# select only lizard data with complete thermal physiology information
tpc <- tpc |> 
  filter(taxa == "lizards") |> 
  rename(ctmin_a = CTmin, ctmax_a = CTmax, topt_a = Topt) |> 
  filter(!is.na(ctmin_a) & !is.na(ctmax_a) & !is.na(topt_a))

# assess elevation for each species based on coordinates
elev <- st_as_sf(tpc, coords = c("lon", "lat"), crs = 4326)
tpc$elev <- get_elev_point(elev, src = "aws")$elevation
tpc$elev <- ifelse(tpc$elev < 0, 0, tpc$elev) # set elev 0 if coordinates in sea

# select columns of interest and add source column
tpc <- tpc |> select(species, lat, lon, elev, topt_a, ctmin_a, ctmax_a) |>
  mutate(source = "b")

# Dominguez-Guerrero et al. 2022 data is available at:
# https://www.nature.com/articles/s41467-022-30535-w#Sec30
# multiple columns where modified from the original file to ease processing

# read Dominguez-Guerrero et al. 2022 data
dg <- read.csv(
  "raw_data/dominguez_guerrero_et_al_2022/dominguez_guerrero_et_al_2022.csv")

# process data and select columns of interest
dg_tpc <- dg |> 
  mutate(species = paste(genus, species, sep = "_")) |> 
  mutate(ctmin_a = as.numeric(ctmin), ctmax_a = as.numeric(ctmax),
         topt_a = as.numeric(tpref), lat = as.numeric(latitude), 
         lon = as.numeric(longitude), elev = as.numeric(elevation)) |>
  filter(!is.na(ctmin_a) & !is.na(ctmax_a) & !is.na(topt_a)) |>
  mutate(source = "dg") |>
  dplyr::select(species, lat, lon, elev, topt_a, ctmin_a, ctmax_a, source) 

# Combine Buckley et al. 2022 with Dominguez-Guerrero et al. 2022 data. 
tphys_a <- rbind(tpc, dg_tpc)

## Life-history and complementary information ----

# SquamBase is available at: 

# https://datadryad.org/stash/dataset/doi:10.5061/dryad.76hdr7t3b
# the original file was modified in the following ways:
# 1. the column "sub-order" was renamed to "suborder" 
# 2. the column "mean female mass (derived from allometric equations; in log10; g)
#    was rename to "mean_female_mass_log10"
# 3. the column "mean hatchling/neonate mass (derived from allometric equations;
#    in log10: g)" was renamed to "mean_hatchling_mass_log10"
# 4. the column "minimum mean brood size" was changed to "min_n"
# 5. the column "maximum mean brood size" was changed to "max_n"
# 6. the column "Maximum Longevity (years)" was changed to "max_l"
# 7. the column "longevity in" was changed to "l_in"
# 8. the column "youngest age at first breeding (female, months)" was changed to
#    "min_m"
# 9. the column "oldest age at first breeding (female, months)" was changed to
#    "max_m"
# 10. the column "Activity time" was changed to "activity"
# 11. the column "Biome" was changed to "biome"
# 12. the column "Mode of reproduction" was changed to "parity"

# read SquamBase dataset
squambase <- read_excel("raw_data/squambase/SquamBase.xlsx")

# re-format the species column
squambase$species <- gsub(" ", "_", squambase$`Species name (Binomial)`)

# filter, process, select columns of interest and add source column
squambase_lh_traits <- squambase |> 
  filter(suborder == "Sauria") |> # filter for lizards only
  group_by(species) |>
  mutate(za = 10^as.numeric(mean_female_mass_log10), # convert from log10
         ze = 10^as.numeric(mean_hatchling_mass_log10)) |> # convert from log10
  mutate(n = mean(as.numeric(min_n), as.numeric(max_n), na.rm = T)) |> 
  select(parity, species, za, ze, n) |> 
  ungroup() |> 
  filter(!is.na(za*ze*n)) |> # filter species with complete info for Za, Ze, N
  mutate(source_lh = "sb")

# process Dominguez-Guerrero et al. 2022 data to obtain life-history info
dg_lh_traits <- dg |> 
  mutate(species = paste(genus, species, sep = "_")) |> 
  mutate(za = as.numeric(adult_mass),
         ze = as.numeric(hatchling_mass),
         n = as.numeric(negg)) |> 
  mutate(l = NA, m = NA, activity = NA, substrate = NA, biome = NA) |> 
  select(parity, species, za, ze, n) |>
  filter(!is.na(za*ze*n)) |>
  mutate(source_lh = "dg") # add source column

# combine life history trait data 
lh_traits <- rbind(squambase_lh_traits, dg_lh_traits) 

# filter data to prioritize "dg" sources of Za, Ze and N over SquamBase
lh_traits <- lh_traits |> 
  mutate(priority = ifelse(source_lh == "dg", 1, 2)) |> 
  group_by(species) |> 
  filter(priority == min(priority)) |> 
  ungroup() |> 
  select(-priority)

## Combine data and add additional information ----

# combine adult thermal physiology with life history traits data 
database <- merge(tphys_a, lh_traits, by = c("species"))

# extract additional information (L, M, Activity, Substrate, Biome) from SquamBase
squambase_add <- squambase |> 
  filter(species %in% database$species) |>
  group_by(species) |>
  mutate(l = ifelse(grepl("ature", l_in), as.numeric(max_l), NA)) |>
  mutate(m = mean(as.numeric(min_m), as.numeric(max_m), na.rm = T)/12) |> # transformed to years
  select(species,l, m, activity, substrate, biome) |> 
  ungroup()

# combine database with additional information
database <- merge(database, squambase_add, by = c("species"), all = TRUE)

# calculate expected gamma from L and M
database$gamma_ex <- rep(NA, nrow(database))
for(i in 1:nrow(database)){
  
  # calculate range of the sum
  range <- round(database$l)[i] - ceiling(database$m[i])
  
  # if there is data
  if(!is.na(range)){
    
    # initial sum value at 0
    sum_value <- 0
    
    # loop to calculate sum
    for(k in 1:(range - 1)){
      sum_value  <- sum_value + (range - k) * 0.01^(k/range)}
    
    # define gamma
    database$gamma_ex[i] <- (database$n[i]/ (range))*sum_value
    
    # if gamma is NaN set to 0
    database$gamma_ex[i] <- ifelse(is.nan(database$gamma_ex[i]), 0,
                                   database$gamma_ex[i])
    
  }else{
    database$gamma_ex[i] <- NA
  }
  
}
  
# read embryonic development and nest choice data
embryo_dev_nest_cond <- read.csv("data/embryo_development_nest_conditions.csv")

# filter relevant data
embryo_dev_nest_cond <- embryo_dev_nest_cond |> 
  select(parity, species, lat, lon, elev, eggs_seen, hatchlings_seen, 
         nest_depth, nest_shade)
  
# combine with database
database <- merge(database, embryo_dev_nest_cond,
                  by = c("parity", "species", "lat", "lon", "elev"), all = TRUE)
  
## Save data ----

# rename for easy understanding
adult_data <- database

# save csv
write.csv(adult_data, "data/adult_data.csv")
  
  
  
  
  
  
  

