
## Assemble embryonic thermal physiology database ##

## Packages ----

library(tidyverse)
library(readxl)
library(elevatr)
library(sf)

# files for Pettersen et al. 2023 can be found at: https://osf.io/jt28v/
# files for the reptile development database can be found at: https://repdevo.com/database/

## Data reading processing ----

# read Topt_e data from Pettersen et al. 2023
topt_e <- read.csv("raw_data/pettersen_et_al_2023/Topt_PBT.csv")

# extract Topt_e data from Pettersen et al. 2023
topt_e <- topt_e |> 
  filter(!is.na(Topt)) |> 
  filter(Group == "Lizard") |> # filter for lizards specifically
  select(Rep_mode, Species, Topt) |> 
  rename(parity = Rep_mode, species = Species, topt_e = Topt)

# read TPC data from Pettersen et al. 2023
tpc <- read.csv("raw_data/pettersen_et_al_2023/TPCs_Toptdata.csv")

# obtain CTmin_e and CTmax_e data TPC data
tpc_e <- tpc |> 
  group_by(Rep_mode, Species) |> 
  mutate(topt_e = mean(T[survival == max(survival, na.rm = T)])) |> 
  mutate(ctmin_e = ifelse(length(T[survival < 0.1 & T < topt_e]) == 0,
                          NA, max(T[survival < 0.1 & T < topt_e]))) |> 
  mutate(ctmax_e = ifelse(length(T[survival < 0.1 & T > topt_e]) == 0,
                          NA, min(T[survival < 0.1 & T > topt_e]))) |>
  select(Rep_mode, Species, ctmin_e, ctmax_e) |> 
  filter(!is.na(mean(c(ctmin_e, ctmax_e), na.rm = T))) |>
  ungroup() |> 
  unique() |>
  rename(parity = Rep_mode, species = Species)

# read reptile development database 
rdd <- read.csv("raw_data/reptile_development_database/Database.csv")

# obtain CTmin_e and CTmax_e from TPC data
rdd <- rdd |> 
  filter(order == "Squamata") |> 
  filter(const_fluct == "Const") |> # constant temp treatment
  filter(egg_embryo_hatchling == "Hatchling") |> # hatchling stage
  filter(trait == "Hatching success") |> # hatching success specifically
  mutate(species = paste(genus_timetree, species_timetree, sep = "_")) |>
  group_by(species) |> 
  mutate(hatch_success = as.numeric(mean)/max(as.numeric(mean))) |>
  mutate(topt_e = mean(T[hatch_success == max(hatch_success, na.rm = T)])) |> 
  mutate(ctmin_e = ifelse(length(T[hatch_success < 0.1 & T < topt_e]) == 0,
                          NA, max(T[hatch_success < 0.1 & T < topt_e]))) |> 
  mutate(ctmax_e = ifelse(length(T[hatch_success < 0.1 & T > topt_e]) == 0,
                          NA, min(T[hatch_success < 0.1 & T > topt_e]))) |>
  select(species, ctmin_e, ctmax_e) |> 
  filter(!is.na(mean(c(ctmin_e, ctmax_e), na.rm = T))) |>
  ungroup() |> 
  unique() 

# merge CTmin_e and CTmax_e data sources
tphys_e <- merge(tpc_e, rdd, by = c("species", "ctmin_e", "ctmax_e"), all = T)

# process merged CTmin_e and CTmax_e data
tphys_e <- tphys_e |> 
  group_by(parity, species) |> 
  mutate(ctmin_e = mean(ctmin_e, na.rm = T)) |>
  mutate(ctmax_e = mean(ctmax_e, na.rm = T)) |> 
  ungroup()

# merge CTmin_e and CTmax_e data with Topt_e data
tphys_e <- merge(tphys_e, topt_e, by = c("species"), all = T)

# process merging of data ources about CTmin_e, CTmax_e and Topt_e
tphys_e <- tphys_e |> 
  mutate(parity = ifelse(is.na(parity.x), parity.y, parity.x)) |>
  group_by(parity, species) |> 
  mutate(ctmin_e = mean(ctmin_e, na.rm = T)) |> 
  mutate(ctmax_e = mean(ctmax_e, na.rm = T)) |> 
  select(parity, species, ctmin_e, topt_e, ctmax_e) |>
  ungroup()

# manually remove snakes and tuatara, and add parity mode to missing species
tphys_e <- tphys_e |> 
  filter(!species %in% c("Naja_atra", "Natrix_natrix", "Python_molurus",
                         "Sphenodon_punctuatus", "Thamnophis_elegans")) |> 
  mutate(parity = ifelse(species == "Scincella_modesta", "Oviparous", parity))

## Save data ----

# rename for easy understanding
embryonic_thermal_phys <- tphys_e

# save csv
write.csv(embryonic_thermal_phys, "data/embryonic_thermal_phys.csv")

