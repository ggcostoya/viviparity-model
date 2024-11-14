
## Summaries of data sets and statistics ##

## Packages ----


## Embryonic thermal physiology summaries ----

# read embryonic thermal physiology data

# total number of species with any data
nrow(tphys_e) # 40

# number of species with any data for each parity mode
table(tphys_e$parity) # Oviparous = 39, Viviparous = 1

# average values of thermal physiology traits 
mean(tphys_e$topt_e, na.rm = T) # Topt_e = 27.57
mean(tphys_e$ctmin_e, na.rm = T) # CTmin_e = 19.41
mean(tphys_e$ctmax_e, na.rm = T) # CTmax_e = 33.83

